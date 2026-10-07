"""Evidence read from alignments: endpoints, local site metrics, and RGEN counts.

Endpoints are 0-based aligned reference bases, and clips never extend them.
`read_ends` says which ends of an alignment count:
- short reads: the 5' end, `reference_start` for forward reads and
  `reference_end - 1` for reverse reads, on the read's strand;
- long reads: both aligned ends, because a long read spans its molecule.
  For Digenome the left end is '+' and the right end '-', the two sides a cut
  leaves; for nDigenome both stay on the read's strand, the DNA strand of the
  molecule.
Counts and depths use primary alignments only; secondary and supplementary
alignments are reported as diagnostics.
"""

from __future__ import annotations

import heapq
import statistics
import struct
from collections import Counter
from dataclasses import dataclass
from typing import Callable, NamedTuple

import pysam

from .regions import MappedContig
from .settings import CallerSettings

INSERTION, DELETION, SOFT_CLIP, HARD_CLIP = 1, 2, 4, 5
REFERENCE_ADVANCING = (0, 2, 3, 7, 8)  # M, D, N, =, X


def open_bam(path: str, label: str) -> pysam.AlignmentFile:
    bam = pysam.AlignmentFile(str(path), "rb")
    if not bam.has_index():
        raise ValueError(f"{label} is not indexed: {path}")
    sort_order = bam.header.to_dict().get("HD", {}).get("SO")
    if sort_order != "coordinate":
        raise ValueError(
            f"{label} must declare coordinate sort order, found: {sort_order or 'missing'}"
        )
    return bam


def check_layout(bam: pysam.AlignmentFile, min_mapq: int, require_paired: bool, label: str) -> None:
    """Fail early on an empty BAM, or on single-end data where pairs are required."""
    checked = paired = 0
    for read in bam.fetch(until_eof=True):
        if not is_counted(read, min_mapq):
            continue
        checked += 1
        paired += int(read.is_paired)
        if checked >= 10000:
            break
    bam.reset()
    if checked == 0:
        raise ValueError(f"{label} contains no eligible primary alignments")
    if require_paired and paired == 0:
        raise ValueError(f"{label} requires paired-end alignments for nDigenome analysis")


def mapped_contigs(bam: pysam.AlignmentFile) -> list[MappedContig]:
    """Contigs with at least one mapped record, in BAM header order."""
    mapped = {stat.contig: stat.mapped for stat in bam.get_index_statistics()}
    return [
        MappedContig(name, length, mapped[name])
        for name, length in zip(bam.references, bam.lengths)
        if mapped.get(name, 0) > 0
    ]


# Single reads.


def is_counted(read: pysam.AlignedSegment, min_mapq: int) -> bool:
    """Primary, mapped, passing QC, not a duplicate, and at least `min_mapq`."""
    return not (
        read.is_unmapped
        or read.is_secondary
        or read.is_supplementary
        or read.is_qcfail
        or read.is_duplicate
        or read.mapping_quality < min_mapq
    )


def strand_of(read: pysam.AlignedSegment) -> str:
    return "-" if read.is_reverse else "+"


class ReadEnd(NamedTuple):
    position: int
    strand: str
    at_start: bool  # the end at reference_start (otherwise reference_end - 1)


def read_ends(read: pysam.AlignedSegment, settings: CallerSettings) -> list[ReadEnd]:
    """The endpoints an alignment shows (see the module docstring)."""
    if read.reference_start is None or read.reference_end is None:
        raise ValueError("Cannot calculate an endpoint for an unmapped alignment")
    left = ReadEnd(read.reference_start, strand_of(read), True)
    right = ReadEnd(read.reference_end - 1, strand_of(read), False)
    if not settings.long_reads:
        return [right] if read.is_reverse else [left]
    if settings.analysis == "digenome":
        left, right = left._replace(strand="+"), right._replace(strand="-")
    # A one-base alignment has a single end, so it must not count twice.
    return [left] if left.position == right.position else [left, right]


def end_clip_length(read: pysam.AlignedSegment, at_start: bool) -> int:
    """The soft or hard clip on one aligned end of `read`."""
    cigar = read.cigartuples or []
    if not cigar:
        return 0
    operation, length = cigar[0] if at_start else cigar[-1]
    return length if operation in (SOFT_CLIP, HARD_CLIP) else 0


def indels_overlapping(read: pysam.AlignedSegment, start: int, end: int) -> list[tuple[int, str, int]]:
    """(reference position, "INS" or "DEL", length) for each CIGAR indel that
    touches [start, end). A deletion counts when any deleted base is inside; an
    insertion when its position is inside."""
    overlapping = []
    position = read.reference_start
    for operation, length in read.cigartuples or []:
        if position >= end:
            break  # CIGAR operations run left to right, so the rest are past the window
        if operation == INSERTION:
            if position >= start:
                overlapping.append((position, "INS", length))
        elif operation == DELETION:
            if position + length > start:
                overlapping.append((position, "DEL", length))
            position += length
        elif operation in REFERENCE_ADVANCING:
            position += length
    return overlapping


# Endpoint scanning.


def count_endpoints(
    bam: pysam.AlignmentFile,
    contig: str,
    start: int,
    end: int,
    settings: CallerSettings,
    min_count: int,
) -> dict[tuple[int, str], int]:
    """Count endpoints in [start, end) on both strands and keep those with at
    least `min_count` reads.

    The BAM is streamed in coordinate order. Every endpoint lies at or after
    its read's start, so an endpoint is final once reads start past it; only
    those near the current read are held in memory.
    """
    found: dict[tuple[int, str], int] = {}
    pending: dict[tuple[int, str], int] = {}
    pending_order: list[tuple[int, str]] = []
    previous_start = -1

    def finish(key: tuple[int, str]) -> None:
        count = pending.pop(key)
        if count >= min_count:
            found[key] = count

    for read in bam.fetch(contig, start, end):
        if not is_counted(read, settings.min_mapq):
            continue
        if read.reference_start < previous_start:
            raise ValueError(
                f"BAM is not coordinate sorted on {contig}: "
                f"{read.reference_start} follows {previous_start}"
            )
        previous_start = read.reference_start
        while pending_order and pending_order[0][0] < read.reference_start:
            finish(heapq.heappop(pending_order))
        for read_end in read_ends(read, settings):
            if not start <= read_end.position < end:
                continue
            key = (read_end.position, read_end.strand)
            if key not in pending:
                pending[key] = 0
                heapq.heappush(pending_order, key)
            pending[key] += 1
    while pending_order:
        finish(heapq.heappop(pending_order))
    return found


# Local site metrics.


@dataclass(frozen=True)
class SiteMetrics:
    """Evidence at one endpoint on one strand. The defaults describe no evidence."""

    endpoint_count: int = 0
    strand_depth: int = 0
    endpoint_fraction: float = 0.0
    support_mean_mapq: float = 0.0
    support_read_count: int = 0
    local_mean_mapq: float = 0.0
    local_read_count: int = 0
    support_mean_nm: float = 0.0
    support_nm_count: int = 0
    local_mean_nm: float = 0.0
    local_nm_count: int = 0
    softclip_fraction: float = 0.0
    softclipped_count: int = 0
    secondary_endpoint_count: int = 0
    indel_position: int | None = None
    indel_type: str = ""
    indel_length: int = 0
    indel_read_count: int = 0
    indel_fraction: float = 0.0


def mean_or_zero(values: list[float]) -> float:
    return statistics.fmean(values) if values else 0.0


def ratio(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else 0.0


def measure_site(
    bam: pysam.AlignmentFile,
    contig: str,
    position: int,
    strand: str,
    settings: CallerSettings,
) -> SiteMetrics:
    """Measure endpoint support at `position` and artifact evidence within
    `cleavage_artifact_window` bases of it.

    Depth counts reads covering `position` that can give endpoints on
    `strand`; the clip metric uses the clip on the end that forms the endpoint.
    """
    min_mapq = settings.min_mapq
    start = max(0, position - settings.cleavage_artifact_window)
    end = position + settings.cleavage_artifact_window + 1
    endpoint_reads = []
    local_mapqs: list[int] = []
    local_nms: list[float] = []
    strand_depth = 0
    secondary_endpoint_count = 0
    indel_reads: set[tuple[str, int, int]] = set()
    indel_counts: Counter[tuple[int, str, int]] = Counter()

    for read in bam.fetch(contig, start, end):
        if read.is_unmapped or read.is_qcfail:
            continue
        ends = read_ends(read, settings)
        this_end = next((e for e in ends if e.position == position and e.strand == strand), None)
        if read.is_secondary or read.is_supplementary:
            if read.mapping_quality >= min_mapq and this_end is not None:
                secondary_endpoint_count += 1
            continue
        if not is_counted(read, min_mapq):
            continue
        local_mapqs.append(read.mapping_quality)
        if read.has_tag("NM"):
            local_nms.append(float(read.get_tag("NM")))
        can_end_here = any(e.strand == strand for e in ends)
        if can_end_here and read.reference_start <= position < read.reference_end:
            strand_depth += 1
        if this_end is not None:
            endpoint_reads.append((read, this_end))
        nearby = indels_overlapping(read, start, end)
        if nearby:
            indel_counts.update(nearby)
            indel_reads.add((read.query_name or "", read.flag, read.reference_start))

    support_mapqs = [read.mapping_quality for read, _end in endpoint_reads]
    support_nms = [float(read.get_tag("NM")) for read, _end in endpoint_reads if read.has_tag("NM")]
    softclipped = sum(1 for read, read_end in endpoint_reads if end_clip_length(read, read_end.at_start) > 0)
    indel_position, indel_type, indel_length = None, "", 0
    if indel_counts:
        (indel_position, indel_type, indel_length), _count = indel_counts.most_common(1)[0]

    return SiteMetrics(
        endpoint_count=len(endpoint_reads),
        strand_depth=strand_depth,
        endpoint_fraction=ratio(len(endpoint_reads), strand_depth),
        support_mean_mapq=mean_or_zero(support_mapqs),
        support_read_count=len(support_mapqs),
        local_mean_mapq=mean_or_zero(local_mapqs),
        local_read_count=len(local_mapqs),
        support_mean_nm=mean_or_zero(support_nms),
        support_nm_count=len(support_nms),
        local_mean_nm=mean_or_zero(local_nms),
        local_nm_count=len(local_nms),
        softclip_fraction=ratio(softclipped, len(endpoint_reads)),
        softclipped_count=softclipped,
        secondary_endpoint_count=secondary_endpoint_count,
        indel_position=indel_position,
        indel_type=indel_type,
        indel_length=indel_length,
        indel_read_count=len(indel_reads),
        indel_fraction=ratio(len(indel_reads), len(local_mapqs)),
    )


def weighted_mean(first: float, first_count: int, second: float, second_count: int) -> float:
    total = first_count + second_count
    if total == 0:
        return 0.0
    return (first * first_count + second * second_count) / total


def combine(forward: SiteMetrics, reverse: SiteMetrics) -> SiteMetrics:
    """Pool the two strands of a Digenome pair.

    The reported indel comes from the strand with the higher indel fraction,
    then read count; on a tie it comes from `forward`, so argument order matters.
    """
    endpoint_count = forward.endpoint_count + reverse.endpoint_count
    strand_depth = forward.strand_depth + reverse.strand_depth
    local_count = forward.local_read_count + reverse.local_read_count
    indel_count = forward.indel_read_count + reverse.indel_read_count
    softclipped = forward.softclipped_count + reverse.softclipped_count
    top_indel = max(
        (forward, reverse),
        key=lambda metrics: (metrics.indel_fraction, metrics.indel_read_count),
    )
    return SiteMetrics(
        endpoint_count=endpoint_count,
        strand_depth=strand_depth,
        endpoint_fraction=ratio(endpoint_count, strand_depth),
        support_mean_mapq=weighted_mean(
            forward.support_mean_mapq, forward.support_read_count,
            reverse.support_mean_mapq, reverse.support_read_count,
        ),
        support_read_count=forward.support_read_count + reverse.support_read_count,
        local_mean_mapq=weighted_mean(
            forward.local_mean_mapq, forward.local_read_count,
            reverse.local_mean_mapq, reverse.local_read_count,
        ),
        local_read_count=local_count,
        support_mean_nm=weighted_mean(
            forward.support_mean_nm, forward.support_nm_count,
            reverse.support_mean_nm, reverse.support_nm_count,
        ),
        support_nm_count=forward.support_nm_count + reverse.support_nm_count,
        local_mean_nm=weighted_mean(
            forward.local_mean_nm, forward.local_nm_count,
            reverse.local_mean_nm, reverse.local_nm_count,
        ),
        local_nm_count=forward.local_nm_count + reverse.local_nm_count,
        softclip_fraction=ratio(softclipped, endpoint_count),
        softclipped_count=softclipped,
        secondary_endpoint_count=forward.secondary_endpoint_count + reverse.secondary_endpoint_count,
        indel_position=top_indel.indel_position,
        indel_type=top_indel.indel_type,
        indel_length=top_indel.indel_length,
        indel_read_count=indel_count,
        indel_fraction=ratio(indel_count, local_count),
    )


# The comparison-only RGEN v1.0 score. It never affects pairing or filtering.


class RgenCounts(NamedTuple):
    forward_endpoints: int
    reverse_endpoints: int
    depth: int


def rgen_counts(bam: pysam.AlignmentFile, contig: str, position: int, settings: CallerSettings) -> RgenCounts:
    """Endpoint counts and unstranded depth at one base, using the standalone
    RGEN v1.0 read filter, which keeps supplementary alignments."""
    forward = reverse = depth = 0
    for read in bam.fetch(contig, position, position + 1):
        if (
            read.is_unmapped
            or read.is_secondary
            or read.is_qcfail
            or read.is_duplicate
            or read.mapping_quality < settings.min_mapq
        ):
            continue
        if read.reference_end is None or not read.reference_start <= position < read.reference_end:
            continue
        depth += 1
        for read_end in read_ends(read, settings):
            if read_end.position == position:
                forward += int(read_end.strand == "+")
                reverse += int(read_end.strand == "-")
    return RgenCounts(forward, reverse, depth)


def reads_covering_either(bam: pysam.AlignmentFile, contig: str, first: int, second: int, min_mapq: int) -> int:
    """Counted reads covering `first`, `second`, or both, each counted once."""
    count = 0
    for read in bam.fetch(contig, min(first, second), max(first, second) + 1):
        if is_counted(read, min_mapq) and any(
            read.reference_start <= position < read.reference_end for position in (first, second)
        ):
            count += 1
    return count


def float32(value: float) -> float:
    """Round to IEEE-754 single precision, as the RGEN executable does."""
    return struct.unpack("<f", struct.pack("<f", float(value)))[0]


def rgen_contribution(first_count: int, first_depth: int, second_count: int, second_depth: int) -> float:
    """((c1 - 1) / d1) * ((c2 - 1) / d2) * (c1 + c2 - 2), step by step in float32."""
    if first_depth <= 0 or second_depth <= 0:
        return 0.0
    value = float32(float32(first_count) - float32(1.0))
    value = float32(value / float32(first_depth))
    value = float32(value * float32(float32(second_count) - float32(1.0)))
    value = float32(value / float32(second_depth))
    return float32(value * float32(float32(first_count + second_count) - float32(2.0)))


def rgen_score(
    counts_at: Callable[[int], RgenCounts | None],
    forward_position: int,
    overhang: int,
) -> float:
    """Sum the five-offset RGEN contributions around both anchors.

    `counts_at` returns None for positions off the contig or blacklisted;
    those contribute nothing.
    """
    forward_anchor = counts_at(forward_position)
    reverse_center = forward_position + overhang - 1
    reverse_anchor = counts_at(reverse_center)
    total = float32(0.0)
    for offset in range(-2, 3):
        from_forward = float32(0.0)
        nearby_reverse = counts_at(reverse_center + offset)
        if forward_anchor is not None and nearby_reverse is not None:
            from_forward = rgen_contribution(
                forward_anchor.forward_endpoints, forward_anchor.depth,
                nearby_reverse.reverse_endpoints, nearby_reverse.depth,
            )
        from_reverse = float32(0.0)
        nearby_forward = counts_at(forward_position + offset)
        if reverse_anchor is not None and nearby_forward is not None:
            from_reverse = rgen_contribution(
                reverse_anchor.reverse_endpoints, reverse_anchor.depth,
                nearby_forward.forward_endpoints, nearby_forward.depth,
            )
        total = float32(float32(from_reverse + from_forward) + total)
    return total
