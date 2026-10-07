"""Genomic regions: the optional BED blacklist and the chunk plan.

Coordinates are 0-based half-open throughout. Each chunk owns a set of
non-overlapping intervals; together the chunks cover every base of every
mapped contig exactly once, blacklisted bases included.
"""

from __future__ import annotations

import gzip
import hashlib
from bisect import bisect_right
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

Segment = tuple[int, int]


@dataclass(frozen=True)
class Blacklist:
    """Merged BED intervals per contig, plus provenance for the QC report."""

    intervals: dict[str, list[Segment]]
    file_name: str
    sha256: str

    @property
    def interval_count(self) -> int:
        return sum(len(segments) for segments in self.intervals.values())

    @property
    def excluded_bases(self) -> int:
        return sum(
            end - start
            for segments in self.intervals.values()
            for start, end in segments
        )

    def contains(self, contig: str, position: int) -> bool:
        segments = self.intervals.get(contig, [])
        # The last segment that starts at or before `position`.
        index = bisect_right(segments, (position, float("inf"))) - 1
        return index >= 0 and position < segments[index][1]

    def subtract(self, contig: str, start: int, end: int) -> list[Segment]:
        """Return the parts of [start, end) outside the blacklist."""
        remaining: list[Segment] = []
        cursor = start
        for blocked_start, blocked_end in self.intervals.get(contig, []):
            if blocked_end <= cursor:
                continue
            if blocked_start >= end:
                break
            if blocked_start > cursor:
                remaining.append((cursor, blocked_start))
            cursor = blocked_end
        if cursor < end:
            remaining.append((cursor, end))
        return remaining

    def provenance(self) -> dict[str, object]:
        return {
            "file": self.file_name,
            "sha256": self.sha256,
            "intervals": self.interval_count,
            "excluded_bases": self.excluded_bases,
        }


def merge_segments(segments: list[Segment]) -> list[Segment]:
    """Merge overlapping or adjacent segments."""
    merged: list[Segment] = []
    for start, end in sorted(segments):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def load_blacklist(path: Path, contig_lengths: dict[str, int]) -> Blacklist:
    """Read a BED or BED.gz file; bad rows, unknown contigs, and
    out-of-range coordinates stop the run."""
    path = Path(path)
    if not path.is_file():
        raise ValueError(f"Genome blacklist does not exist: {path}")
    opener = gzip.open if path.suffix == ".gz" else open
    raw: dict[str, list[Segment]] = {}
    with opener(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            text = line.strip()
            if not text or text.startswith(("#", "track ", "browser ")):
                continue
            contig, start, end = parse_bed_line(text, line_number)
            if contig not in contig_lengths:
                raise ValueError(
                    f"Genome blacklist contig is absent from the BAM: {contig}"
                )
            if not 0 <= start < end <= contig_lengths[contig]:
                raise ValueError(
                    f"Genome blacklist interval is outside {contig}: "
                    f"{start}-{end} of {contig_lengths[contig]}"
                )
            raw.setdefault(contig, []).append((start, end))
    if not raw:
        raise ValueError(f"Genome blacklist contains no BED intervals: {path}")
    return Blacklist(
        intervals={contig: merge_segments(segments) for contig, segments in raw.items()},
        file_name=path.name,
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
    )


def parse_bed_line(text: str, line_number: int) -> tuple[str, int, int]:
    fields = text.split()
    if len(fields) < 3:
        raise ValueError(
            f"Genome blacklist line {line_number} has fewer than three BED columns"
        )
    try:
        return fields[0], int(fields[1]), int(fields[2])
    except ValueError as error:
        raise ValueError(
            f"Genome blacklist line {line_number} has invalid coordinates"
        ) from error


def callable_segments(
    contig: str, start: int, end: int, blacklist: Blacklist | None
) -> list[Segment]:
    if blacklist is None:
        return [(start, end)] if start < end else []
    return blacklist.subtract(contig, start, end)


# Chunk planning.


@dataclass(frozen=True)
class MappedContig:
    name: str
    length: int
    mapped_records: int


@dataclass(frozen=True)
class OwnedInterval:
    contig: str
    start: int
    end: int
    callable_bases: int

    def owns(self, contig: str, position: int) -> bool:
        return contig == self.contig and self.start <= position < self.end


def plan_chunks(
    contigs: list[MappedContig],
    chunk_count: int,
    blacklist: Blacklist | None,
) -> list[list[OwnedInterval]]:
    """Split mapped contigs into `chunk_count` slots of similar callable work.

    Work is a contig's mapped records spread evenly along it, counting only
    bases outside the blacklist. Large contigs are cut at work boundaries and
    small ones share a slot; a slot can be empty when there is too little to
    split. The arithmetic is exact because every chunk task recomputes this
    plan independently and all of them must agree.
    """
    slots: list[list[OwnedInterval]] = [[] for _ in range(chunk_count)]
    segments = {
        contig.name: callable_segments(contig.name, 0, contig.length, blacklist)
        for contig in contigs
    }
    total_work = sum(
        (contig_work(contig, segments[contig.name]) for contig in contigs),
        Fraction(0),
    )
    boundaries = [total_work * index / chunk_count for index in range(1, chunk_count)]

    work_before = Fraction(0)
    for contig in contigs:
        contig_segments = segments[contig.name]
        rate = Fraction(contig.mapped_records, contig.length)
        cuts = cut_points(contig, contig_segments, boundaries, work_before)
        for start, end in zip(cuts, cuts[1:]):
            callable_before = callable_bases_between(contig_segments, 0, start)
            callable_inside = callable_bases_between(contig_segments, start, end)
            midpoint = work_before + rate * (callable_before + Fraction(callable_inside, 2))
            slot = bisect_right(boundaries, midpoint) if total_work else 0
            slots[slot].append(OwnedInterval(contig.name, start, end, callable_inside))
        work_before += contig_work(contig, contig_segments)
    return slots


def contig_work(contig: MappedContig, segments: list[Segment]) -> Fraction:
    bases = sum(end - start for start, end in segments)
    return Fraction(contig.mapped_records * bases, contig.length)


def callable_bases_between(segments: list[Segment], start: int, end: int) -> int:
    return sum(
        max(0, min(end, segment_end) - max(start, segment_start))
        for segment_start, segment_end in segments
    )


def cut_points(
    contig: MappedContig,
    segments: list[Segment],
    boundaries: list[Fraction],
    work_before: Fraction,
) -> list[int]:
    """Coordinates where work boundaries fall inside this contig, plus its ends."""
    rate = Fraction(contig.mapped_records, contig.length)
    work_after = work_before + contig_work(contig, segments)
    cuts = [0]
    for boundary in boundaries:
        if not work_before < boundary < work_after:
            continue
        coordinate = coordinate_at_callable_offset(segments, (boundary - work_before) / rate)
        coordinate = max(cuts[-1] + 1, min(contig.length - 1, coordinate))
        if coordinate < contig.length:
            cuts.append(coordinate)
    cuts.append(contig.length)
    return cuts


def coordinate_at_callable_offset(segments: list[Segment], offset: Fraction) -> int:
    """Walk `offset` callable bases into the segments, rounding half up."""
    for start, end in segments:
        if offset <= end - start:
            return start + (offset.numerator + offset.denominator // 2) // offset.denominator
        offset -= end - start
    return segments[-1][1]


def check_plan(plan: list[list[OwnedInterval]], contig_lengths: dict[str, int]) -> None:
    """Every mapped contig must be owned exactly once, end to end."""
    owned: dict[str, list[tuple[int, int, int]]] = {name: [] for name in contig_lengths}
    for slot, intervals in enumerate(plan):
        for interval in intervals:
            if interval.contig not in owned:
                raise ValueError(f"Chunk {slot} owns an unexpected contig: {interval.contig}")
            owned[interval.contig].append((interval.start, interval.end, slot))
    for contig, length in contig_lengths.items():
        cursor = 0
        for start, end, slot in sorted(owned[contig]):
            if start > cursor:
                raise ValueError(f"Chunk ownership has a gap on {contig}: {cursor}-{start}")
            if start < cursor:
                raise ValueError(
                    f"Chunk ownership overlaps on {contig} at {start}-{min(cursor, end)} "
                    f"(chunk {slot})"
                )
            cursor = end
        if cursor != length:
            raise ValueError(f"Chunk ownership has a gap on {contig}: {cursor}-{length}")


def write_plan(plan: list[list[OwnedInterval]], path: str) -> None:
    """One row per owned interval, for checking how work was split."""
    lines = ["chunk\tcontig\tstart\tend\tcallable_bases\texcluded_bases"]
    for slot, intervals in enumerate(plan):
        for interval in intervals:
            excluded = interval.end - interval.start - interval.callable_bases
            lines.append(
                f"{slot}\t{interval.contig}\t{interval.start}\t{interval.end}\t"
                f"{interval.callable_bases}\t{excluded}"
            )
    Path(path).write_text("\n".join(lines) + "\n")


def plan_to_json(plan: list[list[OwnedInterval]]) -> list[list[list[object]]]:
    return [
        [[interval.contig, interval.start, interval.end, interval.callable_bases] for interval in slot]
        for slot in plan
    ]


def plan_from_json(document: list[list[list[object]]]) -> list[list[OwnedInterval]]:
    return [
        [OwnedInterval(str(contig), int(start), int(end), int(bases)) for contig, start, end, bases in slot]
        for slot in document
    ]
