"""nDigenome calls: isolated strand endpoints (single-strand breaks or nicks).

A focal endpoint is called when at least `ndigenome_min_count` reads end on
it and they make up at least `ndigenome_min_fraction` of same-strand depth
(both `>=`). The strongest opposite-strand endpoint within the window sets
the signal class: POSSIBLE_DSB when it passes the same count and fraction,
AMBIGUOUS when it passes the weaker count or fraction, SSB otherwise.

Each chunk builds complete rows for the endpoints it owns; q-values and
filters are applied later, sample-wide, by the finalizer.
"""

from __future__ import annotations

import pysam

from .bam import SiteMetrics, count_endpoints, endpoint_position, is_counted, measure_site, strand_of
from .output import artifact_columns, one_based
from .regions import Blacklist, OwnedInterval, callable_segments
from .settings import CallerSettings
from .stats import UNCONTROLLED, compare_to_control
from .vcf import known_indels


def call_chunk(
    bam: pysam.AlignmentFile,
    control_bam: pysam.AlignmentFile | None,
    vcf: pysam.VariantFile | None,
    settings: CallerSettings,
    intervals: list[OwnedInterval],
    blacklist: Blacklist | None,
    sample: str,
    control_sample: str,
) -> tuple[list[dict], int]:
    """Return the rows for every owned endpoint, and the number of endpoints
    that reached the count threshold."""
    rows = []
    candidate_count = 0
    for interval in intervals:
        for start, end in callable_segments(interval.contig, interval.start, interval.end, blacklist):
            endpoints = count_endpoints(
                bam, interval.contig, start, end,
                settings.ndigenome_min_mapq, settings.ndigenome_min_count,
            )
            candidate_count += len(endpoints)
            for position, strand in sorted(endpoints):
                row = call_endpoint(
                    bam, control_bam, vcf, settings, blacklist,
                    interval.contig, position, strand, sample, control_sample,
                )
                if row is not None:
                    rows.append(row)
    return rows, candidate_count


def call_endpoint(
    bam: pysam.AlignmentFile,
    control_bam: pysam.AlignmentFile | None,
    vcf: pysam.VariantFile | None,
    settings: CallerSettings,
    blacklist: Blacklist | None,
    contig: str,
    position: int,
    strand: str,
    sample: str,
    control_sample: str,
) -> dict | None:
    window = settings.cleavage_artifact_window
    metrics = measure_site(bam, contig, position, strand, window, settings.ndigenome_min_mapq)
    if metrics.endpoint_fraction < settings.ndigenome_min_fraction:
        return None
    opposite_position, opposite = strongest_opposite(bam, contig, position, strand, settings, blacklist)
    nearby_positions = [position] if opposite_position is None else [position, opposite_position]
    return {
        "sample": sample,
        "analysis": "ndigenome",
        "contig": contig,
        "position_0based": position,
        "position_1based": position + 1,
        "strand": strand,
        "signal_classification": signal_class(opposite, settings),
        "endpoint_count": metrics.endpoint_count,
        "strand_depth": metrics.strand_depth,
        "endpoint_fraction": metrics.endpoint_fraction,
        "opposite_position_0based": opposite_position,
        "opposite_position_1based": one_based(opposite_position),
        "opposite_count": opposite.endpoint_count,
        "opposite_depth": opposite.strand_depth,
        "opposite_fraction": opposite.endpoint_fraction,
        **artifact_columns(metrics, known_indels(vcf, contig, nearby_positions, window)),
        **control_columns(control_bam, control_sample, contig, position, strand, metrics, settings),
    }


def passes_primary(metrics: SiteMetrics, settings: CallerSettings) -> bool:
    return (
        metrics.endpoint_count >= settings.ndigenome_min_count
        and metrics.endpoint_fraction >= settings.ndigenome_min_fraction
    )


def passes_ambiguous(metrics: SiteMetrics, settings: CallerSettings) -> bool:
    return (
        metrics.endpoint_count >= settings.ndigenome_ambiguous_min_count
        or metrics.endpoint_fraction >= settings.ndigenome_ambiguous_min_fraction
    )


def signal_class(opposite: SiteMetrics, settings: CallerSettings) -> str:
    if passes_primary(opposite, settings):
        return "POSSIBLE_DSB"
    if passes_ambiguous(opposite, settings):
        return "AMBIGUOUS"
    return "SSB"


def strongest_opposite(
    bam: pysam.AlignmentFile,
    contig: str,
    position: int,
    strand: str,
    settings: CallerSettings,
    blacklist: Blacklist | None,
) -> tuple[int | None, SiteMetrics]:
    """The opposite-strand endpoint within the window that ranks highest by:
    passes the primary thresholds, passes the ambiguous thresholds, count,
    fraction, closeness, then lower coordinate. Blacklisted endpoints are
    ignored. Returns (None, empty metrics) when there is none."""
    opposite_strand = "-" if strand == "+" else "+"
    min_mapq = settings.ndigenome_min_mapq
    start = max(0, position - settings.ndigenome_opposite_window)
    end = position + settings.ndigenome_opposite_window + 1

    candidates = set()
    for read in bam.fetch(contig, start, end):
        if not is_counted(read, min_mapq) or strand_of(read) != opposite_strand:
            continue
        endpoint = endpoint_position(read)
        if blacklist is not None and blacklist.contains(contig, endpoint):
            continue
        if start <= endpoint < end:
            candidates.add(endpoint)

    best = None
    for candidate in sorted(candidates):
        metrics = measure_site(
            bam, contig, candidate, opposite_strand,
            settings.cleavage_artifact_window, min_mapq,
        )
        rank = (
            passes_primary(metrics, settings),
            passes_ambiguous(metrics, settings),
            metrics.endpoint_count,
            metrics.endpoint_fraction,
            -abs(candidate - position),
            -candidate,
        )
        if best is None or rank > best[0]:
            best = (rank, candidate, metrics)
    if best is None:
        return None, SiteMetrics()
    return best[1], best[2]


def control_columns(
    control_bam: pysam.AlignmentFile | None,
    control_sample: str,
    contig: str,
    position: int,
    strand: str,
    treated: SiteMetrics,
    settings: CallerSettings,
) -> dict:
    """The same endpoint measured in the matched control, if there is one."""
    if control_bam is None:
        return {
            "control_sample": "",
            "control_status": UNCONTROLLED,
            "control_endpoint_count": 0,
            "control_depth": 0,
            "control_fraction": 0.0,
            "control_fold_enrichment": None,
            "control_fisher_p": None,
            "control_fisher_q": None,
        }
    control = measure_site(
        control_bam, contig, position, strand,
        settings.cleavage_artifact_window, settings.ndigenome_min_mapq,
    )
    status, fold, p_value = compare_to_control(
        treated.endpoint_count, treated.strand_depth,
        control.endpoint_count, control.strand_depth,
        settings.cleavage_control_min_depth,
    )
    return {
        "control_sample": control_sample,
        "control_status": status,
        "control_endpoint_count": control.endpoint_count,
        "control_depth": control.strand_depth,
        "control_fraction": control.endpoint_fraction,
        "control_fold_enrichment": fold,
        "control_fisher_p": p_value,
        "control_fisher_q": None,
    }
