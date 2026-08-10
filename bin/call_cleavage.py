#!/usr/bin/env python3
"""Call strand-aware Digenome-seq DSB and nDigenome-seq SSB endpoints."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import heapq
import json
import math
import statistics
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from genome_blacklist import GenomeBlacklist, load_genome_blacklist

try:
    import pysam
except ImportError:  # pragma: no cover - exercised by the command-line guard
    pysam = None


ANALYSES = {"digenome", "ndigenome"}

TSV_FIELDS = [
    "sample",
    "analysis",
    "control_sample",
    "control_status",
    "contig",
    "position_0based",
    "position_1based",
    "strand",
    "forward_position_0based",
    "forward_position_1based",
    "reverse_position_0based",
    "reverse_position_1based",
    "endpoint_count",
    "strand_depth",
    "endpoint_fraction",
    "opposite_position_0based",
    "opposite_position_1based",
    "opposite_count",
    "opposite_depth",
    "opposite_fraction",
    "forward_endpoint_count",
    "forward_depth",
    "forward_fraction",
    "reverse_endpoint_count",
    "reverse_depth",
    "reverse_fraction",
    "digenome_score",
    "signal_classification",
    "classification",
    "filter_status",
    "filter_reasons",
    "support_mean_mapq",
    "local_mean_mapq",
    "support_mean_nm",
    "local_mean_nm",
    "softclip_fraction",
    "secondary_endpoint_count",
    "indel_position_0based",
    "indel_position_1based",
    "indel_type",
    "indel_length",
    "indel_read_count",
    "indel_fraction",
    "known_indel_overlap",
    "control_endpoint_count",
    "control_depth",
    "control_fraction",
    "control_forward_endpoint_count",
    "control_forward_depth",
    "control_forward_fraction",
    "control_reverse_endpoint_count",
    "control_reverse_depth",
    "control_reverse_fraction",
    "control_digenome_score",
    "control_fold_enrichment",
    "control_fisher_p",
    "control_fisher_q",
]

FLOAT_FIELDS = {
    "endpoint_fraction",
    "opposite_fraction",
    "forward_fraction",
    "reverse_fraction",
    "digenome_score",
    "support_mean_mapq",
    "local_mean_mapq",
    "support_mean_nm",
    "local_mean_nm",
    "softclip_fraction",
    "indel_fraction",
    "control_fraction",
    "control_forward_fraction",
    "control_reverse_fraction",
    "control_digenome_score",
    "control_fold_enrichment",
    "control_fisher_p",
    "control_fisher_q",
}


@dataclass(frozen=True)
class GenomicInterval:
    contig: str
    start: int
    end: int
    scan_start: int
    scan_end: int

    def owns(self, contig: str, position: int) -> bool:
        return (
            self.contig == contig
            and self.start <= position < self.end
        )


# Basic alignment and endpoint measurements.

def mean_or_zero(values: Iterable[float]) -> float:
    values = list(values)
    return statistics.fmean(values) if values else 0.0


def weighted_mean(
    first_value: float,
    first_count: int,
    second_value: float,
    second_count: int,
) -> float:
    total = first_count + second_count
    if total == 0:
        return 0.0
    return (
        first_value * first_count + second_value * second_count
    ) / total


def eligible_primary(read: Any, min_mapq: int = 1) -> bool:
    return not (
        read.is_unmapped
        or read.is_secondary
        or read.is_supplementary
        or read.is_qcfail
        or read.is_duplicate
        or read.mapping_quality < min_mapq
    )


def strand_symbol(read: Any) -> str:
    return "-" if read.is_reverse else "+"


def endpoint_position(read: Any) -> int:
    """Return the aligned 5-prime endpoint as a 0-based reference coordinate."""
    if read.reference_start is None or read.reference_end is None:
        raise ValueError("Cannot calculate an endpoint for an unmapped alignment")
    return read.reference_end - 1 if read.is_reverse else read.reference_start


def five_prime_softclip_length(read: Any) -> int:
    cigars = read.cigartuples or []
    if not cigars:
        return 0
    op, length = cigars[-1] if read.is_reverse else cigars[0]
    return length if op in (4, 5) else 0


def cigar_indel_events(read: Any) -> list[tuple[int, str, int]]:
    events: list[tuple[int, str, int]] = []
    reference_pos = read.reference_start
    if reference_pos is None:
        return events
    for operation, length in read.cigartuples or []:
        if operation == 1:
            events.append((reference_pos, "INS", length))
        elif operation == 2:
            events.append((reference_pos, "DEL", length))
            reference_pos += length
        elif operation in (0, 3, 7, 8):
            reference_pos += length
    return events


def scan_candidate_endpoints(
    bam: Any,
    min_count: int,
    min_mapq: int,
    contigs: Iterable[str] | None = None,
    intervals: Iterable[GenomicInterval] | None = None,
    genome_blacklist: GenomeBlacklist | None = None,
) -> list[tuple[str, int, str, int]]:
    """Stream a coordinate-sorted BAM and retain endpoint runs above cutoff."""
    if contigs is not None and intervals is not None:
        raise ValueError("Specify contigs or intervals, not both")

    if intervals is None:
        selected_contigs = (
            list(contigs) if contigs is not None else list(bam.references)
        )
        lengths = dict(zip(bam.references, bam.lengths))
        selected_intervals = [
            GenomicInterval(
                contig=contig,
                start=0,
                end=int(lengths[contig]),
                scan_start=0,
                scan_end=int(lengths[contig]),
            )
            for contig in selected_contigs
            if contig in lengths
        ]
    else:
        selected_intervals = list(intervals)
        selected_contigs = [
            interval.contig for interval in selected_intervals
        ]

    unknown = sorted(set(selected_contigs) - set(bam.references))
    if unknown:
        raise ValueError(
            "Requested contig(s) are absent from the BAM: "
            + ", ".join(unknown)
        )
    candidates: dict[tuple[str, int, str], int] = {}
    for interval in selected_intervals:
        contig = interval.contig
        scan_regions = (
            genome_blacklist.subtract(
                contig,
                interval.scan_start,
                interval.scan_end,
            )
            if genome_blacklist is not None
            else [(interval.scan_start, interval.scan_end)]
        )
        for scan_start, scan_end in scan_regions:
            forward_position: int | None = None
            forward_count = 0
            reverse_counts: dict[int, int] = {}
            reverse_heap: list[int] = []
            previous_start = -1

            def save_forward() -> None:
                if (
                    forward_position is not None
                    and forward_count >= min_count
                ):
                    key = (contig, forward_position, "+")
                    candidates[key] = max(
                        candidates.get(key, 0),
                        forward_count,
                    )

            def save_reverse(position: int) -> None:
                count = reverse_counts.pop(position, 0)
                if count >= min_count:
                    key = (contig, position, "-")
                    candidates[key] = max(candidates.get(key, 0), count)

            for read in bam.fetch(contig, scan_start, scan_end):
                if not eligible_primary(read, min_mapq):
                    continue
                if read.reference_start < previous_start:
                    raise ValueError(
                        f"BAM is not coordinate sorted on {contig}: "
                        f"{read.reference_start} follows {previous_start}"
                    )
                previous_start = read.reference_start

                while (
                    reverse_heap
                    and reverse_heap[0] < read.reference_start
                ):
                    save_reverse(heapq.heappop(reverse_heap))

                endpoint = endpoint_position(read)
                if not scan_start <= endpoint < scan_end:
                    continue
                if read.is_reverse:
                    if endpoint not in reverse_counts:
                        reverse_counts[endpoint] = 0
                        heapq.heappush(reverse_heap, endpoint)
                    reverse_counts[endpoint] += 1
                elif forward_position == endpoint:
                    forward_count += 1
                else:
                    save_forward()
                    forward_position = endpoint
                    forward_count = 1

            save_forward()
            while reverse_heap:
                save_reverse(heapq.heappop(reverse_heap))

    return [
        (*key, candidates[key])
        for key in sorted(candidates)
    ]


def measure_site_metrics(
    bam: Any,
    contig: str,
    position: int,
    strand: str,
    artifact_window: int,
    min_mapq: int,
) -> dict[str, Any]:
    start = max(0, position - artifact_window)
    end = position + artifact_window + 1
    endpoint_reads: list[Any] = []
    local_primary: list[Any] = []
    local_mapqs: list[int] = []
    local_nms: list[float] = []
    strand_depth = 0
    secondary_endpoint_count = 0
    indel_alignments: set[tuple[str, int, int]] = set()
    indel_events: Counter[tuple[int, str, int]] = Counter()

    for read in bam.fetch(contig, start, end):
        if read.is_unmapped or read.is_qcfail:
            continue
        if read.is_secondary or read.is_supplementary:
            if (
                read.mapping_quality >= min_mapq
                and strand_symbol(read) == strand
                and endpoint_position(read) == position
            ):
                secondary_endpoint_count += 1
            continue
        if not eligible_primary(read, min_mapq):
            continue

        local_primary.append(read)
        local_mapqs.append(read.mapping_quality)
        if read.has_tag("NM"):
            local_nms.append(float(read.get_tag("NM")))
        if (
            strand_symbol(read) == strand
            and read.reference_start <= position < read.reference_end
        ):
            strand_depth += 1
        if strand_symbol(read) == strand and endpoint_position(read) == position:
            endpoint_reads.append(read)

        read_has_indel = False
        for event in cigar_indel_events(read):
            event_position, kind, length = event
            event_end = event_position + (length if kind == "DEL" else 0)
            overlaps_window = (
                event_position < end and event_end > start
                if kind == "DEL"
                else start <= event_position < end
            )
            if overlaps_window:
                indel_events[event] += 1
                read_has_indel = True
        if read_has_indel:
            indel_alignments.add(
                (read.query_name or "", read.flag, read.reference_start)
            )

    support_mapqs = [read.mapping_quality for read in endpoint_reads]
    support_nms = [
        float(read.get_tag("NM")) for read in endpoint_reads if read.has_tag("NM")
    ]
    softclipped = sum(
        1 for read in endpoint_reads if five_prime_softclip_length(read) > 0
    )
    top_indel = indel_events.most_common(1)
    indel_position: int | None = None
    indel_type = ""
    indel_length = 0
    if top_indel:
        (indel_position, indel_type, indel_length), _count = top_indel[0]

    return {
        "endpoint_count": len(endpoint_reads),
        "strand_depth": strand_depth,
        "endpoint_fraction": (
            len(endpoint_reads) / strand_depth if strand_depth else 0.0
        ),
        "support_mean_mapq": mean_or_zero(support_mapqs),
        "support_read_count": len(support_mapqs),
        "local_mean_mapq": mean_or_zero(local_mapqs),
        "local_read_count": len(local_mapqs),
        "support_mean_nm": mean_or_zero(support_nms),
        "support_nm_count": len(support_nms),
        "local_mean_nm": mean_or_zero(local_nms),
        "local_nm_count": len(local_nms),
        "softclip_fraction": (
            softclipped / len(endpoint_reads) if endpoint_reads else 0.0
        ),
        "softclipped_count": softclipped,
        "secondary_endpoint_count": secondary_endpoint_count,
        "indel_position": indel_position,
        "indel_type": indel_type,
        "indel_length": indel_length,
        "indel_read_count": len(indel_alignments),
        "indel_fraction": (
            len(indel_alignments) / len(local_primary) if local_primary else 0.0
        ),
    }


def combine_site_metrics(
    forward: dict[str, Any], reverse: dict[str, Any]
) -> dict[str, Any]:
    endpoint_count = forward["endpoint_count"] + reverse["endpoint_count"]
    strand_depth = forward["strand_depth"] + reverse["strand_depth"]
    local_count = forward["local_read_count"] + reverse["local_read_count"]
    indel_count = forward["indel_read_count"] + reverse["indel_read_count"]
    top_indel = max(
        (forward, reverse),
        key=lambda metrics: (
            metrics["indel_fraction"],
            metrics["indel_read_count"],
        ),
    )
    return {
        "endpoint_count": endpoint_count,
        "strand_depth": strand_depth,
        "endpoint_fraction": (
            endpoint_count / strand_depth if strand_depth else 0.0
        ),
        "support_mean_mapq": weighted_mean(
            forward["support_mean_mapq"],
            forward["support_read_count"],
            reverse["support_mean_mapq"],
            reverse["support_read_count"],
        ),
        "support_read_count": (
            forward["support_read_count"] + reverse["support_read_count"]
        ),
        "local_mean_mapq": weighted_mean(
            forward["local_mean_mapq"],
            forward["local_read_count"],
            reverse["local_mean_mapq"],
            reverse["local_read_count"],
        ),
        "local_read_count": local_count,
        "support_mean_nm": weighted_mean(
            forward["support_mean_nm"],
            forward["support_nm_count"],
            reverse["support_mean_nm"],
            reverse["support_nm_count"],
        ),
        "support_nm_count": (
            forward["support_nm_count"] + reverse["support_nm_count"]
        ),
        "local_mean_nm": weighted_mean(
            forward["local_mean_nm"],
            forward["local_nm_count"],
            reverse["local_mean_nm"],
            reverse["local_nm_count"],
        ),
        "local_nm_count": forward["local_nm_count"] + reverse["local_nm_count"],
        "softclip_fraction": (
            (forward["softclipped_count"] + reverse["softclipped_count"])
            / endpoint_count
            if endpoint_count
            else 0.0
        ),
        "softclipped_count": (
            forward["softclipped_count"] + reverse["softclipped_count"]
        ),
        "secondary_endpoint_count": (
            forward["secondary_endpoint_count"]
            + reverse["secondary_endpoint_count"]
        ),
        "indel_position": top_indel["indel_position"],
        "indel_type": top_indel["indel_type"],
        "indel_length": top_indel["indel_length"],
        "indel_read_count": indel_count,
        "indel_fraction": indel_count / local_count if local_count else 0.0,
    }


def find_best_opposite_signal(
    bam: Any,
    contig: str,
    position: int,
    strand: str,
    window: int,
    artifact_window: int,
    min_mapq: int,
    primary_min_count: int,
    primary_min_fraction: float,
    ambiguous_min_count: int,
    ambiguous_min_fraction: float,
    genome_blacklist: GenomeBlacklist | None = None,
) -> tuple[int | None, dict[str, Any]]:
    opposite = "-" if strand == "+" else "+"
    counts: Counter[int] = Counter()
    start = max(0, position - window)
    end = position + window + 1
    for read in bam.fetch(contig, start, end):
        if not eligible_primary(read, min_mapq) or strand_symbol(read) != opposite:
            continue
        endpoint = endpoint_position(read)
        if (
            genome_blacklist is not None
            and genome_blacklist.contains(contig, endpoint)
        ):
            continue
        if start <= endpoint < end:
            counts[endpoint] += 1
    if not counts:
        return None, empty_metrics()

    ranked_signals: list[tuple[tuple[Any, ...], int, dict[str, Any]]] = []
    for candidate_position in sorted(counts):
        candidate_metrics = measure_site_metrics(
            bam,
            contig,
            candidate_position,
            opposite,
            artifact_window,
            min_mapq,
        )
        passes_primary_threshold = (
            candidate_metrics["endpoint_count"] >= primary_min_count
            and candidate_metrics["endpoint_fraction"] >= primary_min_fraction
        )
        passes_ambiguous_threshold = (
            candidate_metrics["endpoint_count"] >= ambiguous_min_count
            or candidate_metrics["endpoint_fraction"]
            >= ambiguous_min_fraction
        )
        rank = (
            passes_primary_threshold,
            passes_ambiguous_threshold,
            candidate_metrics["endpoint_count"],
            candidate_metrics["endpoint_fraction"],
            -abs(candidate_position - position),
            -candidate_position,
        )
        ranked_signals.append((rank, candidate_position, candidate_metrics))

    _rank, opposite_position, opposite_metrics = max(ranked_signals)
    return opposite_position, opposite_metrics


def empty_metrics() -> dict[str, Any]:
    return {
        "endpoint_count": 0,
        "strand_depth": 0,
        "endpoint_fraction": 0.0,
        "support_mean_mapq": 0.0,
        "support_read_count": 0,
        "local_mean_mapq": 0.0,
        "local_read_count": 0,
        "support_mean_nm": 0.0,
        "support_nm_count": 0,
        "local_mean_nm": 0.0,
        "local_nm_count": 0,
        "softclip_fraction": 0.0,
        "softclipped_count": 0,
        "secondary_endpoint_count": 0,
        "indel_position": None,
        "indel_type": "",
        "indel_length": 0,
        "indel_read_count": 0,
        "indel_fraction": 0.0,
    }


def known_indels(
    variant_file: Any | None,
    contig: str,
    positions: Iterable[int],
    window: int,
) -> list[str]:
    if variant_file is None:
        return []
    overlaps: set[str] = set()
    for position in positions:
        start = max(0, position - window)
        end = position + window + 1
        records = variant_file.fetch(contig, start, end)
        for record in records:
            alts = record.alts or ()
            indel_alts = [
                alt for alt in alts if alt and len(alt) != len(record.ref)
            ]
            if indel_alts:
                overlaps.add(
                    f"{contig}:{record.pos}:{record.ref}>{','.join(indel_alts)}"
                )
    return sorted(overlaps)


def digenome_score(
    forward_count: int,
    forward_fraction: float,
    reverse_count: int,
    reverse_fraction: float,
) -> float:
    """Return the published Digenome score using fractional strand ratios."""
    return (
        forward_fraction
        * reverse_fraction
        * (forward_count + reverse_count)
        / 4.0
    )


# One-to-one Digenome endpoint pairing.

@dataclass
class DigenomePairCandidate:
    contig: str
    forward_position: int
    reverse_position: int
    forward_metrics: dict[str, Any]
    reverse_metrics: dict[str, Any]
    score: float
    caller_filter_reasons: list[str]

    @property
    def passes_caller_thresholds(self) -> bool:
        return not self.caller_filter_reasons


@dataclass
class MatchingEdge:
    destination: int
    reverse_edge_index: int
    capacity: int
    cost: tuple[int, int, int, int]


def group_digenome_pair_components(
    candidates: list[DigenomePairCandidate],
) -> list[list[DigenomePairCandidate]]:
    candidates_by_contig: dict[str, list[DigenomePairCandidate]] = {}
    for candidate in candidates:
        candidates_by_contig.setdefault(candidate.contig, []).append(candidate)

    independent_components: list[list[DigenomePairCandidate]] = []
    for contig in sorted(candidates_by_contig):
        contig_candidates = sorted(
            candidates_by_contig[contig],
            key=lambda candidate: (
                candidate.forward_position,
                candidate.reverse_position,
            ),
        )
        candidates_for_forward: dict[int, list[int]] = {}
        candidates_for_reverse: dict[int, list[int]] = {}
        for index, candidate in enumerate(contig_candidates):
            candidates_for_forward.setdefault(
                candidate.forward_position,
                [],
            ).append(index)
            candidates_for_reverse.setdefault(
                candidate.reverse_position,
                [],
            ).append(index)

        unvisited = set(range(len(contig_candidates)))
        while unvisited:
            pending = [min(unvisited)]
            component_indexes: list[int] = []
            while pending:
                index = pending.pop()
                if index not in unvisited:
                    continue
                unvisited.remove(index)
                component_indexes.append(index)
                candidate = contig_candidates[index]
                connected_indexes = (
                    candidates_for_forward[candidate.forward_position]
                    + candidates_for_reverse[candidate.reverse_position]
                )
                pending.extend(
                    connected_index
                    for connected_index in reversed(connected_indexes)
                    if connected_index in unvisited
                )
            independent_components.append(
                [
                    contig_candidates[index]
                    for index in sorted(component_indexes)
                ]
            )
    return independent_components


def select_digenome_pairs(
    candidates: list[DigenomePairCandidate],
) -> list[DigenomePairCandidate]:
    """Choose a deterministic maximum-quality one-to-one endpoint matching."""
    selected: list[DigenomePairCandidate] = []
    independent_components = group_digenome_pair_components(candidates)

    for contig_candidates in independent_components:
        forward_positions = sorted(
            {candidate.forward_position for candidate in contig_candidates}
        )
        reverse_positions = sorted(
            {candidate.reverse_position for candidate in contig_candidates}
        )

        source = 0
        first_forward = 1
        first_reverse = first_forward + len(forward_positions)
        sink = first_reverse + len(reverse_positions)
        graph: list[list[MatchingEdge]] = [
            [] for _node in range(sink + 1)
        ]

        def add_edge(
            start: int,
            destination: int,
            cost: tuple[int, int, int, int],
        ) -> MatchingEdge:
            forward_edge = MatchingEdge(
                destination=destination,
                reverse_edge_index=len(graph[destination]),
                capacity=1,
                cost=cost,
            )
            reverse_edge = MatchingEdge(
                destination=start,
                reverse_edge_index=len(graph[start]),
                capacity=0,
                cost=tuple(-value for value in cost),
            )
            graph[start].append(forward_edge)
            graph[destination].append(reverse_edge)
            return forward_edge

        zero_cost = (0, 0, 0, 0)
        forward_nodes = {
            position: first_forward + index
            for index, position in enumerate(forward_positions)
        }
        reverse_nodes = {
            position: first_reverse + index
            for index, position in enumerate(reverse_positions)
        }
        for position in forward_positions:
            add_edge(source, forward_nodes[position], zero_cost)
        for position in reverse_positions:
            add_edge(reverse_nodes[position], sink, zero_cost)

        candidate_edges: list[
            tuple[MatchingEdge, DigenomePairCandidate]
        ] = []
        score_ratios = [
            float(candidate.score).as_integer_ratio()
            for candidate in contig_candidates
        ]
        # Binary-float denominators are powers of two, so the largest is
        # also their least common multiple.
        common_score_denominator = max(
            denominator for _numerator, denominator in score_ratios
        )
        for coordinate_rank, (candidate, score_ratio) in enumerate(
            zip(contig_candidates, score_ratios)
        ):
            score_numerator, score_denominator = score_ratio
            score_units = score_numerator * (
                common_score_denominator // score_denominator
            )
            reward_cost = (
                -int(candidate.passes_caller_thresholds),
                -1,
                -score_units,
                coordinate_rank,
            )
            edge = add_edge(
                forward_nodes[candidate.forward_position],
                reverse_nodes[candidate.reverse_position],
                reward_cost,
            )
            candidate_edges.append((edge, candidate))

        while True:
            distances: list[tuple[int, int, int, int] | None] = [
                None
            ] * len(graph)
            previous: list[tuple[int, int] | None] = [None] * len(graph)
            distances[source] = zero_cost

            for _iteration in range(len(graph) - 1):
                changed = False
                for node, edges in enumerate(graph):
                    if distances[node] is None:
                        continue
                    for edge_index, edge in enumerate(edges):
                        if edge.capacity == 0:
                            continue
                        new_distance = tuple(
                            left + right
                            for left, right in zip(
                                distances[node],
                                edge.cost,
                            )
                        )
                        if (
                            distances[edge.destination] is None
                            or new_distance < distances[edge.destination]
                        ):
                            distances[edge.destination] = new_distance
                            previous[edge.destination] = (
                                node,
                                edge_index,
                            )
                            changed = True
                if not changed:
                    break

            if distances[sink] is None or distances[sink] >= zero_cost:
                break

            node = sink
            path_nodes: set[int] = set()
            while node != source:
                if node in path_nodes:
                    raise RuntimeError(
                        "Internal Digenome matching path contains a cycle"
                    )
                path_nodes.add(node)
                previous_step = previous[node]
                if previous_step is None:
                    raise RuntimeError(
                        "Internal Digenome matching path is incomplete"
                    )
                previous_node, edge_index = previous_step
                edge = graph[previous_node][edge_index]
                reverse_edge = graph[edge.destination][
                    edge.reverse_edge_index
                ]
                edge.capacity = 0
                reverse_edge.capacity = 1
                node = previous_node

        selected.extend(
            candidate
            for edge, candidate in candidate_edges
            if edge.capacity == 0
        )

    return sorted(
        selected,
        key=lambda candidate: (
            candidate.contig,
            candidate.forward_position,
            candidate.reverse_position,
        ),
    )


def create_base_call_row(
    args: argparse.Namespace,
    contig: str,
    position: int,
    strand: str,
    metrics: dict[str, Any],
) -> dict[str, Any]:
    return {
        "sample": args.sample,
        "analysis": args.analysis,
        "control_sample": "",
        "control_status": "UNCONTROLLED",
        "contig": contig,
        "position_0based": position,
        "position_1based": position + 1,
        "strand": strand,
        "endpoint_count": metrics["endpoint_count"],
        "strand_depth": metrics["strand_depth"],
        "endpoint_fraction": metrics["endpoint_fraction"],
        "support_mean_mapq": metrics["support_mean_mapq"],
        "local_mean_mapq": metrics["local_mean_mapq"],
        "support_mean_nm": metrics["support_mean_nm"],
        "local_mean_nm": metrics["local_mean_nm"],
        "softclip_fraction": metrics["softclip_fraction"],
        "secondary_endpoint_count": metrics["secondary_endpoint_count"],
        "indel_position_0based": metrics["indel_position"],
        "indel_position_1based": (
            metrics["indel_position"] + 1
            if metrics["indel_position"] is not None
            else None
        ),
        "indel_type": metrics["indel_type"],
        "indel_length": metrics["indel_length"],
        "indel_read_count": metrics["indel_read_count"],
        "indel_fraction": metrics["indel_fraction"],
        "known_indel_overlap": "",
        "control_endpoint_count": 0,
        "control_depth": 0,
        "control_fraction": 0.0,
        "control_forward_endpoint_count": 0,
        "control_forward_depth": 0,
        "control_forward_fraction": 0.0,
        "control_reverse_endpoint_count": 0,
        "control_reverse_depth": 0,
        "control_reverse_fraction": 0.0,
        "control_digenome_score": None,
        "control_fold_enrichment": None,
        "control_fisher_p": None,
        "control_fisher_q": None,
        "_caller_filter_reasons": [],
    }


def position_is_owned(
    intervals: Iterable[GenomicInterval] | None,
    contig: str,
    position: int,
    genome_blacklist: GenomeBlacklist | None = None,
) -> bool:
    return (
        (
            intervals is None
            or any(
                interval.owns(contig, position)
                for interval in intervals
            )
        )
        and (
            genome_blacklist is None
            or not genome_blacklist.contains(contig, position)
        )
    )


def call_ndigenome(
    bam: Any,
    variant_file: Any | None,
    args: argparse.Namespace,
    contigs: Iterable[str] | None = None,
    intervals: Iterable[GenomicInterval] | None = None,
    genome_blacklist: GenomeBlacklist | None = None,
) -> tuple[list[dict[str, Any]], int]:
    selected_intervals = (
        list(intervals) if intervals is not None else None
    )
    scanned = scan_candidate_endpoints(
        bam,
        args.ndigenome_min_count,
        args.ndigenome_min_mapq,
        contigs,
        selected_intervals,
        genome_blacklist,
    )
    owned_scanned = [
        candidate
        for candidate in scanned
        if position_is_owned(
            selected_intervals,
            candidate[0],
            candidate[1],
            genome_blacklist,
        )
    ]
    rows: list[dict[str, Any]] = []
    for contig, position, strand, scanned_count in owned_scanned:
        metrics = measure_site_metrics(
            bam,
            contig,
            position,
            strand,
            args.artifact_window,
            args.ndigenome_min_mapq,
        )
        if (
            scanned_count < args.ndigenome_min_count
            or metrics["endpoint_fraction"] < args.ndigenome_min_fraction
        ):
            continue

        opposite_position, opposite = find_best_opposite_signal(
            bam,
            contig,
            position,
            strand,
            args.ndigenome_opposite_window,
            args.artifact_window,
            args.ndigenome_min_mapq,
            args.ndigenome_min_count,
            args.ndigenome_min_fraction,
            args.ndigenome_ambiguous_min_count,
            args.ndigenome_ambiguous_min_fraction,
            genome_blacklist,
        )
        if (
            opposite["endpoint_count"] >= args.ndigenome_min_count
            and opposite["endpoint_fraction"] >= args.ndigenome_min_fraction
        ):
            signal_classification = "POSSIBLE_DSB"
        elif (
            opposite["endpoint_count"] >= args.ndigenome_ambiguous_min_count
            or opposite["endpoint_fraction"]
            >= args.ndigenome_ambiguous_min_fraction
        ):
            signal_classification = "AMBIGUOUS"
        else:
            signal_classification = "SSB"

        row = create_base_call_row(
            args,
            contig,
            position,
            strand,
            metrics,
        )
        row.update(
            {
                "opposite_position_0based": opposite_position,
                "opposite_position_1based": (
                    opposite_position + 1
                    if opposite_position is not None
                    else None
                ),
                "opposite_count": opposite["endpoint_count"],
                "opposite_depth": opposite["strand_depth"],
                "opposite_fraction": opposite["endpoint_fraction"],
                "signal_classification": signal_classification,
            }
        )
        if strand == "+":
            row.update(
                {
                    "forward_position_0based": position,
                    "forward_position_1based": position + 1,
                    "reverse_position_0based": opposite_position,
                    "reverse_position_1based": (
                        opposite_position + 1
                        if opposite_position is not None
                        else None
                    ),
                    "forward_endpoint_count": metrics["endpoint_count"],
                    "forward_depth": metrics["strand_depth"],
                    "forward_fraction": metrics["endpoint_fraction"],
                    "reverse_endpoint_count": opposite["endpoint_count"],
                    "reverse_depth": opposite["strand_depth"],
                    "reverse_fraction": opposite["endpoint_fraction"],
                }
            )
        else:
            row.update(
                {
                    "forward_position_0based": opposite_position,
                    "forward_position_1based": (
                        opposite_position + 1
                        if opposite_position is not None
                        else None
                    ),
                    "reverse_position_0based": position,
                    "reverse_position_1based": position + 1,
                    "forward_endpoint_count": opposite["endpoint_count"],
                    "forward_depth": opposite["strand_depth"],
                    "forward_fraction": opposite["endpoint_fraction"],
                    "reverse_endpoint_count": metrics["endpoint_count"],
                    "reverse_depth": metrics["strand_depth"],
                    "reverse_fraction": metrics["endpoint_fraction"],
                }
            )
        row["digenome_score"] = (
            digenome_score(
                row["forward_endpoint_count"],
                row["forward_fraction"],
                row["reverse_endpoint_count"],
                row["reverse_fraction"],
            )
            if opposite_position is not None
            else 0.0
        )
        row["known_indel_overlap"] = ";".join(
            known_indels(
                variant_file,
                contig,
                [
                    candidate
                    for candidate in (position, opposite_position)
                    if candidate is not None
                ],
                args.artifact_window,
            )
        )
        rows.append(row)
    return rows, len(owned_scanned)


def call_digenome(
    bam: Any,
    variant_file: Any | None,
    args: argparse.Namespace,
    contigs: Iterable[str] | None = None,
    intervals: Iterable[GenomicInterval] | None = None,
    genome_blacklist: GenomeBlacklist | None = None,
) -> tuple[list[dict[str, Any]], int]:
    selected_intervals = (
        list(intervals) if intervals is not None else None
    )
    scan_minimum = min(
        args.digenome_forward_cutoff + 1,
        args.digenome_reverse_cutoff + 1,
    )
    site_metric_cache: dict[tuple[str, int, str], dict[str, Any]] = {}

    def get_site_metrics(
        contig: str,
        position: int,
        strand: str,
    ) -> dict[str, Any]:
        key = (contig, position, strand)
        if key not in site_metric_cache:
            site_metric_cache[key] = measure_site_metrics(
                bam,
                contig,
                position,
                strand,
                args.artifact_window,
                args.digenome_min_mapq,
            )
        return site_metric_cache[key]

    def build_pair_candidates(
        selected_contigs: Iterable[str] | None = None,
        scan_intervals: Iterable[GenomicInterval] | None = None,
    ) -> list[DigenomePairCandidate]:
        scanned = scan_candidate_endpoints(
            bam,
            scan_minimum,
            args.digenome_min_mapq,
            selected_contigs,
            scan_intervals,
            genome_blacklist,
        )
        by_contig: dict[str, dict[str, list[tuple[int, int]]]] = {}
        for contig, position, strand, count in scanned:
            by_contig.setdefault(
                contig,
                {"+": [], "-": []},
            )[strand].append((position, count))

        pair_candidates: list[DigenomePairCandidate] = []
        for contig, strands in by_contig.items():
            reverse = strands["-"]
            reverse_positions = [position for position, _count in reverse]
            for forward_position, _forward_count in strands["+"]:
                expected_reverse = (
                    forward_position - args.digenome_overhang
                )
                lower = expected_reverse - args.digenome_pair_window
                upper = expected_reverse + args.digenome_pair_window
                start = bisect.bisect_left(reverse_positions, lower)
                end = bisect.bisect_right(reverse_positions, upper)
                for reverse_position, _reverse_count in reverse[start:end]:
                    forward_metrics = get_site_metrics(
                        contig,
                        forward_position,
                        "+",
                    )
                    reverse_metrics = get_site_metrics(
                        contig,
                        reverse_position,
                        "-",
                    )
                    score = digenome_score(
                        forward_metrics["endpoint_count"],
                        forward_metrics["endpoint_fraction"],
                        reverse_metrics["endpoint_count"],
                        reverse_metrics["endpoint_fraction"],
                    )
                    caller_reasons: list[str] = []
                    if (
                        forward_metrics["endpoint_count"]
                        <= args.digenome_forward_cutoff
                    ):
                        caller_reasons.append("LOW_FORWARD_COUNT")
                    if (
                        reverse_metrics["endpoint_count"]
                        <= args.digenome_reverse_cutoff
                    ):
                        caller_reasons.append("LOW_REVERSE_COUNT")
                    if (
                        forward_metrics["strand_depth"]
                        <= args.digenome_depth_cutoff
                    ):
                        caller_reasons.append("LOW_FORWARD_DEPTH")
                    if (
                        reverse_metrics["strand_depth"]
                        <= args.digenome_depth_cutoff
                    ):
                        caller_reasons.append("LOW_REVERSE_DEPTH")
                    if (
                        forward_metrics["endpoint_fraction"]
                        <= args.digenome_fraction_cutoff
                    ):
                        caller_reasons.append("LOW_FORWARD_FRACTION")
                    if (
                        reverse_metrics["endpoint_fraction"]
                        <= args.digenome_fraction_cutoff
                    ):
                        caller_reasons.append("LOW_REVERSE_FRACTION")
                    if score <= args.digenome_score_cutoff:
                        caller_reasons.append("LOW_DIGENOME_SCORE")
                    pair_candidates.append(
                        DigenomePairCandidate(
                            contig=contig,
                            forward_position=forward_position,
                            reverse_position=reverse_position,
                            forward_metrics=forward_metrics,
                            reverse_metrics=reverse_metrics,
                            score=score,
                            caller_filter_reasons=caller_reasons,
                        )
                    )
        return pair_candidates

    def closed_interval_candidates(
        owned_interval: GenomicInterval,
    ) -> list[DigenomePairCandidate]:
        contig_length = bam.get_reference_length(owned_interval.contig)
        overhang = args.digenome_overhang
        pair_window = args.digenome_pair_window
        direct_reverse_start = (
            owned_interval.start - overhang - pair_window
        )
        direct_reverse_end = owned_interval.end - overhang + pair_window
        scan_interval = GenomicInterval(
            contig=owned_interval.contig,
            start=owned_interval.start,
            end=owned_interval.end,
            scan_start=max(
                0,
                min(owned_interval.scan_start, direct_reverse_start),
            ),
            scan_end=min(
                contig_length,
                max(owned_interval.scan_end, direct_reverse_end),
            ),
        )
        expansion_step = max(1, abs(overhang) + pair_window)

        for _iteration in range(64):
            candidates = build_pair_candidates(
                scan_intervals=[scan_interval],
            )
            relevant_components = [
                component
                for component in group_digenome_pair_components(candidates)
                if any(
                    owned_interval.owns(
                        candidate.contig,
                        candidate.forward_position,
                    )
                    for candidate in component
                )
            ]
            required_start = scan_interval.scan_start
            required_end = scan_interval.scan_end
            for component in relevant_components:
                forward_positions = {
                    candidate.forward_position
                    for candidate in component
                }
                reverse_positions = {
                    candidate.reverse_position
                    for candidate in component
                }
                required_start = min(
                    required_start,
                    *(
                        position - overhang - pair_window
                        for position in forward_positions
                    ),
                    *(
                        position + overhang - pair_window
                        for position in reverse_positions
                    ),
                )
                required_end = max(
                    required_end,
                    *(
                        position - overhang + pair_window + 1
                        for position in forward_positions
                    ),
                    *(
                        position + overhang + pair_window + 1
                        for position in reverse_positions
                    ),
                )

            target_start = max(0, required_start)
            target_end = min(contig_length, required_end)
            if (
                target_start >= scan_interval.scan_start
                and target_end <= scan_interval.scan_end
            ):
                return candidates

            new_scan_start = scan_interval.scan_start
            new_scan_end = scan_interval.scan_end
            if target_start < scan_interval.scan_start:
                new_scan_start = max(
                    0,
                    min(
                        target_start,
                        scan_interval.scan_start - expansion_step,
                    ),
                )
            if target_end > scan_interval.scan_end:
                new_scan_end = min(
                    contig_length,
                    max(
                        target_end,
                        scan_interval.scan_end + expansion_step,
                    ),
                )
            scan_interval = GenomicInterval(
                contig=owned_interval.contig,
                start=owned_interval.start,
                end=owned_interval.end,
                scan_start=new_scan_start,
                scan_end=new_scan_end,
            )
            expansion_step *= 2

        raise RuntimeError(
            "Digenome matching component expansion did not converge for "
            f"{owned_interval.contig}:{owned_interval.start}-"
            f"{owned_interval.end}"
        )

    pair_candidates: list[DigenomePairCandidate] = []
    selected_candidates: list[DigenomePairCandidate] = []
    if selected_intervals is None:
        pair_candidates = build_pair_candidates(selected_contigs=contigs)
        selected_candidates = select_digenome_pairs(pair_candidates)
    else:
        for owned_interval in selected_intervals:
            interval_candidates = closed_interval_candidates(owned_interval)
            pair_candidates.extend(
                candidate
                for candidate in interval_candidates
                if owned_interval.owns(
                    candidate.contig,
                    candidate.forward_position,
                )
            )
            selected_candidates.extend(
                candidate
                for candidate in select_digenome_pairs(interval_candidates)
                if owned_interval.owns(
                    candidate.contig,
                    candidate.forward_position,
                )
            )

    rows: list[dict[str, Any]] = []
    for candidate in selected_candidates:
        contig = candidate.contig
        forward_position = candidate.forward_position
        if not position_is_owned(
            selected_intervals,
            contig,
            forward_position,
            genome_blacklist,
        ):
            continue
        reverse_position = candidate.reverse_position
        forward_metrics = candidate.forward_metrics
        reverse_metrics = candidate.reverse_metrics
        combined = combine_site_metrics(forward_metrics, reverse_metrics)
        row = create_base_call_row(
            args,
            contig,
            forward_position,
            "both",
            combined,
        )
        row.update(
            {
                "forward_position_0based": forward_position,
                "forward_position_1based": forward_position + 1,
                "reverse_position_0based": reverse_position,
                "reverse_position_1based": reverse_position + 1,
                "opposite_position_0based": reverse_position,
                "opposite_position_1based": reverse_position + 1,
                "opposite_count": reverse_metrics["endpoint_count"],
                "opposite_depth": reverse_metrics["strand_depth"],
                "opposite_fraction": reverse_metrics["endpoint_fraction"],
                "forward_endpoint_count": forward_metrics["endpoint_count"],
                "forward_depth": forward_metrics["strand_depth"],
                "forward_fraction": forward_metrics["endpoint_fraction"],
                "reverse_endpoint_count": reverse_metrics["endpoint_count"],
                "reverse_depth": reverse_metrics["strand_depth"],
                "reverse_fraction": reverse_metrics["endpoint_fraction"],
                "digenome_score": candidate.score,
                "signal_classification": "DSB",
                "_caller_filter_reasons": (
                    candidate.caller_filter_reasons
                ),
            }
        )
        row["known_indel_overlap"] = ";".join(
            known_indels(
                variant_file,
                contig,
                [forward_position, reverse_position],
                args.artifact_window,
            )
        )
        rows.append(row)
    return rows, len(pair_candidates)


def log_combination(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact p-value without a SciPy runtime dependency."""
    row_one = a + b
    row_two = c + d
    column_one = a + c
    total = row_one + row_two
    if total == 0:
        return 1.0

    low = max(0, column_one - row_two)
    high = min(row_one, column_one)

    def probability(cell_a: int) -> float:
        return math.exp(
            log_combination(row_one, cell_a)
            + log_combination(row_two, column_one - cell_a)
            - log_combination(total, column_one)
        )

    observed = probability(a)
    return min(
        1.0,
        sum(
            probability(cell_a)
            for cell_a in range(low, high + 1)
            if probability(cell_a) <= observed * (1.0 + 1e-12)
        ),
    )


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    if not p_values:
        return []
    ordered = sorted(enumerate(p_values), key=lambda item: item[1])
    adjusted = [1.0] * len(p_values)
    running = 1.0
    total = len(p_values)
    for rank_from_end, (original_index, p_value) in enumerate(
        reversed(ordered), start=1
    ):
        rank = total - rank_from_end + 1
        running = min(running, p_value * total / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


# Input compatibility, controls, and final filtering.

def check_bam_layout(
    bam: Any, min_mapq: int, require_paired: bool, label: str
) -> None:
    checked = 0
    paired = 0
    for read in bam.fetch(until_eof=True):
        if not eligible_primary(read, min_mapq):
            continue
        checked += 1
        paired += int(read.is_paired)
        if checked >= 10000:
            break
    bam.reset()
    if checked == 0:
        raise ValueError(f"{label} contains no eligible primary alignments")
    if require_paired and paired == 0:
        raise ValueError(
            f"{label} requires paired-end alignments for nDigenome analysis"
        )


def validate_variant_contigs(
    variant_file: Any,
    bam: Any,
    analyzed_contigs: Iterable[str] | None,
) -> None:
    required_contigs = (
        list(analyzed_contigs)
        if analyzed_contigs is not None
        else list(bam.references)
    )
    variant_contigs = set(variant_file.header.contigs)
    missing_contigs = sorted(set(required_contigs) - variant_contigs)
    if missing_contigs:
        preview = ", ".join(missing_contigs[:10])
        remainder = len(missing_contigs) - 10
        if remainder > 0:
            preview += f", ... ({remainder} more)"
        raise ValueError(
            "Variant VCF is missing analyzed BAM contig(s): "
            f"{preview}. Check reference builds and contig naming."
        )

    bam_lengths = dict(zip(bam.references, bam.lengths))
    length_mismatches: list[str] = []
    for contig in required_contigs:
        variant_length = variant_file.header.contigs[contig].length
        if (
            variant_length is not None
            and variant_length != bam_lengths[contig]
        ):
            length_mismatches.append(
                f"{contig} (BAM {bam_lengths[contig]}, "
                f"VCF {variant_length})"
            )
    if length_mismatches:
        raise ValueError(
            "Variant VCF contig length does not match the BAM: "
            + ", ".join(length_mismatches[:10])
        )


def add_control_evidence(
    rows: list[dict[str, Any]],
    control_bam: Any | None,
    args: argparse.Namespace,
    adjust_q_values: bool = True,
) -> None:
    min_mapq = (
        args.ndigenome_min_mapq
        if args.analysis == "ndigenome"
        else args.digenome_min_mapq
    )
    controlled_indexes: list[int] = []
    for index, row in enumerate(rows):
        if control_bam is None:
            continue
        if args.analysis == "digenome":
            forward = measure_site_metrics(
                control_bam,
                row["contig"],
                row["forward_position_0based"],
                "+",
                args.artifact_window,
                min_mapq,
            )
            reverse = measure_site_metrics(
                control_bam,
                row["contig"],
                row["reverse_position_0based"],
                "-",
                args.artifact_window,
                min_mapq,
            )
            control = combine_site_metrics(forward, reverse)
            row.update(
                {
                    "control_forward_endpoint_count": forward["endpoint_count"],
                    "control_forward_depth": forward["strand_depth"],
                    "control_forward_fraction": forward["endpoint_fraction"],
                    "control_reverse_endpoint_count": reverse["endpoint_count"],
                    "control_reverse_depth": reverse["strand_depth"],
                    "control_reverse_fraction": reverse["endpoint_fraction"],
                    "control_digenome_score": digenome_score(
                        forward["endpoint_count"],
                        forward["endpoint_fraction"],
                        reverse["endpoint_count"],
                        reverse["endpoint_fraction"],
                    ),
                }
            )
        else:
            control = measure_site_metrics(
                control_bam,
                row["contig"],
                row["position_0based"],
                row["strand"],
                args.artifact_window,
                min_mapq,
            )
            if row["strand"] == "+":
                row.update(
                    {
                        "control_forward_endpoint_count": control[
                            "endpoint_count"
                        ],
                        "control_forward_depth": control["strand_depth"],
                        "control_forward_fraction": control[
                            "endpoint_fraction"
                        ],
                    }
                )
            else:
                row.update(
                    {
                        "control_reverse_endpoint_count": control[
                            "endpoint_count"
                        ],
                        "control_reverse_depth": control["strand_depth"],
                        "control_reverse_fraction": control[
                            "endpoint_fraction"
                        ],
                    }
                )

        row.update(
            {
                "control_sample": args.control_sample,
                "control_endpoint_count": control["endpoint_count"],
                "control_depth": control["strand_depth"],
                "control_fraction": control["endpoint_fraction"],
            }
        )

        if control["strand_depth"] < args.control_min_depth:
            row["control_status"] = "INSUFFICIENT_CONTROL_COVERAGE"
            continue

        treated_nonendpoint = max(
            0,
            row["strand_depth"] - row["endpoint_count"],
        )
        control_nonendpoint = max(
            0,
            control["strand_depth"] - control["endpoint_count"],
        )
        fisher_p = fisher_exact_two_sided(
            row["endpoint_count"],
            treated_nonendpoint,
            control["endpoint_count"],
            control_nonendpoint,
        )
        treated_rate = (row["endpoint_count"] + 0.5) / (
            row["strand_depth"] + 1.0
        )
        control_rate = (control["endpoint_count"] + 0.5) / (
            control["strand_depth"] + 1.0
        )
        row.update(
            {
                "control_status": "MATCHED_CONTROL",
                "control_fold_enrichment": treated_rate / control_rate,
                "control_fisher_p": fisher_p,
            }
        )
        controlled_indexes.append(index)

    if adjust_q_values:
        q_values = benjamini_hochberg(
            [rows[index]["control_fisher_p"] for index in controlled_indexes]
        )
        for index, q_value in zip(controlled_indexes, q_values):
            rows[index]["control_fisher_q"] = q_value


def apply_filters_to_row(
    row: dict[str, Any],
    args: argparse.Namespace,
) -> None:
    artifact_reasons: list[str] = []
    if row["softclip_fraction"] >= args.max_softclip_fraction:
        artifact_reasons.append("HIGH_5P_SOFTCLIP")
    if row["indel_fraction"] >= args.max_indel_fraction:
        artifact_reasons.append("NEARBY_INDEL")
    if row["known_indel_overlap"]:
        artifact_reasons.append("KNOWN_INDEL")
    if row["support_mean_mapq"] < args.min_support_mean_mapq:
        artifact_reasons.append("LOW_SUPPORT_MAPQ")
    if row["control_status"] == "INSUFFICIENT_CONTROL_COVERAGE":
        artifact_reasons.append("INSUFFICIENT_CONTROL_COVERAGE")
    if row["control_status"] == "MATCHED_CONTROL":
        if row["control_fraction"] > args.control_max_fraction:
            artifact_reasons.append("HIGH_CONTROL_FRACTION")
        if row["control_fold_enrichment"] < args.control_min_fold:
            artifact_reasons.append("LOW_CONTROL_FOLD")
        if row["control_fisher_q"] > args.control_max_q:
            artifact_reasons.append("CONTROL_Q_FAIL")

    row["classification"] = (
        "ARTIFACT_RISK"
        if artifact_reasons
        else row["signal_classification"]
    )
    filter_reasons = list(row.get("_caller_filter_reasons", []))
    filter_reasons.extend(artifact_reasons)
    if (
        args.analysis == "ndigenome"
        and row["signal_classification"] != "SSB"
    ):
        filter_reasons.append(row["signal_classification"])
    row["filter_status"] = "PASS" if not filter_reasons else "FILTERED"
    row["filter_reasons"] = ";".join(dict.fromkeys(filter_reasons))


def format_value(field: str, value: Any) -> str:
    if value is None:
        return ""
    if field in FLOAT_FIELDS:
        return f"{float(value):.8g}"
    return str(value)


# Final TSV, BED, QC, and MultiQC output.

def parameter_snapshot(args: argparse.Namespace) -> dict[str, Any]:
    internal = {
        "chunk_id",
        "contigs_file",
        "intervals_file",
        "raw_output",
        "summary_output",
    }
    return {
        key: value
        for key, value in vars(args).items()
        if key not in internal
    }


def write_output_stream(
    rows: Iterable[dict[str, Any]],
    candidate_count: int,
    args: argparse.Namespace,
) -> dict[str, Any]:
    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    stem = f"{output_prefix}.{args.analysis}"
    all_path = Path(f"{stem}.all.tsv")
    high_path = Path(f"{stem}.high_confidence.tsv")
    manual_review_path = Path(f"{stem}.manual_review.tsv")
    artifact_path = Path(f"{stem}.artifact.tsv")
    bed_path = Path(f"{stem}.bed")
    qc_path = Path(f"{stem}.qc.json")
    multiqc_path = Path(f"{output_prefix}.{args.analysis}_mqc.tsv")

    classification_counts: Counter[str] = Counter()
    signal_counts: Counter[str] = Counter()
    reported_count = 0
    high_confidence_count = 0
    manual_review_count = 0
    artifact_count = 0

    with (
        all_path.open("w", newline="") as all_handle,
        high_path.open("w", newline="") as high_handle,
        manual_review_path.open("w", newline="") as manual_review_handle,
        artifact_path.open("w", newline="") as artifact_handle,
        bed_path.open("w") as bed_handle,
    ):
        all_writer = csv.DictWriter(
            all_handle,
            fieldnames=TSV_FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        high_writer = csv.DictWriter(
            high_handle,
            fieldnames=TSV_FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        manual_review_writer = csv.DictWriter(
            manual_review_handle,
            fieldnames=TSV_FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        artifact_writer = csv.DictWriter(
            artifact_handle,
            fieldnames=TSV_FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        all_writer.writeheader()
        high_writer.writeheader()
        manual_review_writer.writeheader()
        artifact_writer.writeheader()

        for row in rows:
            reported_count += 1
            classification_counts[row["classification"]] += 1
            signal_counts[row["signal_classification"]] += 1
            formatted = {
                field: format_value(field, row.get(field))
                for field in TSV_FIELDS
            }
            all_writer.writerow(formatted)

            if row["filter_status"] == "PASS":
                high_confidence_count += 1
                high_writer.writerow(formatted)
                score = min(1000, round(row["endpoint_fraction"] * 1000))
                name = (
                    f"{args.sample}|{row['strand']}|"
                    f"{row['classification']}"
                )
                bed_handle.write(
                    f"{row['contig']}\t{row['position_0based']}\t"
                    f"{row['position_0based'] + 1}\t{name}\t{score}\t"
                    f"{'.' if row['strand'] == 'both' else row['strand']}\n"
                )
            elif row["classification"] == "ARTIFACT_RISK":
                artifact_count += 1
                artifact_writer.writerow(formatted)
            else:
                manual_review_count += 1
                manual_review_writer.writerow(formatted)

    qc = {
        "schema_version": 5,
        "sample": args.sample,
        "analysis": args.analysis,
        "control_status": (
            "MATCHED_CONTROL" if args.control_sample else "UNCONTROLLED"
        ),
        "coordinate_convention": {
            "internal": "0-based aligned 5-prime reference coordinate",
            "tsv": "both 0-based and 1-based columns",
            "bed": "0-based half-open",
        },
        "alignment_counting": {
            "keep_multimappers": args.keep_multimappers,
            "primary_alignments": (
                "counted once at the BWA-selected primary placement"
            ),
            "secondary_and_supplementary_alignments": "diagnostic only",
        },
        "genome_blacklist": {
            "enabled": bool(getattr(args, "genome_blacklist", "")),
            "path": getattr(args, "genome_blacklist", ""),
            "sha256": getattr(args, "genome_blacklist_sha256", ""),
            "intervals": getattr(
                args,
                "genome_blacklist_intervals",
                0,
            ),
            "excluded_bases": getattr(
                args,
                "genome_blacklist_excluded_bases",
                0,
            ),
        },
        "warnings": (
            [
                "Multimapper mode counts each read only at its selected "
                "primary placement; support can be diluted across equivalent "
                "reference copies."
            ]
            if args.keep_multimappers
            else []
        ),
        "parameters": parameter_snapshot(args),
        "candidate_endpoints_or_pairs_before_filters": candidate_count,
        "reported_candidates": reported_count,
        "high_confidence_calls": high_confidence_count,
        "manual_review_candidates": manual_review_count,
        "artifact_candidates": artifact_count,
        "high_confidence_ssb": (
            high_confidence_count if args.analysis == "ndigenome" else 0
        ),
        "high_confidence_dsb": (
            high_confidence_count if args.analysis == "digenome" else 0
        ),
        "classifications": dict(classification_counts),
        "signal_classifications": dict(signal_counts),
    }
    with qc_path.open("w") as handle:
        json.dump(qc, handle, indent=2, sort_keys=True)
        handle.write("\n")

    with multiqc_path.open("w") as handle:
        handle.write(f"# id: {args.analysis}_summary\n")
        handle.write(
            f"# section_name: {args.analysis} cleavage summary\n"
        )
        handle.write(
            "# description: Unified strand-aware cleavage endpoint calls\n"
        )
        handle.write("# plot_type: table\n")
        handle.write(
            "Sample\tReported candidates\tHigh-confidence calls\t"
            "Manual-review candidates\tArtifact candidates\t"
            "DSB\tSSB\tPossible DSB\tAmbiguous\n"
        )
        handle.write(
            f"{args.sample}\t{reported_count}\t{high_confidence_count}\t"
            f"{manual_review_count}\t{artifact_count}\t"
            f"{signal_counts.get('DSB', 0)}\t"
            f"{signal_counts.get('SSB', 0)}\t"
            f"{signal_counts.get('POSSIBLE_DSB', 0)}\t"
            f"{signal_counts.get('AMBIGUOUS', 0)}\n"
        )
    return qc


def write_outputs(
    rows: list[dict[str, Any]],
    candidate_count: int,
    args: argparse.Namespace,
) -> dict[str, Any]:
    rows.sort(
        key=lambda row: (
            row["contig"],
            row["position_0based"],
            row["strand"],
        )
    )
    return write_output_stream(rows, candidate_count, args)


# Command-line validation and serial/chunk execution.

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Call Digenome-seq DSB or nDigenome-seq SSB endpoint pileups."
        )
    )
    parser.add_argument("--analysis", choices=sorted(ANALYSES), required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--control-bam")
    parser.add_argument("--control-sample", default="")
    parser.add_argument("--variant-vcf")
    parser.add_argument("--genome-blacklist")
    parser.add_argument("--keep-multimappers", action="store_true")
    parser.add_argument("--intervals-file")
    parser.add_argument("--chunk-id")
    parser.add_argument("--raw-output")
    parser.add_argument("--summary-output")

    parser.add_argument("--artifact-window", type=int, default=10)
    parser.add_argument("--max-softclip-fraction", type=float, default=0.20)
    parser.add_argument("--max-indel-fraction", type=float, default=0.20)
    parser.add_argument("--min-support-mean-mapq", type=float, default=10.0)
    parser.add_argument("--control-min-depth", type=int, default=1)
    parser.add_argument("--control-max-fraction", type=float, default=0.05)
    parser.add_argument("--control-min-fold", type=float, default=5.0)
    parser.add_argument("--control-max-q", type=float, default=0.05)

    parser.add_argument("--ndigenome-min-count", type=int, default=10)
    parser.add_argument(
        "--ndigenome-min-fraction",
        type=float,
        default=0.20,
    )
    parser.add_argument("--ndigenome-min-mapq", type=int, default=1)
    parser.add_argument("--ndigenome-opposite-window", type=int, default=5)
    parser.add_argument(
        "--ndigenome-ambiguous-min-count",
        type=int,
        default=3,
    )
    parser.add_argument(
        "--ndigenome-ambiguous-min-fraction",
        type=float,
        default=0.05,
    )

    parser.add_argument("--digenome-overhang", type=int, default=0)
    parser.add_argument("--digenome-pair-window", type=int, default=2)
    parser.add_argument("--digenome-min-mapq", type=int, default=1)
    parser.add_argument("--digenome-forward-cutoff", type=int, default=5)
    parser.add_argument("--digenome-reverse-cutoff", type=int, default=5)
    parser.add_argument("--digenome-depth-cutoff", type=int, default=10)
    parser.add_argument(
        "--digenome-fraction-cutoff", type=float, default=0.20
    )
    parser.add_argument("--digenome-score-cutoff", type=float, default=2.5)
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if args.analysis not in ANALYSES:
        raise ValueError(
            f"Unknown analysis '{args.analysis}'. "
            f"Expected one of: {', '.join(sorted(ANALYSES))}"
        )
    chunk_fields = (
        getattr(args, "intervals_file", None),
        getattr(args, "chunk_id", None),
        getattr(args, "raw_output", None),
        getattr(args, "summary_output", None),
    )
    if any(chunk_fields) and not all(chunk_fields):
        raise ValueError(
            "Chunk mode requires --intervals-file, --chunk-id, "
            "--raw-output, and --summary-output"
        )
    if args.ndigenome_min_count < 1:
        raise ValueError("--ndigenome-min-count must be at least 1")
    if args.control_min_depth < 1:
        raise ValueError("--control-min-depth must be at least 1")
    for name in (
        "ndigenome_min_fraction",
        "ndigenome_ambiguous_min_fraction",
        "digenome_fraction_cutoff",
        "max_softclip_fraction",
        "max_indel_fraction",
        "control_max_fraction",
        "control_max_q",
    ):
        value = getattr(args, name)
        if not 0 <= value <= 1:
            raise ValueError(
                f"--{name.replace('_', '-')} must be between 0 and 1"
            )
    for name in (
        "artifact_window",
        "ndigenome_opposite_window",
        "ndigenome_ambiguous_min_count",
        "digenome_pair_window",
        "ndigenome_min_mapq",
        "digenome_min_mapq",
        ):
        if getattr(args, name) < 0:
            raise ValueError(
                f"--{name.replace('_', '-')} must be non-negative"
            )
    for name in (
        "digenome_forward_cutoff",
        "digenome_reverse_cutoff",
        "digenome_depth_cutoff",
    ):
        if getattr(args, name) < 0:
            raise ValueError(
                f"--{name.replace('_', '-')} must be non-negative"
            )
    for name in (
        "digenome_score_cutoff",
        "min_support_mean_mapq",
        "control_min_fold",
    ):
        if getattr(args, name) < 0:
            raise ValueError(
                f"--{name.replace('_', '-')} must be non-negative"
            )


def call_candidate_rows(
    args: argparse.Namespace,
    contigs: Iterable[str] | None = None,
    intervals: Iterable[GenomicInterval] | None = None,
    adjust_q_values: bool = True,
) -> tuple[list[dict[str, Any]], int]:
    if pysam is None:
        raise RuntimeError(
            "pysam is required. Run this command in the cleavage container."
        )
    validate_args(args)
    min_mapq = (
        args.ndigenome_min_mapq
        if args.analysis == "ndigenome"
        else args.digenome_min_mapq
    )
    require_paired = args.analysis == "ndigenome"
    if contigs is not None and intervals is not None:
        raise ValueError("Specify contigs or intervals, not both")
    selected_intervals = (
        list(intervals) if intervals is not None else None
    )
    analyzed_contigs = (
        list(dict.fromkeys(interval.contig for interval in selected_intervals))
        if selected_intervals is not None
        else contigs
    )

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        if not bam.has_index():
            raise ValueError(f"BAM is not indexed: {args.bam}")
        sort_order = bam.header.to_dict().get("HD", {}).get("SO")
        if sort_order != "coordinate":
            raise ValueError(
                "BAM must declare coordinate sort order, found: "
                f"{sort_order or 'missing'}"
            )
        if selected_intervals is not None:
            validate_intervals_against_bam(selected_intervals, bam)
        reference_lengths = dict(zip(bam.references, bam.lengths))
        genome_blacklist_path = getattr(
            args,
            "genome_blacklist",
            None,
        )
        genome_blacklist = (
            load_genome_blacklist(
                genome_blacklist_path,
                reference_lengths,
            )
            if genome_blacklist_path
            else None
        )
        args.genome_blacklist = genome_blacklist_path or ""
        args.genome_blacklist_sha256 = (
            genome_blacklist.sha256
            if genome_blacklist is not None
            else ""
        )
        args.genome_blacklist_intervals = (
            genome_blacklist.interval_count
            if genome_blacklist is not None
            else 0
        )
        args.genome_blacklist_excluded_bases = (
            genome_blacklist.excluded_bases
            if genome_blacklist is not None
            else 0
        )
        check_bam_layout(bam, min_mapq, require_paired, "BAM")

        control_bam = None
        variant_file = None
        try:
            if args.control_bam and Path(args.control_bam).stat().st_size > 0:
                control_bam = pysam.AlignmentFile(args.control_bam, "rb")
                if not control_bam.has_index():
                    raise ValueError(
                        f"Control BAM is not indexed: {args.control_bam}"
                    )
                control_sort = (
                    control_bam.header.to_dict().get("HD", {}).get("SO")
                )
                if control_sort != "coordinate":
                    raise ValueError(
                        "Control BAM must declare coordinate sort order, "
                        f"found: {control_sort or 'missing'}"
                    )
                if (
                    control_bam.references != bam.references
                    or control_bam.lengths != bam.lengths
                ):
                    raise ValueError(
                        "Control BAM reference names and lengths do not match "
                        "the treated BAM"
                    )
                check_bam_layout(
                    control_bam,
                    min_mapq,
                    require_paired,
                    "Control BAM",
                )

            if args.variant_vcf and Path(args.variant_vcf).stat().st_size > 0:
                if not any(
                    Path(f"{args.variant_vcf}{suffix}").is_file()
                    for suffix in (".tbi", ".csi")
                ):
                    raise ValueError(
                        "Variant VCF must have a .tbi or .csi index: "
                        f"{args.variant_vcf}"
                    )
                variant_file = pysam.VariantFile(args.variant_vcf)
                validate_variant_contigs(
                    variant_file,
                    bam,
                    analyzed_contigs,
                )

            if args.analysis == "ndigenome":
                rows, candidate_count = call_ndigenome(
                    bam,
                    variant_file,
                    args,
                    contigs,
                    selected_intervals,
                    genome_blacklist,
                )
            else:
                rows, candidate_count = call_digenome(
                    bam,
                    variant_file,
                    args,
                    contigs,
                    selected_intervals,
                    genome_blacklist,
                )
            add_control_evidence(
                rows,
                control_bam,
                args,
                adjust_q_values=adjust_q_values,
            )
        finally:
            if control_bam is not None:
                control_bam.close()
            if variant_file is not None:
                variant_file.close()
    return rows, candidate_count


def read_intervals_file(path: str) -> list[GenomicInterval]:
    intervals: list[GenomicInterval] = []
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_fields = {
            "contig",
            "owner_start",
            "owner_end",
            "scan_start",
            "scan_end",
        }
        if reader.fieldnames is None or not required_fields.issubset(
            reader.fieldnames
        ):
            raise ValueError(
                f"Interval file has an invalid header: {path}"
            )
        for line_number, row in enumerate(reader, start=2):
            try:
                interval = GenomicInterval(
                    contig=row["contig"],
                    start=int(row["owner_start"]),
                    end=int(row["owner_end"]),
                    scan_start=int(row["scan_start"]),
                    scan_end=int(row["scan_end"]),
                )
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid interval coordinates in {path} "
                    f"line {line_number}"
                ) from exc
            if not interval.contig:
                raise ValueError(
                    f"Missing interval contig in {path} line {line_number}"
                )
            if not (
                0 <= interval.scan_start
                <= interval.start
                < interval.end
                <= interval.scan_end
            ):
                raise ValueError(
                    f"Invalid owned/padded interval in {path} "
                    f"line {line_number}"
                )
            intervals.append(interval)
    if not intervals:
        raise ValueError(f"Interval file is empty: {path}")

    previous_by_contig: dict[str, GenomicInterval] = {}
    for interval in sorted(
        intervals,
        key=lambda item: (item.contig, item.start, item.end),
    ):
        previous = previous_by_contig.get(interval.contig)
        if previous is not None and interval.start < previous.end:
            raise ValueError(
                "Interval file contains overlapping ownership ranges on "
                f"{interval.contig}: {previous.start}-{previous.end} and "
                f"{interval.start}-{interval.end}"
            )
        previous_by_contig[interval.contig] = interval
    return intervals


def validate_intervals_against_bam(
    intervals: Iterable[GenomicInterval],
    bam: Any,
) -> None:
    lengths = dict(zip(bam.references, bam.lengths))
    for interval in intervals:
        if interval.contig not in lengths:
            raise ValueError(
                f"Requested contig is absent from the BAM: {interval.contig}"
            )
        contig_length = int(lengths[interval.contig])
        if interval.scan_end > contig_length:
            raise ValueError(
                f"Interval exceeds BAM contig length for {interval.contig}: "
                f"{interval.scan_end} > {contig_length}"
            )


def mapped_bam_contig_lengths(bam_path: str) -> dict[str, int]:
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        if not bam.has_index():
            raise ValueError(f"BAM is not indexed: {bam_path}")
        lengths = dict(zip(bam.references, bam.lengths))
        return {
            index_stats.contig: int(lengths[index_stats.contig])
            for index_stats in bam.get_index_statistics()
            if index_stats.mapped > 0
        }


def run_chunk_calling(args: argparse.Namespace) -> dict[str, Any]:
    validate_args(args)
    intervals = read_intervals_file(args.intervals_file)
    contigs = list(dict.fromkeys(interval.contig for interval in intervals))
    rows, candidate_count = call_candidate_rows(
        args,
        intervals=intervals,
        adjust_q_values=False,
    )
    rows.sort(
        key=lambda row: (
            row["contig"],
            row["position_0based"],
            row["strand"],
        )
    )

    raw_path = Path(args.raw_output)
    summary_path = Path(args.summary_output)
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(raw_path, "wt", compresslevel=1) as handle:
        for row in rows:
            row["_chunk_id"] = args.chunk_id
            handle.write(json.dumps(row, sort_keys=True))
            handle.write("\n")

    expected_contig_lengths = mapped_bam_contig_lengths(args.bam)
    summary = {
        "schema_version": 2,
        "sample": args.sample,
        "analysis": args.analysis,
        "chunk_id": args.chunk_id,
        "contigs": contigs,
        "intervals": [
            {
                "contig": interval.contig,
                "start": interval.start,
                "end": interval.end,
                "scan_start": interval.scan_start,
                "scan_end": interval.scan_end,
            }
            for interval in intervals
        ],
        "expected_contigs": list(expected_contig_lengths),
        "expected_contig_lengths": expected_contig_lengths,
        "candidate_count": candidate_count,
        "reported_rows": len(rows),
        "parameters": parameter_snapshot(args),
    }
    with summary_path.open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return summary


def run_serial_calling(args: argparse.Namespace) -> dict[str, Any]:
    rows, candidate_count = call_candidate_rows(args)
    for row in rows:
        apply_filters_to_row(row, args)

    return write_outputs(rows, candidate_count, args)


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    try:
        if args.intervals_file:
            chunk = run_chunk_calling(args)
            print(
                f"{args.analysis} chunk {args.chunk_id}: "
                f"{chunk['reported_rows']} row(s) from "
                f"{len(chunk['intervals'])} interval(s)."
            )
            return
        qc = run_serial_calling(args)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
    print(
        f"{args.analysis}: reported {qc['reported_candidates']} candidate(s), "
        f"{qc['high_confidence_calls']} high-confidence call(s)."
    )


if __name__ == "__main__":
    main()
