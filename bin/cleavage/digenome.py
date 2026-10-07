"""Digenome calls: forward and reverse endpoints paired into double-strand breaks.

A forward endpoint f and a reverse endpoint r can pair when
|r - (f - overhang)| <= pair_window. Each endpoint joins at most one pair; the
pairs are chosen by a deterministic min-cost matching that maximizes, in
order: pairs passing every caller threshold, the number of pairs, the total
pair score, and stable coordinate order.

Work is split in two:
- Each chunk measures the endpoints it owns that have at least one possible
  partner (`endpoint_records`).
- The finalizer sees every endpoint of a contig at once, builds all candidate
  pairs, and runs the matching (`pair_rows`). No chain of competing pairs can
  be cut by a chunk boundary, so chunked and serial runs agree.

Count, depth, fraction, and pair-score cutoffs are strict `>`.
"""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from dataclasses import asdict, dataclass, replace

import pysam

from .bam import (
    SiteMetrics,
    combine,
    count_endpoints,
    measure_site,
    ratio,
    reads_covering_either,
    rgen_counts,
    rgen_score,
)
from .output import artifact_columns
from .regions import Blacklist, OwnedInterval, callable_segments, merge_segments
from .settings import CallerSettings
from .stats import UNCONTROLLED, compare_to_control
from .vcf import known_indels


def pair_score(forward: SiteMetrics, reverse: SiteMetrics) -> float:
    """digenome_pair_score: drives pairing and filtering (unlike the RGEN score)."""
    return (
        forward.endpoint_fraction
        * reverse.endpoint_fraction
        * (forward.endpoint_count + reverse.endpoint_count)
        / 4.0
    )


# Chunk side: measure endpoints.


def endpoint_records(
    bam: pysam.AlignmentFile,
    control_bam: pysam.AlignmentFile | None,
    vcf: pysam.VariantFile | None,
    settings: CallerSettings,
    intervals: list[OwnedInterval],
    blacklist: Blacklist | None,
) -> list[dict]:
    """One record per owned, non-blacklisted endpoint with enough reads and at
    least one possible partner. The scan reaches |overhang| + pair_window past
    the owned range only to find partners owned by other chunks."""
    scan_minimum = min(settings.digenome_forward_cutoff, settings.digenome_reverse_cutoff) + 1
    reach = abs(settings.digenome_overhang) + settings.digenome_pair_window
    rgen_caches: tuple[dict, dict] = ({}, {})
    records = []
    for contig in dict.fromkeys(interval.contig for interval in intervals):
        owned = [interval for interval in intervals if interval.contig == contig]
        length = bam.get_reference_length(contig)
        scan_ranges = merge_segments(
            [(max(0, interval.start - reach), min(length, interval.end + reach)) for interval in owned]
        )
        endpoints: dict[tuple[int, str], int] = {}
        for start, end in scan_ranges:
            for segment_start, segment_end in callable_segments(contig, start, end, blacklist):
                endpoints.update(
                    count_endpoints(bam, contig, segment_start, segment_end, settings, scan_minimum)
                )
        forward = sorted(position for position, strand in endpoints if strand == "+")
        reverse = sorted(position for position, strand in endpoints if strand == "-")
        for position, strand in sorted(endpoints):
            if not any(interval.owns(contig, position) for interval in owned):
                continue
            partners = partner_positions(position, strand, forward, reverse, settings)
            if not partners:
                continue
            records.append(
                endpoint_record(
                    bam, control_bam, vcf, settings, blacklist, contig, position, strand, partners, rgen_caches
                )
            )
    return records


def partner_positions(
    position: int, strand: str, forward: list[int], reverse: list[int], settings: CallerSettings
) -> list[int]:
    """Opposite-strand endpoints within the pair window of this one."""
    overhang, window = settings.digenome_overhang, settings.digenome_pair_window
    if strand == "+":
        low, high, others = position - overhang - window, position - overhang + window, reverse
    else:
        low, high, others = position + overhang - window, position + overhang + window, forward
    return others[bisect_left(others, low):bisect_right(others, high)]


def endpoint_record(
    bam: pysam.AlignmentFile,
    control_bam: pysam.AlignmentFile | None,
    vcf: pysam.VariantFile | None,
    settings: CallerSettings,
    blacklist: Blacklist | None,
    contig: str,
    position: int,
    strand: str,
    partners: list[int],
    rgen_caches: tuple[dict, dict],
) -> dict:
    record = {
        "contig": contig,
        "position_0based": position,
        "strand": strand,
        "metrics": asdict(measure_site(bam, contig, position, strand, settings)),
        "control_metrics": None,
        "known_indels": known_indels(vcf, contig, [position], settings.cleavage_artifact_window),
        "rgen_score": None,
        "control_rgen_score": None,
        "pair_depths": {},
    }
    if control_bam is not None:
        record["control_metrics"] = asdict(measure_site(control_bam, contig, position, strand, settings))
    # A long read spanning the cut covers both endpoints of a pair, so the pair's
    # combined depth is measured here, per partner, with each read counted once.
    if settings.long_reads and strand == "+":
        record["pair_depths"] = {
            str(partner): [
                reads_covering_either(bam, contig, position, partner, settings.min_mapq),
                None if control_bam is None
                else reads_covering_either(control_bam, contig, position, partner, settings.min_mapq),
            ]
            for partner in partners
        }
    # The RGEN score depends only on the forward endpoint, so it is computed
    # here rather than for each candidate pair.
    if strand == "+":
        record["rgen_score"] = rgen_score_at(bam, contig, position, settings, blacklist, rgen_caches[0])
        if control_bam is not None:
            record["control_rgen_score"] = rgen_score_at(
                control_bam, contig, position, settings, blacklist, rgen_caches[1]
            )
    return record


def rgen_score_at(
    bam: pysam.AlignmentFile,
    contig: str,
    forward_position: int,
    settings: CallerSettings,
    blacklist: Blacklist | None,
    cache: dict,
) -> float:
    length = bam.get_reference_length(contig)

    def counts_at(position: int):
        if not 0 <= position < length:
            return None
        if blacklist is not None and blacklist.contains(contig, position):
            return None
        if (contig, position) not in cache:
            cache[contig, position] = rgen_counts(bam, contig, position, settings)
        return cache[contig, position]

    return rgen_score(counts_at, forward_position, settings.digenome_overhang)


# Finalizer side: pair endpoints.


@dataclass
class PairCandidate:
    contig: str
    forward_position: int
    reverse_position: int
    forward_metrics: SiteMetrics
    reverse_metrics: SiteMetrics
    pair_score: float
    caller_filter_reasons: list[str]

    @property
    def passes_caller_thresholds(self) -> bool:
        return not self.caller_filter_reasons


def pair_rows(
    records: list[dict],
    settings: CallerSettings,
    sample: str,
    control_sample: str,
) -> tuple[list[dict], int]:
    """Pair one contig's endpoint records. Returns the rows for the selected
    pairs, in forward-position order, and the number of candidate pairs."""
    forward = {record["position_0based"]: record for record in records if record["strand"] == "+"}
    reverse = {record["position_0based"]: record for record in records if record["strand"] == "-"}
    if len(forward) + len(reverse) != len(records):
        raise ValueError(f"Duplicate Digenome endpoint records on {records[0]['contig']}")
    candidates = candidate_pairs(forward, reverse, settings)
    rows = [
        pair_row(candidate, forward[candidate.forward_position], reverse[candidate.reverse_position],
                 settings, sample, control_sample)
        for candidate in select_pairs(candidates)
    ]
    return rows, len(candidates)


def candidate_pairs(forward: dict[int, dict], reverse: dict[int, dict], settings: CallerSettings) -> list[PairCandidate]:
    reverse_positions = sorted(reverse)
    candidates = []
    for forward_position in sorted(forward):
        expected = forward_position - settings.digenome_overhang
        low = bisect_left(reverse_positions, expected - settings.digenome_pair_window)
        high = bisect_right(reverse_positions, expected + settings.digenome_pair_window)
        forward_metrics = SiteMetrics(**forward[forward_position]["metrics"])
        for reverse_position in reverse_positions[low:high]:
            reverse_metrics = SiteMetrics(**reverse[reverse_position]["metrics"])
            score = pair_score(forward_metrics, reverse_metrics)
            candidates.append(
                PairCandidate(
                    contig=forward[forward_position]["contig"],
                    forward_position=forward_position,
                    reverse_position=reverse_position,
                    forward_metrics=forward_metrics,
                    reverse_metrics=reverse_metrics,
                    pair_score=score,
                    caller_filter_reasons=caller_filter_reasons(forward_metrics, reverse_metrics, score, settings),
                )
            )
    return candidates


def caller_filter_reasons(
    forward: SiteMetrics, reverse: SiteMetrics, score: float, settings: CallerSettings
) -> list[str]:
    """Digenome thresholds; each must be strictly exceeded."""
    checks = [
        (forward.endpoint_count, settings.digenome_forward_cutoff, "LOW_FORWARD_COUNT"),
        (reverse.endpoint_count, settings.digenome_reverse_cutoff, "LOW_REVERSE_COUNT"),
        (forward.strand_depth, settings.digenome_depth_cutoff, "LOW_FORWARD_DEPTH"),
        (reverse.strand_depth, settings.digenome_depth_cutoff, "LOW_REVERSE_DEPTH"),
        (forward.endpoint_fraction, settings.digenome_fraction_cutoff, "LOW_FORWARD_FRACTION"),
        (reverse.endpoint_fraction, settings.digenome_fraction_cutoff, "LOW_REVERSE_FRACTION"),
        (score, settings.digenome_pair_score_cutoff, "LOW_DIGENOME_PAIR_SCORE"),
    ]
    return [reason for value, cutoff, reason in checks if value <= cutoff]


def pair_row(
    candidate: PairCandidate,
    forward_record: dict,
    reverse_record: dict,
    settings: CallerSettings,
    sample: str,
    control_sample: str,
) -> dict:
    forward, reverse = candidate.forward_metrics, candidate.reverse_metrics
    combined = combine(forward, reverse)
    control_depth = None
    if settings.long_reads:
        treated_depth, control_depth = forward_record["pair_depths"][str(candidate.reverse_position)]
        combined = counted_once(combined, treated_depth)
    known = sorted(set(forward_record["known_indels"]) | set(reverse_record["known_indels"]))
    return {
        "sample": sample,
        "analysis": "digenome",
        "contig": candidate.contig,
        "forward_position_0based": candidate.forward_position,
        "forward_position_1based": candidate.forward_position + 1,
        "reverse_position_0based": candidate.reverse_position,
        "reverse_position_1based": candidate.reverse_position + 1,
        "digenome_pair_score": candidate.pair_score,
        "rgen_digenome_score": forward_record["rgen_score"],
        "forward_endpoint_count": forward.endpoint_count,
        "forward_depth": forward.strand_depth,
        "forward_fraction": forward.endpoint_fraction,
        "reverse_endpoint_count": reverse.endpoint_count,
        "reverse_depth": reverse.strand_depth,
        "reverse_fraction": reverse.endpoint_fraction,
        "combined_endpoint_count": combined.endpoint_count,
        "combined_depth": combined.strand_depth,
        "combined_fraction": combined.endpoint_fraction,
        **artifact_columns(combined, known),
        **pair_control_columns(forward_record, reverse_record, combined, control_depth, settings, control_sample),
        "caller_filter_reasons": candidate.caller_filter_reasons,
    }


def counted_once(combined: SiteMetrics, depth: int) -> SiteMetrics:
    """Pooled long-read metrics with the combined depth replaced by `depth`, the
    reads covering either endpoint each counted once. Summing the two sides'
    depths would count every spanning read twice."""
    return replace(combined, strand_depth=depth, endpoint_fraction=ratio(combined.endpoint_count, depth))


def pair_control_columns(
    forward_record: dict,
    reverse_record: dict,
    treated: SiteMetrics,
    control_depth: int | None,
    settings: CallerSettings,
    control_sample: str,
) -> dict:
    """Both endpoints measured in the matched control, if there is one.
    `control_depth` is the long-read combined depth (see counted_once)."""
    if forward_record["control_metrics"] is None:
        return {
            "control_sample": "",
            "control_status": UNCONTROLLED,
            "control_forward_endpoint_count": 0,
            "control_forward_depth": 0,
            "control_forward_fraction": 0.0,
            "control_reverse_endpoint_count": 0,
            "control_reverse_depth": 0,
            "control_reverse_fraction": 0.0,
            "control_combined_endpoint_count": 0,
            "control_combined_depth": 0,
            "control_combined_fraction": 0.0,
            "control_digenome_pair_score": None,
            "control_rgen_digenome_score": None,
            "control_fold_enrichment": None,
            "control_fisher_p": None,
            "control_fisher_q": None,
        }
    forward = SiteMetrics(**forward_record["control_metrics"])
    reverse = SiteMetrics(**reverse_record["control_metrics"])
    control = combine(forward, reverse)
    if control_depth is not None:
        control = counted_once(control, control_depth)
    status, fold, p_value = compare_to_control(
        treated.endpoint_count, treated.strand_depth,
        control.endpoint_count, control.strand_depth,
        settings.cleavage_control_min_depth,
    )
    return {
        "control_sample": control_sample,
        "control_status": status,
        "control_forward_endpoint_count": forward.endpoint_count,
        "control_forward_depth": forward.strand_depth,
        "control_forward_fraction": forward.endpoint_fraction,
        "control_reverse_endpoint_count": reverse.endpoint_count,
        "control_reverse_depth": reverse.strand_depth,
        "control_reverse_fraction": reverse.endpoint_fraction,
        "control_combined_endpoint_count": control.endpoint_count,
        "control_combined_depth": control.strand_depth,
        "control_combined_fraction": control.endpoint_fraction,
        "control_digenome_pair_score": pair_score(forward, reverse),
        "control_rgen_digenome_score": forward_record["control_rgen_score"],
        "control_fold_enrichment": fold,
        "control_fisher_p": p_value,
        "control_fisher_q": None,
    }


# One-to-one matching.


@dataclass
class MatchingEdge:
    destination: int
    reverse_edge_index: int
    capacity: int
    cost: tuple[int, int, int, int]


def connected_components(candidates: list[PairCandidate]) -> list[list[PairCandidate]]:
    """Group candidates that share an endpoint, directly or through a chain.
    Each component is matched on its own."""
    by_contig: dict[str, list[PairCandidate]] = {}
    for candidate in candidates:
        by_contig.setdefault(candidate.contig, []).append(candidate)

    components: list[list[PairCandidate]] = []
    for contig in sorted(by_contig):
        ordered = sorted(
            by_contig[contig],
            key=lambda candidate: (candidate.forward_position, candidate.reverse_position),
        )
        sharing_forward: dict[int, list[int]] = {}
        sharing_reverse: dict[int, list[int]] = {}
        for index, candidate in enumerate(ordered):
            sharing_forward.setdefault(candidate.forward_position, []).append(index)
            sharing_reverse.setdefault(candidate.reverse_position, []).append(index)

        unvisited = set(range(len(ordered)))
        while unvisited:
            pending = [min(unvisited)]
            members: list[int] = []
            while pending:
                index = pending.pop()
                if index not in unvisited:
                    continue
                unvisited.remove(index)
                members.append(index)
                candidate = ordered[index]
                neighbors = (
                    sharing_forward[candidate.forward_position]
                    + sharing_reverse[candidate.reverse_position]
                )
                pending.extend(
                    neighbor for neighbor in reversed(neighbors) if neighbor in unvisited
                )
            components.append([ordered[index] for index in sorted(members)])
    return components


def select_pairs(candidates: list[PairCandidate]) -> list[PairCandidate]:
    """Choose a deterministic best one-to-one set of pairs.

    Each component is solved as min-cost flow: source -> forward endpoints ->
    reverse endpoints -> sink, one unit per edge, with lexicographic costs
    (fails a threshold, -1 per pair, -score in exact integer units, coordinate
    rank). Augmenting along shortest paths (Bellman-Ford) until no path
    lowers the cost gives the optimum.
    """
    selected: list[PairCandidate] = []
    for component in connected_components(candidates):
        selected.extend(match_component(component))
    return sorted(
        selected,
        key=lambda candidate: (candidate.contig, candidate.forward_position, candidate.reverse_position),
    )


def match_component(component: list[PairCandidate]) -> list[PairCandidate]:
    forward_positions = sorted({candidate.forward_position for candidate in component})
    reverse_positions = sorted({candidate.reverse_position for candidate in component})
    source = 0
    first_forward = 1
    first_reverse = first_forward + len(forward_positions)
    sink = first_reverse + len(reverse_positions)
    graph: list[list[MatchingEdge]] = [[] for _node in range(sink + 1)]

    def add_edge(start: int, destination: int, cost: tuple[int, int, int, int]) -> MatchingEdge:
        forward_edge = MatchingEdge(destination, len(graph[destination]), 1, cost)
        backward_edge = MatchingEdge(start, len(graph[start]), 0, tuple(-value for value in cost))
        graph[start].append(forward_edge)
        graph[destination].append(backward_edge)
        return forward_edge

    zero_cost = (0, 0, 0, 0)
    forward_nodes = {position: first_forward + index for index, position in enumerate(forward_positions)}
    reverse_nodes = {position: first_reverse + index for index, position in enumerate(reverse_positions)}
    for node in forward_nodes.values():
        add_edge(source, node, zero_cost)
    for node in reverse_nodes.values():
        add_edge(node, sink, zero_cost)

    # Scores are binary floats; their denominators are powers of two, so the
    # largest one is a common denominator and the scores become exact integers.
    score_ratios = [float(candidate.pair_score).as_integer_ratio() for candidate in component]
    common_denominator = max(denominator for _numerator, denominator in score_ratios)
    candidate_edges = []
    for coordinate_rank, (candidate, (numerator, denominator)) in enumerate(zip(component, score_ratios)):
        score_units = numerator * (common_denominator // denominator)
        cost = (-int(candidate.passes_caller_thresholds), -1, -score_units, coordinate_rank)
        edge = add_edge(forward_nodes[candidate.forward_position], reverse_nodes[candidate.reverse_position], cost)
        candidate_edges.append((edge, candidate))

    while augment_shortest_path(graph, source, sink, zero_cost):
        pass
    return [candidate for edge, candidate in candidate_edges if edge.capacity == 0]


def augment_shortest_path(graph: list[list[MatchingEdge]], source: int, sink: int, zero_cost) -> bool:
    """Push one unit along the cheapest source-to-sink path if it lowers the
    total cost. Returns False when no such path exists."""
    distances: list[tuple[int, int, int, int] | None] = [None] * len(graph)
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
                distance = tuple(left + right for left, right in zip(distances[node], edge.cost))
                if distances[edge.destination] is None or distance < distances[edge.destination]:
                    distances[edge.destination] = distance
                    previous[edge.destination] = (node, edge_index)
                    changed = True
        if not changed:
            break

    if distances[sink] is None or distances[sink] >= zero_cost:
        return False
    node = sink
    visited: set[int] = set()
    while node != source:
        if node in visited:
            raise RuntimeError("Internal Digenome matching path contains a cycle")
        visited.add(node)
        step = previous[node]
        if step is None:
            raise RuntimeError("Internal Digenome matching path is incomplete")
        previous_node, edge_index = step
        edge = graph[previous_node][edge_index]
        edge.capacity = 0
        graph[edge.destination][edge.reverse_edge_index].capacity = 1
        node = previous_node
    return True
