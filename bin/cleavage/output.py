"""Output columns, filters, and the files written for each sample.

Every evaluated row goes to <sample>.<analysis>.all.tsv. Rows that pass
every filter also go to .high_confidence.tsv and the BED; filtered rows with
artifact evidence go to .artifact.tsv; the remaining filtered rows go to
.manual_review.tsv. The `tier` column names that file.
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from typing import Iterable

from .bam import SiteMetrics
from .settings import CallerSettings
from .stats import INSUFFICIENT_CONTROL, MATCHED_CONTROL

ARTIFACT_COLUMNS = [
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
]

DIGENOME_COLUMNS = [
    "sample", "analysis", "contig",
    "forward_position_0based", "forward_position_1based",
    "reverse_position_0based", "reverse_position_1based",
    "tier", "filter_status", "filter_reasons",
    "digenome_pair_score", "rgen_digenome_score",
    "forward_endpoint_count", "forward_depth", "forward_fraction",
    "reverse_endpoint_count", "reverse_depth", "reverse_fraction",
    "combined_endpoint_count", "combined_depth", "combined_fraction",
    *ARTIFACT_COLUMNS,
    "control_sample", "control_status",
    "control_forward_endpoint_count", "control_forward_depth", "control_forward_fraction",
    "control_reverse_endpoint_count", "control_reverse_depth", "control_reverse_fraction",
    "control_combined_endpoint_count", "control_combined_depth", "control_combined_fraction",
    "control_digenome_pair_score", "control_rgen_digenome_score",
    "control_fold_enrichment", "control_fisher_p", "control_fisher_q",
]

NDIGENOME_COLUMNS = [
    "sample", "analysis", "contig", "position_0based", "position_1based", "strand",
    "signal_classification", "tier", "filter_status", "filter_reasons",
    "endpoint_count", "strand_depth", "endpoint_fraction",
    "opposite_position_0based", "opposite_position_1based",
    "opposite_count", "opposite_depth", "opposite_fraction",
    *ARTIFACT_COLUMNS,
    "control_sample", "control_status",
    "control_endpoint_count", "control_depth", "control_fraction",
    "control_fold_enrichment", "control_fisher_p", "control_fisher_q",
]

TIERS = ["high_confidence", "manual_review", "artifact"]


def one_based(position: int | None) -> int | None:
    return None if position is None else position + 1


def artifact_columns(metrics: SiteMetrics, known_indels: list[str]) -> dict:
    return {
        "support_mean_mapq": metrics.support_mean_mapq,
        "local_mean_mapq": metrics.local_mean_mapq,
        "support_mean_nm": metrics.support_mean_nm,
        "local_mean_nm": metrics.local_mean_nm,
        "softclip_fraction": metrics.softclip_fraction,
        "secondary_endpoint_count": metrics.secondary_endpoint_count,
        "indel_position_0based": metrics.indel_position,
        "indel_position_1based": one_based(metrics.indel_position),
        "indel_type": metrics.indel_type,
        "indel_length": metrics.indel_length,
        "indel_read_count": metrics.indel_read_count,
        "indel_fraction": metrics.indel_fraction,
        "known_indel_overlap": ";".join(known_indels),
    }


# Filters.


def artifact_reasons(row: dict, settings: CallerSettings) -> list[str]:
    """Shared artifact filters. Comparisons follow AGENTS.md: clipping and
    indel fractions fail at >=, support MAPQ below, control fraction above,
    fold below, and q above their limits."""
    reasons = []
    if row["softclip_fraction"] >= settings.cleavage_max_softclip_fraction:
        reasons.append("HIGH_5P_SOFTCLIP")
    if row["indel_fraction"] >= settings.cleavage_max_indel_fraction:
        reasons.append("NEARBY_INDEL")
    if row["known_indel_overlap"]:
        reasons.append("KNOWN_INDEL")
    if row["support_mean_mapq"] < settings.cleavage_min_support_mean_mapq:
        reasons.append("LOW_SUPPORT_MAPQ")
    if row["control_status"] == INSUFFICIENT_CONTROL:
        reasons.append("INSUFFICIENT_CONTROL_COVERAGE")
    if row["control_status"] == MATCHED_CONTROL:
        if settings.analysis == "digenome":
            control_fraction = row["control_combined_fraction"]
        else:
            control_fraction = row["control_fraction"]
        if control_fraction > settings.cleavage_control_max_fraction:
            reasons.append("HIGH_CONTROL_FRACTION")
        if row["control_fold_enrichment"] < settings.cleavage_control_min_fold:
            reasons.append("LOW_CONTROL_FOLD")
        if row["control_fisher_q"] > settings.cleavage_control_max_q:
            reasons.append("CONTROL_Q_FAIL")
    return reasons


def apply_filters(row: dict, settings: CallerSettings) -> None:
    """Set filter_status, filter_reasons, and tier.

    Digenome rows carry their caller threshold failures in
    `caller_filter_reasons`. nDigenome rows that are not SSB are filtered with
    their signal class as the reason.
    """
    artifacts = artifact_reasons(row, settings)
    reasons = list(row.get("caller_filter_reasons", [])) + artifacts
    if settings.analysis == "ndigenome" and row["signal_classification"] != "SSB":
        reasons.append(row["signal_classification"])
    reasons = list(dict.fromkeys(reasons))
    row["filter_status"] = "FILTERED" if reasons else "PASS"
    row["filter_reasons"] = ";".join(reasons)
    if not reasons:
        row["tier"] = "high_confidence"
    elif artifacts:
        row["tier"] = "artifact"
    else:
        row["tier"] = "manual_review"


# Files.


def format_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.8g}"
    return str(value)


def bed_line(row: dict, settings: CallerSettings) -> str:
    """One BED6 record. Digenome calls sit at the forward endpoint with no strand."""
    if settings.analysis == "digenome":
        start, strand, fraction = row["forward_position_0based"], ".", row["combined_fraction"]
        name = f"{row['sample']}|both|DSB"
    else:
        start, strand, fraction = row["position_0based"], row["strand"], row["endpoint_fraction"]
        name = f"{row['sample']}|{strand}|{row['signal_classification']}"
    score = min(1000, round(fraction * 1000))
    return f"{row['contig']}\t{start}\t{start + 1}\t{name}\t{score}\t{strand}\n"


def write_outputs(
    rows: Iterable[dict],
    prefix: str,
    settings: CallerSettings,
    report: dict,
) -> dict:
    """Write the TSVs, BED, QC JSON, and MultiQC table; return the QC report.

    `rows` must already be filtered and in output order. `report` holds the
    QC fields known before writing (sample, settings, chunks, candidates, ...).
    """
    stem = f"{prefix}.{settings.analysis}"
    columns = DIGENOME_COLUMNS if settings.analysis == "digenome" else NDIGENOME_COLUMNS
    tiers: Counter[str] = Counter({tier: 0 for tier in TIERS})
    signals: Counter[str] = Counter()
    handles = {name: open(f"{stem}.{name}.tsv", "w", newline="") for name in ["all", *TIERS]}
    try:
        writers = {
            name: csv.writer(handle, delimiter="\t", lineterminator="\n")
            for name, handle in handles.items()
        }
        for writer in writers.values():
            writer.writerow(columns)
        with open(f"{stem}.bed", "w") as bed:
            for row in rows:
                values = [format_value(row.get(column)) for column in columns]
                writers["all"].writerow(values)
                writers[row["tier"]].writerow(values)
                tiers[row["tier"]] += 1
                if settings.analysis == "ndigenome":
                    signals[row["signal_classification"]] += 1
                if row["tier"] == "high_confidence":
                    bed.write(bed_line(row, settings))
    finally:
        for handle in handles.values():
            handle.close()

    qc = {**report, "rows": sum(tiers.values()), "tiers": dict(tiers)}
    if settings.analysis == "ndigenome":
        qc["signal_classifications"] = dict(sorted(signals.items()))
    Path(f"{stem}.qc.json").write_text(json.dumps(qc, indent=2) + "\n")
    write_multiqc_table(f"{prefix}.{settings.analysis}_mqc.tsv", qc, settings)
    return qc


def write_multiqc_table(path: str, qc: dict, settings: CallerSettings) -> None:
    header = ["Sample", "Reported candidates", "High-confidence calls", "Manual-review candidates", "Artifact candidates"]
    values = [qc["sample"], qc["rows"], *(qc["tiers"][tier] for tier in TIERS)]
    if settings.analysis == "ndigenome":
        header += ["SSB", "Possible DSB", "Ambiguous"]
        signals = qc["signal_classifications"]
        values += [signals.get(name, 0) for name in ("SSB", "POSSIBLE_DSB", "AMBIGUOUS")]
    lines = [
        f"# id: {settings.analysis}_summary",
        f"# section_name: {settings.analysis} cleavage summary",
        "# description: Cleavage calls per sample",
        "# plot_type: table",
        "\t".join(header),
        "\t".join(str(value) for value in values),
    ]
    Path(path).write_text("\n".join(lines) + "\n")
