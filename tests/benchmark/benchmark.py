#!/usr/bin/env python3
"""Run the old (b9455e2) and new cleavage callers on one BAM and compare them.

Both run their chunk steps one after another with the default settings, so
the timings are comparable step by step. Prints wall time, CPU time, and peak
memory for each step, then checks that every output value of the new caller
matches the old one after mapping old columns to the new layout. For
Digenome it also reports how often the pair-score cutoff agrees with the
standalone RGEN cutoff of 2.5.

See tests/benchmark/README.md for how to run it on Longleaf.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from collections import Counter
from dataclasses import fields
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import pysam  # noqa: E402
from cleavage.output import DIGENOME_COLUMNS, NDIGENOME_COLUMNS  # noqa: E402
from cleavage.settings import CallerSettings  # noqa: E402

TIERS = ["all", "high_confidence", "manual_review", "artifact"]
# New Digenome column -> old column, where the refactor renamed it.
RENAMED = {
    "combined_endpoint_count": "endpoint_count",
    "combined_depth": "strand_depth",
    "combined_fraction": "endpoint_fraction",
    "control_combined_endpoint_count": "control_endpoint_count",
    "control_combined_depth": "control_depth",
    "control_combined_fraction": "control_fraction",
}
# ru_maxrss is in KiB on Linux and in bytes on macOS.
RSS_UNITS_PER_MB = 1024 * 1024 if sys.platform == "darwin" else 1024
# The standalone RGEN tool's score cutoff (digenome -s 2.5), and the
# pair-score cutoffs tried against it: 0.50, 0.55, ..., 3.00.
RGEN_SCORE_CUTOFF = 2.5
CUTOFF_GRID = [step / 20 for step in range(10, 61)]
# Caller thresholds other than the pair score.
OTHER_CALLER_REASONS = {
    "LOW_FORWARD_COUNT", "LOW_REVERSE_COUNT", "LOW_FORWARD_DEPTH",
    "LOW_REVERSE_DEPTH", "LOW_FORWARD_FRACTION", "LOW_REVERSE_FRACTION",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--legacy-bin", required=True, help="bin/ extracted from commit b9455e2")
    parser.add_argument("--analysis", required=True, choices=["digenome", "ndigenome"])
    parser.add_argument("--bam", required=True)
    parser.add_argument("--control-bam")
    parser.add_argument("--vcf")
    parser.add_argument("--blacklist")
    parser.add_argument("--keep-multimappers", action="store_true")
    parser.add_argument("--chunks", type=int, default=8)
    parser.add_argument("--out", required=True, help="working directory for both runs")
    return parser.parse_args()


def run(label: str, command: list[str], cwd: Path, timings: list[dict], env: dict | None = None) -> None:
    """Run one step and record its wall time, CPU time, and peak memory."""
    started = time.monotonic()
    process = subprocess.Popen(command, cwd=cwd, env=env)
    _pid, status, usage = os.wait4(process.pid, 0)
    if status != 0:
        raise SystemExit(f"{label} failed: {' '.join(command)}")
    timings.append({
        "step": label,
        "wall_seconds": round(time.monotonic() - started, 1),
        "cpu_seconds": round(usage.ru_utime + usage.ru_stime, 1),
        "peak_rss_mb": round(usage.ru_maxrss / RSS_UNITS_PER_MB),
    })


def default_settings(analysis: str, keep_multimappers: bool) -> dict:
    """The pipeline's defaults from nextflow_schema.json, with the same
    --keep_multimappers overrides as callerSettings() in main.nf."""
    schema = json.loads((ROOT / "nextflow_schema.json").read_text())
    defaults = {
        name: definition.get("default")
        for group in schema["$defs"].values()
        for name, definition in group["properties"].items()
    }
    settings = {field.name: defaults.get(field.name) for field in fields(CallerSettings)}
    settings.update(analysis=analysis, keep_multimappers=keep_multimappers)
    if keep_multimappers:
        settings.update(digenome_min_mapq=0, ndigenome_min_mapq=0, cleavage_min_support_mean_mapq=0)
    return settings


def run_legacy(args, directory: Path, timings: list[dict]) -> str:
    legacy = Path(args.legacy_bin).resolve()
    # The old caller defaulted to RGEN's 2.5; compare at the pipeline's default.
    cutoff = default_settings(args.analysis, args.keep_multimappers)["digenome_pair_score_cutoff"]
    optional = ["--digenome-pair-score-cutoff", str(cutoff)]
    if args.control_bam:
        optional += ["--control-bam", args.control_bam, "--control-sample", "Control"]
    if args.vcf:
        optional += ["--variant-vcf", args.vcf]
    if args.blacklist:
        optional += ["--genome-blacklist", args.blacklist]
    if args.keep_multimappers:
        optional += ["--keep-multimappers", "--digenome-min-mapq", "0", "--ndigenome-min-mapq", "0",
                     "--min-support-mean-mapq", "0"]
    blacklist = ["--genome-blacklist", args.blacklist] if args.blacklist else []
    run("legacy plan", [sys.executable, str(legacy / "plan_cleavage_chunks.py"), "--bam", args.bam,
                        "--chunks", str(args.chunks), "--padding", "10", *blacklist,
                        "--output-dir", "plan", "--plan", "plan.tsv"], directory, timings)
    raw, summaries = [], []
    for intervals in sorted((directory / "plan").glob("*.intervals.tsv")):
        chunk = intervals.name.removesuffix(".intervals.tsv")
        raw.append(f"{chunk}.raw.jsonl.gz")
        summaries.append(f"{chunk}.chunk.json")
        run(f"legacy {chunk}", [sys.executable, str(legacy / "call_cleavage.py"), "--analysis", args.analysis,
                                "--bam", args.bam, "--sample", "Sample", "--output-prefix", "Sample", *optional,
                                "--intervals-file", str(intervals), "--chunk-id", chunk,
                                "--raw-output", raw[-1], "--summary-output", summaries[-1]], directory, timings)
    run("legacy finalize", [sys.executable, str(legacy / "finalize_cleavage_chunks.py"), "--analysis", args.analysis,
                            "--sample", "Sample", "--output-prefix", "Sample",
                            *[f"--raw-fragment={path}" for path in raw],
                            *[f"--chunk-summary={path}" for path in summaries]], directory, timings)
    return str(directory / "Sample")


def run_new(args, directory: Path, timings: list[dict]) -> str:
    settings = directory / "settings.json"
    settings.write_text(json.dumps(default_settings(args.analysis, args.keep_multimappers)))
    env = {**os.environ, "PYTHONPATH": str(ROOT / "bin")}
    optional = []
    if args.control_bam:
        optional += ["--control-bam", args.control_bam, "--control-sample", "Control"]
    if args.vcf:
        optional += ["--vcf", args.vcf]
    if args.blacklist:
        optional += ["--blacklist", args.blacklist]
    summaries = []
    for index in range(args.chunks):
        prefix = f"chunk_{index:03d}"
        summaries.append(f"{prefix}.json")
        run(f"new chunk {index}", [sys.executable, "-m", "cleavage", "call", "--settings", str(settings),
                                   "--bam", args.bam, "--sample", "Sample", "--chunk", str(index),
                                   "--chunks", str(args.chunks), "--out-prefix", prefix, *optional],
            directory, timings, env)
    control = ["--control-sample", "Control"] if args.control_bam else []
    run("new finalize", [sys.executable, "-m", "cleavage", "finalize", "--settings", str(settings),
                         "--sample", "Sample", *control, "--out-prefix", "Sample", *summaries],
        directory, timings, env)
    return str(directory / "Sample")


def legacy_tier(row: dict[str, str]) -> str:
    if row["filter_status"] == "PASS":
        return "high_confidence"
    return "artifact" if row["classification"] == "ARTIFACT_RISK" else "manual_review"


def legacy_rows_in_new_layout(path: Path, analysis: str, contig_order: dict[str, int]) -> list[list[str]]:
    """Old rows mapped to the new columns and sorted in the new order."""
    columns = DIGENOME_COLUMNS if analysis == "digenome" else NDIGENOME_COLUMNS
    position = "forward_position_0based" if analysis == "digenome" else "position_0based"
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rows.sort(key=lambda row: (contig_order[row["contig"]], int(row[position]), row["strand"]))
    mapped = []
    for row in rows:
        row["tier"] = legacy_tier(row)
        mapped.append([row[RENAMED.get(column, column)] for column in columns])
    return mapped


def legacy_bed_in_new_order(path: Path, contig_order: dict[str, int]) -> list[str]:
    lines = path.read_text().splitlines()
    return sorted(lines, key=lambda line: (contig_order[line.split("\t")[0]], int(line.split("\t")[1])))


def compare(args, legacy_prefix: str, new_prefix: str) -> list[str]:
    """Names of the output files whose values differ."""
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        order = {name: index for index, name in enumerate(bam.references)}
    differences = []
    for tier in TIERS:
        expected = legacy_rows_in_new_layout(Path(f"{legacy_prefix}.{args.analysis}.{tier}.tsv"), args.analysis, order)
        with open(f"{new_prefix}.{args.analysis}.{tier}.tsv", newline="") as handle:
            observed = list(csv.reader(handle, delimiter="\t"))[1:]
        if observed != expected:
            differences.append(f"{tier}.tsv")
    bed = Path(f"{new_prefix}.{args.analysis}.bed").read_text().splitlines()
    if bed != legacy_bed_in_new_order(Path(f"{legacy_prefix}.{args.analysis}.bed"), order):
        differences.append("bed")
    return differences


def score_agreement(prefix: str, cutoff: float) -> str:
    """How often pair score > `cutoff` agrees with RGEN score > 2.5, over the
    pairs that pass every other caller threshold ("Cutoff provenance" in
    docs/cleavage_algorithm.md)."""
    with open(f"{prefix}.digenome.all.tsv", newline="") as handle:
        pairs = [
            (float(row["digenome_pair_score"]), float(row["rgen_digenome_score"]) > RGEN_SCORE_CUTOFF)
            for row in csv.DictReader(handle, delimiter="\t")
            if not OTHER_CALLER_REASONS & set(row["filter_reasons"].split(";"))
        ]
    if not pairs:
        return "No pairs pass the count, depth, and fraction cutoffs; no score agreement to report."

    def agreeing(threshold: float) -> int:
        return sum((pair > threshold) == rgen_passes for pair, rgen_passes in pairs)

    table = Counter((pair > cutoff, rgen_passes) for pair, rgen_passes in pairs)
    current = agreeing(cutoff)
    agreement = {threshold: agreeing(threshold) for threshold in CUTOFF_GRID}
    best = max(agreement.values())
    best_cutoffs = [threshold for threshold, count in agreement.items() if count == best]
    return "\n".join([
        f"Pair score versus RGEN score, {len(pairs)} pairs passing the count, depth, and fraction cutoffs:",
        f"{'':<16}{f'RGEN > {RGEN_SCORE_CUTOFF}':>12}{f'RGEN <= {RGEN_SCORE_CUTOFF}':>13}",
        f"{f'pair > {cutoff}':<16}{table[True, True]:>12}{table[True, False]:>13}",
        f"{f'pair <= {cutoff}':<16}{table[False, True]:>12}{table[False, False]:>13}",
        f"agreement at {cutoff}: {current} of {len(pairs)} ({current / len(pairs):.1%})",
        f"best agreement over cutoffs {CUTOFF_GRID[0]} to {CUTOFF_GRID[-1]}: {best} of {len(pairs)} "
        f"({best / len(pairs):.1%}), at {best_cutoffs[0]} to {best_cutoffs[-1]}",
    ])


def main() -> None:
    args = parse_args()
    out = Path(args.out).resolve()
    timings: list[dict] = []
    legacy_dir, new_dir = out / "legacy", out / "new"
    legacy_dir.mkdir(parents=True, exist_ok=True)
    new_dir.mkdir(parents=True, exist_ok=True)
    legacy_prefix = run_legacy(args, legacy_dir, timings)
    new_prefix = run_new(args, new_dir, timings)
    differences = compare(args, legacy_prefix, new_prefix)

    print(f"{'step':<24}{'wall s':>10}{'cpu s':>10}{'peak MB':>10}")
    for row in timings:
        print(f"{row['step']:<24}{row['wall_seconds']:>10}{row['cpu_seconds']:>10}{row['peak_rss_mb']:>10}")
    for implementation in ("legacy", "new"):
        steps = [row for row in timings if row["step"].startswith(implementation)]
        print(
            f"{implementation} total: {sum(row['wall_seconds'] for row in steps):.1f} s wall, "
            f"max {max(row['peak_rss_mb'] for row in steps)} MB"
        )
    (out / "timings.json").write_text(json.dumps(timings, indent=2) + "\n")
    if args.analysis == "digenome":
        cutoff = default_settings(args.analysis, args.keep_multimappers)["digenome_pair_score_cutoff"]
        print(score_agreement(new_prefix, cutoff))
    print("outputs match" if not differences else f"OUTPUTS DIFFER: {', '.join(differences)}")
    sys.exit(1 if differences else 0)


if __name__ == "__main__":
    main()
