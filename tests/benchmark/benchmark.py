#!/usr/bin/env python3
"""Run the old (b9455e2) and new cleavage callers on one BAM and compare them.

Both run their chunk steps one after another with the default settings, so
the timings are comparable step by step. Prints wall time, CPU time, and peak
memory for each step, then checks that every output value of the new caller
matches the old one after mapping old columns to the new layout.

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
from dataclasses import fields
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "tests"))
from test_golden import expected_bed, expected_rows  # noqa: E402  (also puts bin/ on sys.path)

import pysam  # noqa: E402
from cleavage.settings import CallerSettings  # noqa: E402

TIERS = ["all", "high_confidence", "manual_review", "artifact"]
# ru_maxrss is in KiB on Linux and in bytes on macOS.
RSS_UNITS_PER_MB = 1024 * 1024 if sys.platform == "darwin" else 1024


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
    optional = []
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


def compare(args, legacy_prefix: str, new_prefix: str) -> list[str]:
    """Names of the output files whose values differ."""
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        order = {name: index for index, name in enumerate(bam.references)}
    differences = []
    for tier in TIERS:
        expected = expected_rows(Path(f"{legacy_prefix}.{args.analysis}.{tier}.tsv"), args.analysis, order)
        with open(f"{new_prefix}.{args.analysis}.{tier}.tsv", newline="") as handle:
            observed = list(csv.reader(handle, delimiter="\t"))[1:]
        if observed != expected:
            differences.append(f"{tier}.tsv")
    bed = Path(f"{new_prefix}.{args.analysis}.bed").read_text().splitlines()
    if bed != expected_bed(Path(f"{legacy_prefix}.{args.analysis}.bed"), order):
        differences.append("bed")
    return differences


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
    print("outputs match" if not differences else f"OUTPUTS DIFFER: {', '.join(differences)}")
    sys.exit(1 if differences else 0)


if __name__ == "__main__":
    main()
