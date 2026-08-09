#!/usr/bin/env python3
"""Validate and normalize the Digenome-seq pipeline samplesheet."""

from __future__ import annotations

import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

REQUIRED_COLUMNS = ["sample", "fastq_1", "fastq_2"]
OPTIONAL_COLUMNS = ["lane", "control", "variant_vcf"]
OUTPUT_COLUMNS = [
    "sample",
    "lane",
    "fastq_1",
    "fastq_2",
    "control",
    "variant_vcf",
    "variant_index",
    "is_control",
]
SAMPLE_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
FASTQ_RE = re.compile(r"\.(fastq|fq)\.gz$", re.IGNORECASE)
VCF_RE = re.compile(r"\.vcf\.gz$", re.IGNORECASE)
VALID_ANALYSES = {"digenome", "ndigenome"}
ALLOWED_COLUMNS = set(REQUIRED_COLUMNS + OPTIONAL_COLUMNS)


def die(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def clean(value: object) -> str:
    return str(value or "").strip()


def resolve_existing(path_text: str) -> str:
    return str(Path(path_text).expanduser().resolve())


def find_vcf_index(vcf_path: str) -> str:
    for suffix in (".tbi", ".csi"):
        candidate = Path(vcf_path + suffix)
        if candidate.exists():
            return str(candidate.resolve())
    return ""


def validate_samplesheet(input_csv: Path, output_csv: Path, analysis: str) -> None:
    if analysis not in VALID_ANALYSES:
        die(
            f"Unknown analysis '{analysis}'. "
            f"Expected one of: {', '.join(sorted(VALID_ANALYSES))}"
        )
    if not input_csv.exists():
        die(f"Samplesheet does not exist: {input_csv}")

    rows: list[dict[str, str]] = []
    seen_sample_lanes: set[tuple[str, str]] = set()
    sample_layout: dict[str, str] = {}
    sample_metadata: dict[str, tuple[str, str, str]] = {}
    auto_lane_counts: dict[str, int] = defaultdict(int)
    staged_basenames: dict[str, dict[str, str]] = defaultdict(dict)
    seen_fastq_paths: dict[str, tuple[int, str, str]] = {}
    errors: list[str] = []

    with input_csv.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            die("Samplesheet is empty or missing a header row.")

        normalized_headers = [clean(name) for name in reader.fieldnames]
        if len(set(normalized_headers)) != len(normalized_headers):
            die("Samplesheet contains duplicate column names.")
        reader.fieldnames = normalized_headers
        has_lane_column = "lane" in reader.fieldnames
        unknown_columns = sorted(set(reader.fieldnames) - ALLOWED_COLUMNS)
        if unknown_columns:
            die(
                "Samplesheet contains unknown column(s): "
                + ", ".join(unknown_columns)
                + ". Allowed columns are: "
                + ", ".join(REQUIRED_COLUMNS + OPTIONAL_COLUMNS)
            )

        missing = [col for col in REQUIRED_COLUMNS if col not in reader.fieldnames]
        if missing:
            die(
                "Samplesheet is missing required column(s): "
                + ", ".join(missing)
                + ". Expected columns: "
                + ", ".join(REQUIRED_COLUMNS)
            )

        for line_number, row in enumerate(reader, start=2):
            sample = clean(row.get("sample"))
            fastq_1 = clean(row.get("fastq_1"))
            fastq_2 = clean(row.get("fastq_2"))
            control = clean(row.get("control"))
            variant_vcf = clean(row.get("variant_vcf"))

            if not sample:
                errors.append(f"line {line_number}: sample is blank")
                continue
            if has_lane_column:
                lane = clean(row.get("lane"))
            else:
                auto_lane_counts[sample] += 1
                lane = f"L{auto_lane_counts[sample]:03d}"
            if not SAMPLE_RE.fullmatch(sample):
                errors.append(
                    f"line {line_number}: sample '{sample}' contains unsupported characters. "
                    "Use only letters, numbers, dots, underscores, and hyphens."
                )
            if control and not SAMPLE_RE.fullmatch(control):
                errors.append(
                    f"line {line_number}: control '{control}' contains "
                    "unsupported characters"
                )
            if not lane:
                errors.append(
                    f"line {line_number}: lane is blank. Either provide a value "
                    "for every row or omit the lane column to auto-number rows."
                )
            if not fastq_1:
                errors.append(f"line {line_number}: fastq_1 is blank")

            key = (sample, lane)
            if key in seen_sample_lanes:
                errors.append(
                    f"line {line_number}: duplicate sample/lane combination: {sample}/{lane}"
                )
            seen_sample_lanes.add(key)

            layout = "PE" if fastq_2 else "SE"
            if sample in sample_layout and sample_layout[sample] != layout:
                errors.append(
                    f"line {line_number}: sample '{sample}' mixes single-end and paired-end rows"
                )
            sample_layout[sample] = layout
            if analysis == "ndigenome" and layout != "PE":
                errors.append(
                    f"line {line_number}: nDigenome analysis currently requires paired-end data"
                )

            resolved_fastqs: dict[str, str] = {}
            for label, fq in (("fastq_1", fastq_1), ("fastq_2", fastq_2)):
                if not fq:
                    resolved_fastqs[label] = ""
                    continue
                if not FASTQ_RE.search(fq):
                    errors.append(
                        f"line {line_number}: {label} should end with .fastq.gz or .fq.gz: {fq}"
                    )
                if not os.path.exists(fq):
                    errors.append(f"line {line_number}: {label} file does not exist: {fq}")
                resolved = resolve_existing(fq)
                resolved_fastqs[label] = resolved
                previous_use = seen_fastq_paths.get(resolved)
                if previous_use:
                    previous_line, previous_label, previous_sample = previous_use
                    errors.append(
                        f"line {line_number}: {label} reuses FASTQ '{resolved}' "
                        f"already used as {previous_label} for sample "
                        f"'{previous_sample}' on line {previous_line}"
                    )
                else:
                    seen_fastq_paths[resolved] = (
                        line_number,
                        label,
                        sample,
                    )
                basename = Path(resolved).name
                previous = staged_basenames[sample].get(basename)
                if previous and previous != resolved:
                    errors.append(
                        f"line {line_number}: sample '{sample}' has FASTQs from different "
                        f"directories with the same staged basename '{basename}'"
                    )
                staged_basenames[sample][basename] = resolved

            resolved_vcf = ""
            variant_index = ""
            if variant_vcf:
                if not VCF_RE.search(variant_vcf):
                    errors.append(
                        f"line {line_number}: variant_vcf must be a bgzip-compressed .vcf.gz file"
                    )
                if not os.path.exists(variant_vcf):
                    errors.append(
                        f"line {line_number}: variant_vcf file does not exist: {variant_vcf}"
                    )
                else:
                    resolved_vcf = resolve_existing(variant_vcf)
                    variant_index = find_vcf_index(resolved_vcf)
                    if not variant_index:
                        errors.append(
                            f"line {line_number}: variant_vcf is not indexed with .tbi or .csi: "
                            f"{resolved_vcf}"
                        )

            metadata = (control, resolved_vcf, variant_index)
            if sample in sample_metadata and sample_metadata[sample] != metadata:
                errors.append(
                    f"line {line_number}: sample '{sample}' has inconsistent control "
                    "or variant_vcf values across lanes"
                )
            sample_metadata[sample] = metadata

            if sample and lane and fastq_1:
                rows.append(
                    {
                        "sample": sample,
                        "lane": lane,
                        "fastq_1": resolved_fastqs.get("fastq_1", ""),
                        "fastq_2": resolved_fastqs.get("fastq_2", ""),
                        "control": control,
                        "variant_vcf": resolved_vcf,
                        "variant_index": variant_index,
                        "is_control": "false",
                    }
                )

    if not rows:
        errors.append("Samplesheet contains no usable data rows.")

    referenced_controls = {
        control
        for control, _vcf, _index in sample_metadata.values()
        if control
    }
    for sample, (control, _vcf, _index) in sample_metadata.items():
        if not control:
            continue
        if control == sample:
            errors.append(
                f"sample '{sample}' cannot use itself as its control"
            )
        elif control not in sample_metadata:
            errors.append(
                f"sample '{sample}' references control '{control}', but no sample "
                f"named '{control}' exists in the samplesheet"
            )

    for control_sample in sorted(referenced_controls):
        if control_sample not in sample_metadata:
            continue
        nested_control = sample_metadata[control_sample][0]
        if nested_control:
            errors.append(
                f"sample '{control_sample}' is referenced as a control but also "
                f"declares control '{nested_control}'. Control rows must leave "
                "the control column blank."
            )

    if errors:
        print(
            "Samplesheet validation failed with the following problem(s):",
            file=sys.stderr,
        )
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
        sys.exit(1)

    for row in rows:
        row["is_control"] = (
            "true" if row["sample"] in referenced_controls else "false"
        )

    rows.sort(
        key=lambda row: (
            0 if row["is_control"] == "true" else 1,
            row["sample"],
            row["lane"],
            row["fastq_1"],
            row["fastq_2"],
        )
    )

    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    sample_counts = defaultdict(int)
    for row in rows:
        sample_counts[row["sample"]] += 1
    print(
        f"Validated {len(rows)} FASTQ row(s) across {len(sample_counts)} "
        f"sample(s) for {analysis} analysis."
    )


def main() -> None:
    if len(sys.argv) not in (3, 4):
        die(
            "Usage: validate_samplesheet.py "
            "<input.csv> <output.valid.csv> [digenome|ndigenome]"
        )
    analysis = sys.argv[3].lower() if len(sys.argv) == 4 else "digenome"
    validate_samplesheet(Path(sys.argv[1]), Path(sys.argv[2]), analysis)


if __name__ == "__main__":
    main()
