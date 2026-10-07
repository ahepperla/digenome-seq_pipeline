"""Validate and normalize the Digenome-seq pipeline samplesheet."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

REQUIRED_COLUMNS = ["sample", "fastq_1", "fastq_2"]
OPTIONAL_COLUMNS = ["control", "variant_vcf"]
SAMPLE_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
FASTQ_RE = re.compile(r"\.(fastq|fq)\.gz$", re.IGNORECASE)
VCF_RE = re.compile(r"\.vcf\.gz$", re.IGNORECASE)
VALID_ANALYSES = {"digenome", "ndigenome"}
ALLOWED_COLUMNS = set(REQUIRED_COLUMNS + OPTIONAL_COLUMNS)


def _clean(value: object) -> str:
    """Strip whitespace from a value and convert to string."""
    return str(value or "").strip()


def _resolve_existing(path_text: str) -> str:
    """Expand user home and resolve to absolute path."""
    return str(Path(path_text).expanduser().resolve())


def _find_vcf_index(vcf_path: str) -> str:
    """Find VCF index file (.tbi or .csi), preferring .tbi."""
    for suffix in (".tbi", ".csi"):
        candidate = Path(vcf_path + suffix)
        if candidate.exists():
            return str(candidate)
    return ""


def validate_samplesheet(input_csv: Path, analysis: str) -> list[dict]:
    """
    Return one record per sample, or raise ValueError listing every problem.

    Args:
        input_csv: Path to input samplesheet CSV.
        analysis: "digenome" or "ndigenome"; anything else raises ValueError.

    Returns:
        List of sample records (one per sample, sorted by sample name).
        Each record contains: sample, single_end, fastq_1, fastq_2, control,
        variant_vcf, variant_index, is_control.

    Raises:
        ValueError: If input CSV is invalid or any validation fails.
    """
    if analysis not in VALID_ANALYSES:
        raise ValueError(
            f"Unknown analysis '{analysis}'. "
            f"Expected one of: {', '.join(sorted(VALID_ANALYSES))}"
        )

    # Read and validate header
    if not input_csv.exists():
        raise ValueError(f"Samplesheet does not exist: {input_csv}")

    with input_csv.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("Samplesheet is empty or missing a header row.")

        # Normalize header names
        normalized_headers = [_clean(name) for name in reader.fieldnames]

        # Check for duplicate column names
        if len(set(normalized_headers)) != len(normalized_headers):
            raise ValueError("Samplesheet contains duplicate column names.")

        reader.fieldnames = normalized_headers

        # Check for unknown columns
        unknown_columns = sorted(set(reader.fieldnames) - ALLOWED_COLUMNS)
        if unknown_columns:
            raise ValueError(
                f"Samplesheet contains unknown column(s): "
                f"{', '.join(unknown_columns)}"
            )

        # Check for missing required columns
        missing = [col for col in REQUIRED_COLUMNS if col not in reader.fieldnames]
        if missing:
            raise ValueError(
                f"Samplesheet is missing required column(s): "
                f"{', '.join(missing)}"
            )

        # Collect all problems before raising
        problems: list[str] = []
        sample_rows: dict[str, list[dict]] = {}  # sample -> list of row data
        sample_layout: dict[str, str] = {}  # sample -> "SE" or "PE"
        sample_metadata: dict[str, tuple[str, str, str]] = {}  # sample -> (control, vcf, index)
        seen_fastq_paths: dict[str, tuple[int, str, str]] = {}  # resolved path -> (line, label, sample)

        for line_number, row in enumerate(reader, start=2):
            sample = _clean(row.get("sample", ""))
            fastq_1 = _clean(row.get("fastq_1", ""))
            fastq_2 = _clean(row.get("fastq_2", ""))
            control = _clean(row.get("control", ""))
            variant_vcf = _clean(row.get("variant_vcf", ""))

            # Check for blank sample (skip rest of row if blank)
            if not sample:
                problems.append(f"line {line_number}: sample is blank")
                continue

            # Check sample name format
            if not SAMPLE_RE.fullmatch(sample):
                problems.append(
                    f"line {line_number}: sample '{sample}' contains unsupported "
                    "characters. Use only letters, numbers, dots, underscores, and hyphens."
                )

            # Check control name format
            if control and not SAMPLE_RE.fullmatch(control):
                problems.append(
                    f"line {line_number}: control '{control}' contains unsupported characters"
                )

            # Check fastq_1 is not blank
            if not fastq_1:
                problems.append(f"line {line_number}: fastq_1 is blank")
                continue

            # Check layout (PE vs SE)
            layout = "PE" if fastq_2 else "SE"
            if sample in sample_layout and sample_layout[sample] != layout:
                problems.append(
                    f"line {line_number}: sample '{sample}' mixes single-end and paired-end rows"
                )
            sample_layout[sample] = layout

            # nDigenome requires paired-end
            if analysis == "ndigenome" and layout == "SE":
                problems.append(
                    f"line {line_number}: nDigenome requires paired-end data"
                )

            # Validate FASTQ files
            resolved_fastqs: dict[str, str] = {}
            for label, fq in (("fastq_1", fastq_1), ("fastq_2", fastq_2)):
                if not fq:
                    resolved_fastqs[label] = ""
                    continue

                # Check FASTQ suffix
                if not FASTQ_RE.search(fq):
                    problems.append(
                        f"line {line_number}: {label} must end with .fastq.gz or .fq.gz"
                    )

                # Check file exists
                if not Path(fq).expanduser().exists():
                    problems.append(
                        f"line {line_number}: {label} file does not exist: {fq}"
                    )

                resolved = _resolve_existing(fq)
                resolved_fastqs[label] = resolved

                # Check for FASTQ reuse across rows
                previous_use = seen_fastq_paths.get(resolved)
                if previous_use:
                    previous_line, previous_label, previous_sample = previous_use
                    problems.append(
                        f"line {line_number}: {label} reuses FASTQ already used as "
                        f"{previous_label} for sample '{previous_sample}' on line {previous_line}"
                    )
                else:
                    seen_fastq_paths[resolved] = (line_number, label, sample)

            # Check fastq_1 and fastq_2 are not identical
            if (
                resolved_fastqs.get("fastq_1")
                and resolved_fastqs.get("fastq_2")
                and resolved_fastqs["fastq_1"] == resolved_fastqs["fastq_2"]
            ):
                problems.append(
                    f"line {line_number}: fastq_2 reuses FASTQ from fastq_1"
                )

            # Handle VCF validation
            resolved_vcf = ""
            variant_index = ""
            if variant_vcf:
                # Check VCF suffix
                if not VCF_RE.search(variant_vcf):
                    problems.append(
                        f"line {line_number}: variant_vcf must be a bgzip-compressed .vcf.gz file"
                    )

                # Check file exists
                if not Path(variant_vcf).expanduser().exists():
                    problems.append(
                        f"line {line_number}: variant_vcf file does not exist"
                    )
                else:
                    resolved_vcf = _resolve_existing(variant_vcf)
                    variant_index = _find_vcf_index(resolved_vcf)
                    if not variant_index:
                        problems.append(
                            f"line {line_number}: variant_vcf is not indexed with .tbi or .csi"
                        )

            # Store metadata and check consistency across lanes
            metadata = (control, resolved_vcf, variant_index)
            if sample in sample_metadata and sample_metadata[sample] != metadata:
                problems.append(
                    f"line {line_number}: inconsistent control or variant_vcf values across lanes"
                )
            sample_metadata[sample] = metadata

            # Only add valid rows
            if sample and fastq_1:
                if sample not in sample_rows:
                    sample_rows[sample] = []
                sample_rows[sample].append(
                    {
                        "fastq_1": resolved_fastqs.get("fastq_1", ""),
                        "fastq_2": resolved_fastqs.get("fastq_2", ""),
                    }
                )

    # Check for no usable rows
    if not sample_rows:
        problems.append("Samplesheet contains no usable data rows.")

    # Cross-sample validation
    referenced_controls = {
        control
        for control, _vcf, _index in sample_metadata.values()
        if control
    }

    # Check for self-control and missing controls
    for sample, (control, _vcf, _index) in sample_metadata.items():
        if not control:
            continue
        if control == sample:
            problems.append(f"sample '{sample}' cannot use itself as its control")
        elif control not in sample_metadata:
            problems.append(
                f"sample '{sample}' references control '{control}', but no sample "
                f"named '{control}' exists in the samplesheet"
            )

    # Check for control chains (control sample must not have its own control)
    for control_sample in referenced_controls:
        if control_sample not in sample_metadata:
            continue
        nested_control = sample_metadata[control_sample][0]
        if nested_control:
            problems.append(
                f"Control rows must leave the control column blank"
            )

    if problems:
        raise ValueError("Samplesheet validation failed:\n" + "\n".join(f"  - {p}" for p in problems))

    # Build output records (one per sample)
    records: list[dict] = []
    for sample in sorted(sample_rows.keys()):
        control, resolved_vcf, variant_index = sample_metadata[sample]
        is_control = sample in referenced_controls

        fastq_1_list = []
        fastq_2_list = []
        for row_data in sample_rows[sample]:
            fastq_1_list.append(row_data["fastq_1"])
            if row_data["fastq_2"]:
                fastq_2_list.append(row_data["fastq_2"])

        records.append(
            {
                "sample": sample,
                "single_end": sample_layout[sample] == "SE",
                "fastq_1": fastq_1_list,
                "fastq_2": fastq_2_list,
                "control": control,
                "variant_vcf": resolved_vcf,
                "variant_index": variant_index,
                "is_control": is_control,
            }
        )

    return records


def write_samples_json(samples: list[dict], output_json: Path) -> None:
    """
    Write the records as an indented JSON list (indent=2, trailing newline).

    Args:
        samples: List of sample records from validate_samplesheet.
        output_json: Path to output JSON file.
    """
    with output_json.open("w") as handle:
        json.dump(samples, handle, indent=2)
        handle.write("\n")
