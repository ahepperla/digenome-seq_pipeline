"""Check the input samplesheet and group its rows into one record per sample.

Columns: sample, fastq_1, fastq_2 (required); control, variant_vcf (optional).
Each row is one FASTQ pair, or one single-end FASTQ when fastq_2 is blank.
Rows that share a sample name are lanes of that sample. Every problem in the
sheet is reported at once.
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

REQUIRED_COLUMNS = ["sample", "fastq_1", "fastq_2"]
OPTIONAL_COLUMNS = ["control", "variant_vcf"]
NAME_PATTERN = re.compile(r"^[A-Za-z0-9_.-]+$")
FASTQ_SUFFIX = re.compile(r"\.(fastq|fq)\.gz$", re.IGNORECASE)
VCF_SUFFIX = re.compile(r"\.vcf\.gz$", re.IGNORECASE)


def validate_samplesheet(input_csv: Path, analysis: str) -> list[dict]:
    """Return one record per sample, or raise ValueError listing every problem."""
    if analysis not in ("digenome", "ndigenome"):
        raise ValueError(f"Unknown analysis '{analysis}'. Expected digenome or ndigenome")
    problems: list[str] = []
    samples: dict[str, dict] = {}
    fastq_owners: dict[str, str] = {}
    for line_number, row in read_rows(Path(input_csv)):
        check_row(line_number, row, analysis, samples, fastq_owners, problems)
    if not samples:
        problems.append("Samplesheet contains no usable data rows")
    problems += control_problems(samples)
    if problems:
        raise ValueError(
            "Samplesheet validation failed:\n" + "\n".join(f"  - {problem}" for problem in problems)
        )
    controls = {record["control"] for record in samples.values()}
    for name, record in samples.items():
        record["is_control"] = name in controls
    return [samples[name] for name in sorted(samples)]


def write_samples_json(samples: list[dict], output_json: Path) -> None:
    Path(output_json).write_text(json.dumps(samples, indent=2) + "\n")


def read_rows(input_csv: Path) -> list[tuple[int, dict[str, str]]]:
    """Read the sheet after checking its header. Values are stripped; the
    header is line 1."""
    if not input_csv.is_file():
        raise ValueError(f"Samplesheet does not exist: {input_csv}")
    with input_csv.open(newline="") as handle:
        reader = csv.reader(handle)
        header = [name.strip() for name in next(reader, [])]
        check_header(header)
        return [
            (line_number, {name: value.strip() for name, value in zip(header, values)})
            for line_number, values in enumerate(reader, start=2)
            if any(value.strip() for value in values)
        ]


def check_header(header: list[str]) -> None:
    if not header:
        raise ValueError("Samplesheet is empty or missing a header row")
    if len(set(header)) != len(header):
        raise ValueError("Samplesheet contains duplicate column names")
    unknown = sorted(set(header) - set(REQUIRED_COLUMNS + OPTIONAL_COLUMNS))
    if unknown:
        raise ValueError(
            f"Samplesheet contains unknown column(s): {', '.join(unknown)}. "
            f"Allowed columns are: {', '.join(REQUIRED_COLUMNS + OPTIONAL_COLUMNS)}"
        )
    missing = [name for name in REQUIRED_COLUMNS if name not in header]
    if missing:
        raise ValueError(f"Samplesheet is missing required column(s): {', '.join(missing)}")


def check_row(
    line_number: int,
    row: dict[str, str],
    analysis: str,
    samples: dict[str, dict],
    fastq_owners: dict[str, str],
    problems: list[str],
) -> None:
    """Check one row and add it to `samples` (keyed by sample name)."""
    sample, control = row.get("sample", ""), row.get("control", "")
    fastq_1, fastq_2 = row.get("fastq_1", ""), row.get("fastq_2", "")

    def problem(message: str) -> None:
        problems.append(f"line {line_number}: {message}")

    if not sample:
        problem("sample is blank")
        return
    for label, name in (("sample", sample), ("control", control)):
        if name and not NAME_PATTERN.fullmatch(name):
            problem(
                f"{label} '{name}' contains unsupported characters. "
                "Use only letters, numbers, dots, underscores, and hyphens."
            )
    if not fastq_1:
        problem("fastq_1 is blank")
        return
    single_end = not fastq_2
    if analysis == "ndigenome" and single_end:
        problem("nDigenome requires paired-end data")

    fastqs = {}
    for label, path in (("fastq_1", fastq_1), ("fastq_2", fastq_2)):
        if path:
            fastqs[label] = check_fastq(label, path, f"line {line_number} ({sample})", fastq_owners, problem)
    vcf, index = check_vcf(row.get("variant_vcf", ""), problem)

    record = samples.setdefault(sample, {
        "sample": sample,
        "single_end": single_end,
        "fastq_1": [],
        "fastq_2": [],
        "control": control,
        "variant_vcf": vcf,
        "variant_index": index,
    })
    if record["single_end"] != single_end:
        problem(f"sample '{sample}' mixes single-end and paired-end rows")
    if (record["control"], record["variant_vcf"], record["variant_index"]) != (control, vcf, index):
        problem(f"sample '{sample}' has inconsistent control or variant_vcf values across lanes")
    record["fastq_1"].append(fastqs["fastq_1"])
    if "fastq_2" in fastqs:
        record["fastq_2"].append(fastqs["fastq_2"])


def check_fastq(label: str, path: str, owner: str, fastq_owners: dict[str, str], problem) -> str:
    """Check one FASTQ path and return it resolved. A file may be used once
    in the whole sheet, which also rejects fastq_1 == fastq_2."""
    if not FASTQ_SUFFIX.search(path):
        problem(f"{label} must end with .fastq.gz or .fq.gz: {path}")
    if not Path(path).expanduser().is_file():
        problem(f"{label} file does not exist: {path}")
    resolved = str(Path(path).expanduser().resolve())
    if resolved in fastq_owners:
        problem(f"{label} reuses FASTQ '{resolved}' already used on {fastq_owners[resolved]}")
    else:
        fastq_owners[resolved] = f"{owner} as {label}"
    return resolved


def check_vcf(path: str, problem) -> tuple[str, str]:
    """Return the resolved VCF and its .tbi or .csi index, or ("", "")."""
    if not path:
        return "", ""
    if not VCF_SUFFIX.search(path):
        problem(f"variant_vcf must be a bgzip-compressed .vcf.gz file: {path}")
    if not Path(path).expanduser().is_file():
        problem(f"variant_vcf file does not exist: {path}")
        return "", ""
    vcf = str(Path(path).expanduser().resolve())
    for suffix in (".tbi", ".csi"):
        if Path(vcf + suffix).is_file():
            return vcf, vcf + suffix
    problem(f"variant_vcf is not indexed with .tbi or .csi: {vcf}")
    return vcf, ""


def control_problems(samples: dict[str, dict]) -> list[str]:
    """A control must be another sample in the sheet, and must not have a
    control of its own."""
    problems = []
    for name, record in samples.items():
        control = record["control"]
        if not control:
            continue
        if control == name:
            problems.append(f"sample '{name}' cannot use itself as its control")
        elif control not in samples:
            problems.append(
                f"sample '{name}' references control '{control}', but no sample "
                f"named '{control}' exists in the samplesheet"
            )
        elif samples[control]["control"]:
            problems.append(
                f"sample '{control}' is the control of '{name}' but declares control "
                f"'{samples[control]['control']}'. Control rows must leave the control column blank."
            )
    return problems
