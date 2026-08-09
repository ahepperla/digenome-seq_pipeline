#!/usr/bin/env python3
"""Create deterministic DSB/SSB FASTQs and a smoke-test samplesheet."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

HERE = Path(__file__).resolve().parent
READ_LENGTH = 60
DSB_POSITION = 250
SSB_POSITION = 300
PILEUP_COUNT = 11


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def read_reference() -> str:
    return "".join(
        line.strip()
        for line in (HERE / "tiny.fa").read_text().splitlines()
        if not line.startswith(">")
    )


def paired_reads(
    reference: str,
    name: str,
    forward_start: int,
    reverse_end: int,
) -> tuple[tuple[str, str], tuple[str, str]]:
    read1 = reference[forward_start : forward_start + READ_LENGTH]
    read2_reference = reference[reverse_end - READ_LENGTH : reverse_end]
    if len(read1) != READ_LENGTH or len(read2_reference) != READ_LENGTH:
        raise ValueError(f"Read coordinates are outside tiny.fa for {name}")
    return (
        (f"@{name}/1", read1),
        (f"@{name}/2", reverse_complement(read2_reference)),
    )


def write_fastq(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w") as handle:
        for name, sequence in records:
            handle.write(
                f"{name}\n{sequence}\n+\n{'I' * len(sequence)}\n"
            )


def gzip_file(source: Path, destination: Path) -> None:
    with source.open("rb") as input_handle, gzip.open(destination, "wb") as output:
        output.write(input_handle.read())


def main() -> None:
    reference = read_reference()
    read1_records: list[tuple[str, str]] = []
    read2_records: list[tuple[str, str]] = []

    # Forward reads pile up at the DSB while mate endpoints vary.
    for index in range(PILEUP_COUNT):
        read1, read2 = paired_reads(
            reference,
            f"dsb_forward_{index + 1}",
            DSB_POSITION,
            380 + index * 4,
        )
        read1_records.append(read1)
        read2_records.append(read2)

    # Reverse reads pile up at the same DSB while mate starts vary.
    for index in range(PILEUP_COUNT):
        read1, read2 = paired_reads(
            reference,
            f"dsb_reverse_{index + 1}",
            40 + index * 7,
            DSB_POSITION + 1,
        )
        read1_records.append(read1)
        read2_records.append(read2)

    # A separate one-strand pileup provides a high-confidence SSB in nDigenome.
    for index in range(PILEUP_COUNT):
        read1, read2 = paired_reads(
            reference,
            f"ssb_forward_{index + 1}",
            SSB_POSITION,
            440 + index * 5,
        )
        read1_records.append(read1)
        read2_records.append(read2)

    plain_r1 = HERE / "tiny_R1.fastq"
    plain_r2 = HERE / "tiny_R2.fastq"
    write_fastq(plain_r1, read1_records)
    write_fastq(plain_r2, read2_records)

    r1 = HERE / "tiny_R1.fastq.gz"
    r2 = HERE / "tiny_R2.fastq.gz"
    gzip_file(plain_r1, r1)
    gzip_file(plain_r2, r2)
    with (HERE / "tiny_samplesheet.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample", "fastq_1", "fastq_2"])
        writer.writerow(["Tiny", r1.resolve(), r2.resolve()])
    print(
        f"{HERE / 'tiny_samplesheet.csv'} "
        f"({len(read1_records)} pairs; DSB={DSB_POSITION}, SSB={SSB_POSITION})"
    )


if __name__ == "__main__":
    main()
