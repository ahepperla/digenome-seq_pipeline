#!/usr/bin/env python3
"""Write the tiny smoke-test inputs to tests/fixtures/generated/.

The paired FASTQs align uniquely to tiny.fa and contain:
- a DSB at 250: 11 forward reads start there and 11 reverse reads end there;
- an SSB at 300: 11 forward reads start there and their mates vary.
Every fragment has distinct coordinates, so duplicate marking keeps them all.

Usage: python3 tests/fixtures/build_smoke_fixture.py
"""

from __future__ import annotations

import gzip
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "generated"
READ_LENGTH = 60
DSB_POSITION = 250
SSB_POSITION = 300
PILEUP_COUNT = 11


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def reference() -> str:
    lines = (HERE / "tiny.fa").read_text().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


def fragments() -> list[tuple[str, int, int]]:
    """(name, forward read start, reverse read end + 1) for every fragment."""
    dsb_forward = [(f"dsb_forward_{i + 1}", DSB_POSITION, 380 + i * 4) for i in range(PILEUP_COUNT)]
    dsb_reverse = [(f"dsb_reverse_{i + 1}", 40 + i * 7, DSB_POSITION + 1) for i in range(PILEUP_COUNT)]
    ssb = [(f"ssb_forward_{i + 1}", SSB_POSITION, 440 + i * 5) for i in range(PILEUP_COUNT)]
    return dsb_forward + dsb_reverse + ssb


def read_pair(sequence: str, forward_start: int, reverse_end: int) -> tuple[str, str]:
    read1 = sequence[forward_start:forward_start + READ_LENGTH]
    read2 = reverse_complement(sequence[reverse_end - READ_LENGTH:reverse_end])
    if len(read1) != READ_LENGTH or len(read2) != READ_LENGTH:
        raise ValueError("Read coordinates are outside tiny.fa")
    return read1, read2


def build(directory: Path = OUTPUT) -> Path:
    """Write tiny_R1/R2.fastq.gz and samplesheet.csv; return the samplesheet."""
    directory.mkdir(parents=True, exist_ok=True)
    sequence = reference()
    read1_path, read2_path = directory / "tiny_R1.fastq.gz", directory / "tiny_R2.fastq.gz"
    with gzip.open(read1_path, "wt") as read1_file, gzip.open(read2_path, "wt") as read2_file:
        for name, forward_start, reverse_end in fragments():
            read1, read2 = read_pair(sequence, forward_start, reverse_end)
            read1_file.write(f"@{name}/1\n{read1}\n+\n{'I' * READ_LENGTH}\n")
            read2_file.write(f"@{name}/2\n{read2}\n+\n{'I' * READ_LENGTH}\n")
    samplesheet = directory / "samplesheet.csv"
    samplesheet.write_text(f"sample,fastq_1,fastq_2\nTiny,{read1_path},{read2_path}\n")
    return samplesheet


if __name__ == "__main__":
    print(f"Wrote {build()} ({len(fragments())} read pairs)")
