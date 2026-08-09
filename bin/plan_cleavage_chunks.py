#!/usr/bin/env python3
"""Balance BAM contigs into deterministic cleavage-calling chunks."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

try:
    import pysam
except ImportError:  # pragma: no cover - exercised by the command-line guard
    pysam = None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Group contigs into read-balanced chunks for parallel cleavage "
            "calling."
        )
    )
    parser.add_argument("--bam", required=True)
    parser.add_argument("--chunks", type=int, required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--plan", required=True)
    return parser


def plan_chunks(
    bam_path: str,
    requested_chunks: int,
) -> list[dict[str, object]]:
    if pysam is None:
        raise RuntimeError(
            "pysam is required. Run this command in the cleavage container."
        )
    if requested_chunks < 1:
        raise ValueError("--chunks must be at least 1")

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        if not bam.has_index():
            raise ValueError(f"BAM is not indexed: {bam_path}")
        order = {
            contig: index for index, contig in enumerate(bam.references)
        }
        weighted = [
            (stat.contig, int(stat.mapped), order[stat.contig])
            for stat in bam.get_index_statistics()
            if stat.mapped > 0
        ]

    if not weighted:
        raise ValueError("BAM index reports no mapped alignments")

    chunk_count = min(requested_chunks, len(weighted))
    chunks = [
        {"chunk_id": index, "mapped_records": 0, "contigs": []}
        for index in range(chunk_count)
    ]
    for contig, mapped, contig_order in sorted(
        weighted, key=lambda item: (-item[1], item[2])
    ):
        target = min(
            chunks,
            key=lambda chunk: (
                int(chunk["mapped_records"]),
                int(chunk["chunk_id"]),
            ),
        )
        target["mapped_records"] = (
            int(target["mapped_records"]) + mapped
        )
        target["contigs"].append((contig_order, contig))

    for chunk in chunks:
        chunk["contigs"] = [
            contig
            for _order, contig in sorted(chunk["contigs"])
        ]
    return chunks


def write_plan(
    chunks: list[dict[str, object]],
    output_dir: Path,
    plan_path: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    width = max(3, len(str(len(chunks) - 1)))
    with plan_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["chunk_id", "contig_count", "mapped_records", "contigs"]
        )
        for chunk in chunks:
            chunk_id = f"chunk_{int(chunk['chunk_id']):0{width}d}"
            contigs = list(chunk["contigs"])
            contig_path = output_dir / f"{chunk_id}.contigs.txt"
            contig_path.write_text("".join(f"{contig}\n" for contig in contigs))
            writer.writerow(
                [
                    chunk_id,
                    len(contigs),
                    chunk["mapped_records"],
                    ",".join(contigs),
                ]
            )


def main() -> None:
    args = build_parser().parse_args()
    try:
        chunks = plan_chunks(args.bam, args.chunks)
        write_plan(
            chunks,
            Path(args.output_dir),
            Path(args.plan),
        )
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    print(
        f"Planned {len(chunks)} cleavage chunk(s) from {args.bam}."
    )


if __name__ == "__main__":
    main()
