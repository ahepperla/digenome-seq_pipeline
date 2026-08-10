#!/usr/bin/env python3
"""Plan deterministic, read-balanced genomic intervals for cleavage calling."""

from __future__ import annotations

import argparse
import csv
import sys
from fractions import Fraction
from pathlib import Path

try:
    import pysam
except ImportError:  # pragma: no cover - exercised by the command-line guard
    pysam = None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Split mapped BAM contigs into coordinate intervals balanced by "
            "mapped-record counts."
        )
    )
    parser.add_argument("--bam", required=True)
    parser.add_argument("--chunks", type=int, required=True)
    parser.add_argument("--padding", type=int, required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--plan", required=True)
    return parser


def _rounded_fraction(value: Fraction) -> int:
    return (value.numerator + value.denominator // 2) // value.denominator


def plan_chunks(
    bam_path: str,
    requested_chunks: int,
    padding: int = 0,
) -> list[dict[str, object]]:
    if pysam is None:
        raise RuntimeError(
            "pysam is required. Run this command in the cleavage container."
        )
    if requested_chunks < 1:
        raise ValueError("--chunks must be at least 1")
    if padding < 0:
        raise ValueError("--padding must be non-negative")

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        if not bam.has_index():
            raise ValueError(f"BAM is not indexed: {bam_path}")
        lengths = dict(zip(bam.references, bam.lengths))
        weighted = [
            {
                "contig": stat.contig,
                "length": int(lengths[stat.contig]),
                "mapped_records": int(stat.mapped),
            }
            for stat in bam.get_index_statistics()
            if stat.mapped > 0
        ]

    if not weighted:
        raise ValueError("BAM index reports no mapped alignments")

    total_mapped = sum(int(item["mapped_records"]) for item in weighted)
    total_bases = sum(int(item["length"]) for item in weighted)
    chunk_count = min(requested_chunks, total_mapped, total_bases)
    boundaries = [
        Fraction(total_mapped * index, chunk_count)
        for index in range(1, chunk_count)
    ]
    chunks: list[dict[str, object]] = [
        {
            "chunk_id": index,
            "mapped_records": Fraction(0),
            "intervals": [],
        }
        for index in range(chunk_count)
    ]

    cumulative_mapped = 0
    boundary_index = 0
    for item in weighted:
        contig = str(item["contig"])
        contig_length = int(item["length"])
        mapped_records = int(item["mapped_records"])
        contig_start_weight = Fraction(cumulative_mapped)
        contig_end_weight = Fraction(cumulative_mapped + mapped_records)
        cuts = [0]

        while (
            boundary_index < len(boundaries)
            and boundaries[boundary_index] <= contig_start_weight
        ):
            boundary_index += 1
        scan_boundary_index = boundary_index
        while (
            scan_boundary_index < len(boundaries)
            and boundaries[scan_boundary_index] < contig_end_weight
        ):
            offset = (
                boundaries[scan_boundary_index] - contig_start_weight
            ) * contig_length / mapped_records
            coordinate = max(
                cuts[-1] + 1,
                min(contig_length - 1, _rounded_fraction(offset)),
            )
            if coordinate < contig_length:
                cuts.append(coordinate)
            scan_boundary_index += 1
        cuts.append(contig_length)
        boundary_index = scan_boundary_index

        for start, end in zip(cuts, cuts[1:]):
            interval_weight = Fraction(
                mapped_records * (end - start),
                contig_length,
            )
            midpoint_weight = contig_start_weight + Fraction(
                mapped_records * (start + end),
                2 * contig_length,
            )
            chunk_index = min(
                chunk_count - 1,
                int(midpoint_weight * chunk_count // total_mapped),
            )
            chunks[chunk_index]["mapped_records"] += interval_weight
            chunks[chunk_index]["intervals"].append(
                {
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "scan_start": max(0, start - padding),
                    "scan_end": min(contig_length, end + padding),
                }
            )
        cumulative_mapped += mapped_records

    nonempty_chunks = [
        chunk for chunk in chunks if chunk["intervals"]
    ]
    for chunk_id, chunk in enumerate(nonempty_chunks):
        chunk["chunk_id"] = chunk_id
        chunk["mapped_records"] = _rounded_fraction(
            chunk["mapped_records"]
        )
    return nonempty_chunks


def write_plan(
    chunks: list[dict[str, object]],
    output_dir: Path,
    plan_path: Path,
) -> None:
    if not chunks:
        raise ValueError("No cleavage chunks were planned")
    output_dir.mkdir(parents=True, exist_ok=True)
    width = max(3, len(str(len(chunks) - 1)))
    with plan_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "chunk_id",
                "interval_count",
                "contig_count",
                "estimated_mapped_records",
                "owned_bases",
                "intervals",
            ]
        )
        for chunk in chunks:
            chunk_id = f"chunk_{int(chunk['chunk_id']):0{width}d}"
            intervals = list(chunk["intervals"])
            interval_path = output_dir / f"{chunk_id}.intervals.tsv"
            with interval_path.open("w", newline="") as interval_handle:
                interval_writer = csv.writer(
                    interval_handle,
                    delimiter="\t",
                    lineterminator="\n",
                )
                interval_writer.writerow(
                    [
                        "contig",
                        "owner_start",
                        "owner_end",
                        "scan_start",
                        "scan_end",
                    ]
                )
                for interval in intervals:
                    interval_writer.writerow(
                        [
                            interval["contig"],
                            interval["start"],
                            interval["end"],
                            interval["scan_start"],
                            interval["scan_end"],
                        ]
                    )

            contigs = list(
                dict.fromkeys(str(interval["contig"]) for interval in intervals)
            )
            writer.writerow(
                [
                    chunk_id,
                    len(intervals),
                    len(contigs),
                    chunk["mapped_records"],
                    sum(
                        int(interval["end"]) - int(interval["start"])
                        for interval in intervals
                    ),
                    ",".join(
                        f"{interval['contig']}:"
                        f"{interval['start']}-{interval['end']}"
                        for interval in intervals
                    ),
                ]
            )


def main() -> None:
    args = build_parser().parse_args()
    try:
        chunks = plan_chunks(args.bam, args.chunks, args.padding)
        write_plan(
            chunks,
            Path(args.output_dir),
            Path(args.plan),
        )
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    print(
        f"Planned {len(chunks)} cleavage interval chunk(s) from {args.bam}."
    )


if __name__ == "__main__":
    main()
