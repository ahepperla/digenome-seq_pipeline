#!/usr/bin/env python3
"""Plan deterministic, read-balanced genomic intervals for cleavage calling."""

from __future__ import annotations

import argparse
import bisect
import csv
import sys
from fractions import Fraction
from pathlib import Path

from genome_blacklist import load_genome_blacklist

try:
    import pysam
except ImportError:  # pragma: no cover - exercised by the command-line guard
    pysam = None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Split mapped BAM contigs into coordinate intervals balanced by "
            "estimated callable mapped-record counts."
        )
    )
    parser.add_argument("--bam", required=True)
    parser.add_argument("--chunks", type=int, required=True)
    parser.add_argument("--padding", type=int, required=True)
    parser.add_argument("--genome-blacklist")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--plan", required=True)
    return parser


def _rounded_fraction(value: Fraction) -> int:
    return (value.numerator + value.denominator // 2) // value.denominator


def _coordinate_for_callable_offset(
    callable_segments: list[tuple[int, int]],
    callable_offset: Fraction,
) -> int:
    if callable_offset <= 0:
        return callable_segments[0][0]
    remaining = callable_offset
    for start, end in callable_segments:
        segment_length = end - start
        if remaining <= segment_length:
            return start + _rounded_fraction(remaining)
        remaining -= segment_length
    return callable_segments[-1][1]


def plan_chunks(
    bam_path: str,
    requested_chunks: int,
    padding: int = 0,
    genome_blacklist_path: str | None = None,
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
    genome_blacklist = (
        load_genome_blacklist(genome_blacklist_path, lengths)
        if genome_blacklist_path
        else None
    )

    if not weighted:
        raise ValueError("BAM index reports no mapped alignments")

    for item in weighted:
        contig = str(item["contig"])
        contig_length = int(item["length"])
        callable_segments = (
            genome_blacklist.subtract(contig, 0, contig_length)
            if genome_blacklist is not None
            else [(0, contig_length)]
        )
        callable_bases = sum(
            end - start for start, end in callable_segments
        )
        item["callable_segments"] = callable_segments
        item["callable_bases"] = callable_bases
        item["callable_mapped_records"] = Fraction(
            int(item["mapped_records"]) * callable_bases,
            contig_length,
        )

    total_mapped = sum(int(item["mapped_records"]) for item in weighted)
    total_callable_bases = sum(
        int(item["callable_bases"]) for item in weighted
    )
    total_callable_weight = sum(
        (
            item["callable_mapped_records"]
            for item in weighted
        ),
        Fraction(0),
    )
    if total_callable_weight > 0:
        chunk_count = min(
            requested_chunks,
            total_mapped,
            total_callable_bases,
        )
        boundaries = [
            total_callable_weight * index / chunk_count
            for index in range(1, chunk_count)
        ]
    else:
        chunk_count = 1
        boundaries = []
    chunks: list[dict[str, object]] = [
        {
            "chunk_id": index,
            "mapped_records": Fraction(0),
            "callable_mapped_records": Fraction(0),
            "intervals": [],
        }
        for index in range(chunk_count)
    ]

    cumulative_callable_weight = Fraction(0)
    boundary_index = 0
    for item in weighted:
        contig = str(item["contig"])
        contig_length = int(item["length"])
        mapped_records = int(item["mapped_records"])
        callable_segments = list(item["callable_segments"])
        contig_start_weight = cumulative_callable_weight
        contig_end_weight = (
            contig_start_weight + item["callable_mapped_records"]
        )
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
            callable_offset = (
                boundaries[scan_boundary_index] - contig_start_weight
            ) * contig_length / mapped_records
            coordinate = _coordinate_for_callable_offset(
                callable_segments,
                callable_offset,
            )
            coordinate = max(
                cuts[-1] + 1,
                min(contig_length - 1, coordinate),
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
            excluded_bases = (
                genome_blacklist.overlap_bases(contig, start, end)
                if genome_blacklist is not None
                else 0
            )
            callable_bases = end - start - excluded_bases
            callable_weight = Fraction(
                mapped_records * callable_bases,
                contig_length,
            )
            callable_bases_before = start - (
                genome_blacklist.overlap_bases(contig, 0, start)
                if genome_blacklist is not None
                else 0
            )
            midpoint_weight = (
                contig_start_weight
                + Fraction(
                    mapped_records * callable_bases_before,
                    contig_length,
                )
                + callable_weight / 2
            )
            chunk_index = bisect.bisect_right(
                boundaries,
                midpoint_weight,
            )
            chunks[chunk_index]["mapped_records"] += interval_weight
            chunks[chunk_index][
                "callable_mapped_records"
            ] += callable_weight
            chunks[chunk_index]["intervals"].append(
                {
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "scan_start": max(0, start - padding),
                    "scan_end": min(contig_length, end + padding),
                    "callable_bases": callable_bases,
                    "excluded_bases": excluded_bases,
                }
            )
        cumulative_callable_weight = contig_end_weight

    nonempty_chunks = [
        chunk for chunk in chunks if chunk["intervals"]
    ]
    for chunk_id, chunk in enumerate(nonempty_chunks):
        chunk["chunk_id"] = chunk_id
        chunk["mapped_records"] = _rounded_fraction(
            chunk["mapped_records"]
        )
        chunk["callable_mapped_records"] = _rounded_fraction(
            chunk["callable_mapped_records"]
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
                "estimated_callable_mapped_records",
                "owned_bases",
                "callable_bases",
                "excluded_bases",
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
                        "callable_bases",
                        "excluded_bases",
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
                            interval["callable_bases"],
                            interval["excluded_bases"],
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
                    chunk["callable_mapped_records"],
                    sum(
                        int(interval["end"]) - int(interval["start"])
                        for interval in intervals
                    ),
                    sum(
                        int(interval["callable_bases"])
                        for interval in intervals
                    ),
                    sum(
                        int(interval["excluded_bases"])
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
        chunks = plan_chunks(
            args.bam,
            args.chunks,
            args.padding,
            args.genome_blacklist,
        )
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
