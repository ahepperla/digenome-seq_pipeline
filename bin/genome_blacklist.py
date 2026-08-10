#!/usr/bin/env python3
"""Parse and query optional genome blacklist BED intervals."""

from __future__ import annotations

import gzip
import hashlib
from bisect import bisect_right
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO


@dataclass(frozen=True)
class GenomeBlacklist:
    intervals: dict[str, tuple[tuple[int, int], ...]]
    starts: dict[str, tuple[int, ...]]
    sha256: str
    interval_count: int
    excluded_bases: int

    def contains(self, contig: str, position: int) -> bool:
        contig_intervals = self.intervals.get(contig, ())
        if not contig_intervals:
            return False
        index = bisect_right(self.starts[contig], position) - 1
        return (
            index >= 0
            and position < contig_intervals[index][1]
        )

    def subtract(
        self,
        contig: str,
        start: int,
        end: int,
    ) -> list[tuple[int, int]]:
        if start >= end:
            return []
        callable_intervals: list[tuple[int, int]] = []
        cursor = start
        for blocked_start, blocked_end in self.intervals.get(contig, ()):
            if blocked_end <= cursor:
                continue
            if blocked_start >= end:
                break
            if blocked_start > cursor:
                callable_intervals.append(
                    (cursor, min(blocked_start, end))
                )
            cursor = max(cursor, blocked_end)
            if cursor >= end:
                break
        if cursor < end:
            callable_intervals.append((cursor, end))
        return callable_intervals

    def overlap_bases(self, contig: str, start: int, end: int) -> int:
        callable_bases = sum(
            interval_end - interval_start
            for interval_start, interval_end
            in self.subtract(contig, start, end)
        )
        return max(0, end - start - callable_bases)


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def _merge_intervals(
    intervals: Iterable[tuple[int, int]],
) -> tuple[tuple[int, int], ...]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return tuple(merged)


def load_genome_blacklist(
    path: str,
    reference_lengths: dict[str, int],
) -> GenomeBlacklist:
    blacklist_path = Path(path)
    if not blacklist_path.is_file():
        raise ValueError(f"Genome blacklist does not exist: {path}")
    raw_intervals: dict[str, list[tuple[int, int]]] = {}
    with _open_text(blacklist_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if (
                not stripped
                or stripped.startswith("#")
                or stripped.startswith("track ")
                or stripped.startswith("browser ")
            ):
                continue
            fields = stripped.split()
            if len(fields) < 3:
                raise ValueError(
                    f"Genome blacklist line {line_number} has fewer than "
                    "three BED columns"
                )
            contig = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"Genome blacklist line {line_number} has invalid "
                    "coordinates"
                ) from exc
            if contig not in reference_lengths:
                raise ValueError(
                    f"Genome blacklist contig is absent from the BAM: "
                    f"{contig}"
                )
            contig_length = int(reference_lengths[contig])
            if not 0 <= start < end <= contig_length:
                raise ValueError(
                    f"Genome blacklist interval is outside {contig}: "
                    f"{start}-{end} of {contig_length}"
                )
            raw_intervals.setdefault(contig, []).append((start, end))

    merged = {
        contig: _merge_intervals(intervals)
        for contig, intervals in raw_intervals.items()
    }
    interval_count = sum(len(intervals) for intervals in merged.values())
    if interval_count == 0:
        raise ValueError(
            f"Genome blacklist contains no BED intervals: {path}"
        )
    excluded_bases = sum(
        end - start
        for intervals in merged.values()
        for start, end in intervals
    )
    digest = hashlib.sha256(blacklist_path.read_bytes()).hexdigest()
    return GenomeBlacklist(
        intervals=merged,
        starts={
            contig: tuple(start for start, _end in intervals)
            for contig, intervals in merged.items()
        },
        sha256=digest,
        interval_count=interval_count,
        excluded_bases=excluded_bases,
    )
