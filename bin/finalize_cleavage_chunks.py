#!/usr/bin/env python3
"""Merge cleavage-call chunks and apply sample-wide statistics and filters."""

from __future__ import annotations

import argparse
import gzip
import json
import sqlite3
import sys
import tempfile
from pathlib import Path
from typing import Any, Iterable

import call_cleavage as cleavage


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Merge interval-chunk cleavage calls, calculate global q-values, "
            "and write final pipeline outputs."
        )
    )
    parser.add_argument("--sample", required=True)
    parser.add_argument(
        "--analysis", choices=sorted(cleavage.ANALYSES), required=True
    )
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument(
        "--raw-fragment", action="append", required=True
    )
    parser.add_argument(
        "--chunk-summary", action="append", required=True
    )
    return parser


def load_summaries(
    paths: Iterable[str],
    sample: str,
    analysis: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    summaries = [json.loads(Path(path).read_text()) for path in paths]
    if not summaries:
        raise ValueError("No chunk summaries were supplied")

    chunk_ids: set[str] = set()
    reference_parameters = dict(summaries[0]["parameters"])
    expected_contigs = summaries[0].get("expected_contigs")
    if not expected_contigs:
        raise ValueError(
            "Chunk summaries do not declare the expected BAM contigs"
        )
    if len(expected_contigs) != len(set(expected_contigs)):
        raise ValueError("Expected BAM contig list contains duplicates")
    expected_contig_lengths = summaries[0].get(
        "expected_contig_lengths"
    )
    if (
        not isinstance(expected_contig_lengths, dict)
        or set(expected_contig_lengths) != set(expected_contigs)
    ):
        raise ValueError(
            "Chunk summaries do not declare all expected BAM contig lengths"
        )

    assigned_intervals: dict[
        str, list[tuple[int, int, str]]
    ] = {
        contig: [] for contig in expected_contigs
    }
    for summary in summaries:
        if summary.get("sample") != sample:
            raise ValueError(
                "Chunk summary sample mismatch: "
                f"{summary.get('sample')} != {sample}"
            )
        if summary.get("analysis") != analysis:
            raise ValueError(
                "Chunk summary analysis mismatch: "
                f"{summary.get('analysis')} != {analysis}"
            )
        chunk_id = summary.get("chunk_id")
        if not chunk_id or chunk_id in chunk_ids:
            raise ValueError(
                f"Missing or duplicate chunk identifier: {chunk_id}"
            )
        chunk_ids.add(chunk_id)
        if summary.get("parameters") != reference_parameters:
            raise ValueError(
                "Chunk summaries were generated with different parameters"
            )
        if summary.get("expected_contigs") != expected_contigs:
            raise ValueError(
                "Chunk summaries disagree about the expected BAM contigs"
            )
        if summary.get(
            "expected_contig_lengths"
        ) != expected_contig_lengths:
            raise ValueError(
                "Chunk summaries disagree about expected BAM contig lengths"
            )
        chunk_contigs = summary.get("contigs")
        if not isinstance(chunk_contigs, list) or not chunk_contigs:
            raise ValueError(
                f"Chunk {chunk_id} has no declared contigs"
            )
        chunk_intervals = summary.get("intervals")
        if not isinstance(chunk_intervals, list) or not chunk_intervals:
            raise ValueError(
                f"Chunk {chunk_id} has no declared intervals"
            )
        interval_contigs: list[str] = []
        for interval in chunk_intervals:
            if not isinstance(interval, dict):
                raise ValueError(
                    f"Chunk {chunk_id} contains an invalid interval"
                )
            contig = interval.get("contig")
            if contig not in expected_contig_lengths:
                raise ValueError(
                    f"Chunk {chunk_id} contains unexpected contig: {contig}"
                )
            try:
                start = int(interval["start"])
                end = int(interval["end"])
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError(
                    f"Chunk {chunk_id} contains invalid interval coordinates"
                ) from exc
            contig_length = int(expected_contig_lengths[contig])
            if not 0 <= start < end <= contig_length:
                raise ValueError(
                    f"Chunk {chunk_id} interval is outside {contig}: "
                    f"{start}-{end} of {contig_length}"
                )
            assigned_intervals[contig].append((start, end, chunk_id))
            interval_contigs.append(contig)
        if list(dict.fromkeys(interval_contigs)) != chunk_contigs:
            raise ValueError(
                f"Chunk {chunk_id} contig and interval declarations disagree"
            )

    for contig in expected_contigs:
        intervals = sorted(assigned_intervals[contig])
        if not intervals:
            raise ValueError(
                f"Chunk interval coverage is missing: {contig}:0-"
                f"{expected_contig_lengths[contig]}"
            )
        cursor = 0
        previous_chunk = ""
        for start, end, chunk_id in intervals:
            if start > cursor:
                raise ValueError(
                    f"Chunk interval coverage has a gap on {contig}: "
                    f"{cursor}-{start}"
                )
            if start < cursor:
                raise ValueError(
                    f"Chunk interval ownership overlaps on {contig} at "
                    f"{start}-{min(cursor, end)} between "
                    f"{previous_chunk} and {chunk_id}"
                )
            cursor = end
            previous_chunk = chunk_id
        contig_length = int(expected_contig_lengths[contig])
        if cursor < contig_length:
            raise ValueError(
                f"Chunk interval coverage has a gap on {contig}: "
                f"{cursor}-{contig_length}"
            )
    return summaries, reference_parameters


def create_database(path: str) -> sqlite3.Connection:
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    connection.execute(
        """
        CREATE TABLE calls (
            row_id INTEGER PRIMARY KEY,
            contig TEXT NOT NULL,
            position_0based INTEGER NOT NULL,
            strand TEXT NOT NULL,
            fisher_p REAL,
            fisher_q REAL,
            payload TEXT NOT NULL,
            UNIQUE (contig, position_0based, strand)
        )
        """
    )
    return connection


def load_fragments(
    connection: sqlite3.Connection,
    paths: Iterable[str],
    summaries: Iterable[dict[str, Any]],
) -> int:
    owned_intervals = {
        summary["chunk_id"]: [
            (
                interval["contig"],
                int(interval["start"]),
                int(interval["end"]),
            )
            for interval in summary["intervals"]
        ]
        for summary in summaries
    }
    inserted = 0
    batch: list[tuple[str, int, str, float | None, str]] = []
    for path in paths:
        with gzip.open(path, "rt") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip():
                    continue
                try:
                    row = json.loads(line)
                except json.JSONDecodeError as exc:
                    raise ValueError(
                        f"Invalid JSON in {path} line {line_number}: {exc}"
                    ) from exc
                chunk_id = row.get("_chunk_id")
                if chunk_id not in owned_intervals:
                    raise ValueError(
                        f"Raw row in {path} line {line_number} has unknown "
                        f"chunk identifier: {chunk_id}"
                    )
                contig = row["contig"]
                position = int(row["position_0based"])
                if not any(
                    owned_contig == contig
                    and start <= position < end
                    for owned_contig, start, end
                    in owned_intervals[chunk_id]
                ):
                    raise ValueError(
                        f"Raw row in {path} line {line_number} is outside "
                        f"chunk {chunk_id} ownership: {contig}:{position}"
                    )
                batch.append(
                    (
                        contig,
                        position,
                        row["strand"],
                        row.get("control_fisher_p"),
                        json.dumps(row, sort_keys=True),
                    )
                )
                if len(batch) >= 1000:
                    connection.executemany(
                        """
                        INSERT INTO calls (
                            contig, position_0based, strand, fisher_p, payload
                        ) VALUES (?, ?, ?, ?, ?)
                        """,
                        batch,
                    )
                    inserted += len(batch)
                    batch.clear()
    if batch:
        connection.executemany(
            """
            INSERT INTO calls (
                contig, position_0based, strand, fisher_p, payload
            ) VALUES (?, ?, ?, ?, ?)
            """,
            batch,
        )
        inserted += len(batch)
    connection.commit()
    return inserted


def assign_global_q_values(connection: sqlite3.Connection) -> None:
    controlled = connection.execute(
        "SELECT COUNT(*) FROM calls WHERE fisher_p IS NOT NULL"
    ).fetchone()[0]
    if controlled == 0:
        return

    running = 1.0
    updates: list[tuple[float, int]] = []
    query = connection.execute(
        """
        SELECT row_id, fisher_p
        FROM calls
        WHERE fisher_p IS NOT NULL
        ORDER BY fisher_p DESC, row_id DESC
        """
    )
    for offset, (row_id, fisher_p) in enumerate(query):
        rank = controlled - offset
        running = min(running, float(fisher_p) * controlled / rank)
        updates.append((min(1.0, running), row_id))
        if len(updates) >= 1000:
            connection.executemany(
                "UPDATE calls SET fisher_q = ? WHERE row_id = ?",
                updates,
            )
            updates.clear()
    if updates:
        connection.executemany(
            "UPDATE calls SET fisher_q = ? WHERE row_id = ?",
            updates,
        )
    connection.commit()


def finalized_rows(
    connection: sqlite3.Connection,
    args: argparse.Namespace,
) -> Iterable[dict[str, Any]]:
    query = connection.execute(
        """
        SELECT payload, fisher_q
        FROM calls
        ORDER BY contig, position_0based, strand
        """
    )
    for payload, fisher_q in query:
        row = json.loads(payload)
        row["control_fisher_q"] = fisher_q
        cleavage.apply_filters_to_row(row, args)
        yield row


def finalize_chunks(args: argparse.Namespace) -> dict[str, Any]:
    if len(args.raw_fragment) != len(args.chunk_summary):
        raise ValueError(
            "Chunk input-count mismatch: "
            f"{len(args.raw_fragment)} raw fragment(s), "
            f"{len(args.chunk_summary)} summary file(s)"
        )

    summaries, parameters = load_summaries(
        sorted(args.chunk_summary),
        args.sample,
        args.analysis,
    )
    expected_rows = sum(
        int(summary["reported_rows"]) for summary in summaries
    )
    candidate_count = sum(
        int(summary["candidate_count"]) for summary in summaries
    )

    parameters.update(
        {
            "sample": args.sample,
            "analysis": args.analysis,
            "output_prefix": args.output_prefix,
            "cleavage_chunks": len(summaries),
        }
    )
    caller_args = argparse.Namespace(**parameters)
    cleavage.validate_args(caller_args)

    with tempfile.NamedTemporaryFile(
        prefix=f"{args.sample}.cleavage.",
        suffix=".sqlite",
        dir=".",
        delete=False,
    ) as temporary:
        database_path = temporary.name

    connection = create_database(database_path)
    try:
        observed_rows = load_fragments(
            connection,
            sorted(args.raw_fragment),
            summaries,
        )
        if observed_rows != expected_rows:
            raise ValueError(
                "Chunk row-count mismatch: "
                f"expected {expected_rows}, observed {observed_rows}"
            )
        assign_global_q_values(connection)
        return cleavage.write_output_stream(
            finalized_rows(connection, caller_args),
            candidate_count,
            caller_args,
        )
    finally:
        connection.close()
        Path(database_path).unlink(missing_ok=True)


def main() -> None:
    args = build_parser().parse_args()
    try:
        qc = finalize_chunks(args)
    except (OSError, ValueError, RuntimeError, sqlite3.Error) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    print(
        f"{args.analysis}: finalized {qc['reported_candidates']} candidate(s), "
        f"{qc['high_confidence_calls']} high-confidence call(s)."
    )


if __name__ == "__main__":
    main()
