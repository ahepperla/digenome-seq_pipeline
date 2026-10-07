"""Merge one sample's chunks and write its final outputs.

1. Check that all `chunks` summaries are present, agree on their inputs, and
   that the plan owns every mapped contig exactly once.
2. Stream-merge the sorted chunk files. nDigenome rows pass through;
   Digenome endpoint records are paired one contig at a time.
3. Compute Benjamini-Hochberg q-values over every controlled row in the
   sample, then apply the filters and write the files.
"""

from __future__ import annotations

import gzip
import heapq
import itertools
import json
import math
import tempfile
from pathlib import Path
from typing import Iterator

from . import digenome
from .call import record_key
from .output import apply_filters, write_outputs
from .regions import OwnedInterval, check_plan, plan_from_json, write_plan
from .settings import CallerSettings
from .stats import benjamini_hochberg

MULTIMAPPER_WARNING = (
    "Multimapper mode counts each read only at its selected primary placement; "
    "support can be diluted across equivalent reference copies."
)


def finalize(
    settings: CallerSettings,
    summary_paths: list[str],
    sample: str,
    control_sample: str,
    out_prefix: str,
) -> dict:
    chunks = load_summaries(summary_paths, settings, sample)
    summaries = [summary for _path, summary in chunks]
    first = summaries[0]
    contig_lengths = {name: length for name, length in first["contigs"]}
    contig_order = {name: index for index, name in enumerate(contig_lengths)}
    plan = plan_from_json(first["plan"])
    check_plan(plan, contig_lengths)
    write_plan(plan, f"{out_prefix}.cleavage_chunks.tsv")

    merged = heapq.merge(
        *(read_chunk(summary, path, plan[summary["chunk"]]) for path, summary in chunks),
        key=lambda record: record_key(record, contig_order),
    )
    with tempfile.TemporaryFile("w+") as staged:
        if settings.analysis == "ndigenome":
            p_values = stage_rows(unique(merged, contig_order), staged)
            candidates = sum(summary["candidates"] for summary in summaries)
        else:
            p_values, candidates = stage_pairs(unique(merged, contig_order), settings, sample, control_sample, staged)
        q_values = iter(benjamini_hochberg(p_values))
        staged.seek(0)
        report = {
            "schema_version": 1,
            "sample": sample,
            "analysis": settings.analysis,
            "control_sample": control_sample,
            "chunks": len(summaries),
            "settings": settings.to_dict(),
            "genome_blacklist": first["genome_blacklist"],
            "candidates_before_filters": candidates,
            "warnings": [MULTIMAPPER_WARNING] if settings.keep_multimappers else [],
        }
        return write_outputs(filtered_rows(staged, q_values, settings), out_prefix, settings, report)


def load_summaries(paths: list[str], settings: CallerSettings, sample: str) -> list[tuple[str, dict]]:
    """Read the chunk summaries as (path, summary) pairs ordered by chunk index,
    and check they belong together: one per index, same sample, settings,
    contigs, plan, and blacklist. Paths may arrive in any order."""
    chunks = sorted(
        ((path, json.loads(Path(path).read_text())) for path in paths),
        key=lambda item: item[1]["chunk"],
    )
    summaries = [summary for _path, summary in chunks]
    count = len(summaries)
    if [summary["chunk"] for summary in summaries] != list(range(count)):
        raise ValueError(f"Expected chunks 0..{count - 1}, found {[s['chunk'] for s in summaries]}")
    first = summaries[0]
    for summary in summaries:
        if summary["sample"] != sample:
            raise ValueError(f"Chunk {summary['chunk']} belongs to sample {summary['sample']}, not {sample}")
        if summary["chunks"] != count:
            raise ValueError(f"Chunk {summary['chunk']} expects {summary['chunks']} chunks, found {count}")
        if summary["settings"] != settings.to_dict():
            raise ValueError(f"Chunk {summary['chunk']} was called with different settings")
        for field in ("contigs", "plan", "genome_blacklist"):
            if summary[field] != first[field]:
                raise ValueError(f"Chunks disagree about {field}")
    return chunks


def read_chunk(summary: dict, summary_path: str, owned: list[OwnedInterval]) -> Iterator[dict]:
    """Yield a chunk's records, checking each is owned by that chunk and that
    none are missing."""
    records_path = str(Path(summary_path).with_suffix("")) + ".jsonl.gz"
    count = 0
    with gzip.open(records_path, "rt") as handle:
        for line in handle:
            record = json.loads(line)
            if not any(interval.owns(record["contig"], record["position_0based"]) for interval in owned):
                raise ValueError(
                    f"Chunk {summary['chunk']} reported {record['contig']}:{record['position_0based']}, "
                    "outside the ranges it owns"
                )
            count += 1
            yield record
    if count != summary["records"]:
        raise ValueError(f"Chunk {summary['chunk']} has {count} records, its summary says {summary['records']}")


def unique(records: Iterator[dict], contig_order: dict[str, int]) -> Iterator[dict]:
    """Pass records through, failing on a repeated contig/position/strand."""
    previous = None
    for record in records:
        key = record_key(record, contig_order)
        if key == previous:
            raise ValueError(f"Duplicate call at {record['contig']}:{record['position_0based']} {record['strand']}")
        previous = key
        yield record


def stage_rows(rows: Iterator[dict], staged) -> list[float]:
    """Write rows to the staging file; return their control p-values in order."""
    p_values = []
    for row in rows:
        p_value = row["control_fisher_p"]
        if p_value is not None:
            if math.isnan(p_value):
                raise ValueError(f"NaN control p-value on {row['contig']}")
            p_values.append(p_value)
        staged.write(json.dumps(row) + "\n")
    return p_values


def stage_pairs(records: Iterator[dict], settings, sample, control_sample, staged) -> tuple[list[float], int]:
    """Pair each contig's endpoint records and stage the resulting rows."""
    p_values: list[float] = []
    candidates = 0
    for _contig, contig_records in itertools.groupby(records, key=lambda record: record["contig"]):
        rows, contig_candidates = digenome.pair_rows(list(contig_records), settings, sample, control_sample)
        candidates += contig_candidates
        p_values += stage_rows(iter(rows), staged)
    return p_values, candidates


def filtered_rows(staged, q_values: Iterator[float], settings: CallerSettings) -> Iterator[dict]:
    """Read staged rows back, attach the sample-wide q-values, and filter."""
    for line in staged:
        row = json.loads(line)
        if row["control_fisher_p"] is not None:
            row["control_fisher_q"] = next(q_values)
        apply_filters(row, settings)
        yield row
