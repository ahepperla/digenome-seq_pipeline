"""Call one chunk of one sample: chunk `index` of `count`.

Every chunk task recomputes the same deterministic plan from the BAM index
and the blacklist, then measures only the endpoints its slot owns. It writes:
- <prefix>.jsonl.gz: nDigenome rows, or Digenome endpoint records, sorted by
  contig (BAM header order), position, and strand. No q-values or filters yet.
- <prefix>.json: a summary the finalizer uses to check that every chunk ran
  with the same inputs and that ownership covers the genome.
"""

from __future__ import annotations

import gzip
import json
from contextlib import ExitStack
from pathlib import Path

from . import digenome, ndigenome
from .bam import check_layout, mapped_contigs, open_bam
from .output import record_key
from .regions import load_blacklist, plan_chunks, plan_to_json
from .settings import CallerSettings
from .vcf import open_vcf


def call_chunk(
    settings: CallerSettings,
    *,
    bam_path: str,
    sample: str,
    index: int,
    count: int,
    out_prefix: str,
    control_path: str | None = None,
    control_sample: str = "",
    vcf_path: str | None = None,
    blacklist_path: str | None = None,
) -> dict:
    if not 0 <= index < count:
        raise ValueError(f"Chunk index {index} is outside 0..{count - 1}")
    # Short-read nDigenome needs paired-end libraries; long reads are single molecules.
    require_paired = settings.analysis == "ndigenome" and not settings.long_reads
    with ExitStack() as stack:
        bam = stack.enter_context(open_bam(bam_path, "BAM"))
        check_layout(bam, settings.min_mapq, require_paired, "BAM")
        lengths = dict(zip(bam.references, bam.lengths))
        blacklist = load_blacklist(Path(blacklist_path), lengths) if blacklist_path else None
        contigs = mapped_contigs(bam)
        if not contigs:
            raise ValueError("BAM index reports no mapped alignments")
        plan = plan_chunks(contigs, count, blacklist)
        intervals = plan[index]

        control_bam = None
        if control_path:
            control_bam = stack.enter_context(open_bam(control_path, "Control BAM"))
            if control_bam.references != bam.references or control_bam.lengths != bam.lengths:
                raise ValueError("Control BAM reference names and lengths do not match the treated BAM")
            check_layout(control_bam, settings.min_mapq, require_paired, "Control BAM")
        vcf = None
        if vcf_path:
            analyzed = list(dict.fromkeys(interval.contig for interval in intervals))
            vcf = stack.enter_context(open_vcf(vcf_path, lengths, analyzed))

        if settings.analysis == "ndigenome":
            records, candidates = ndigenome.call_chunk(
                bam, control_bam, vcf, settings, intervals, blacklist, sample, control_sample
            )
        else:
            records = digenome.endpoint_records(bam, control_bam, vcf, settings, intervals, blacklist)
            candidates = None  # Digenome candidates are pairs, counted by the finalizer.

    order = {contig.name: position for position, contig in enumerate(contigs)}
    records.sort(key=lambda record: record_key(record, order))
    with gzip.open(f"{out_prefix}.jsonl.gz", "wt", compresslevel=1) as handle:
        for record in records:
            handle.write(json.dumps(record) + "\n")
    summary = {
        "sample": sample,
        "analysis": settings.analysis,
        "chunk": index,
        "chunks": count,
        "settings": settings.to_dict(),
        "contigs": [[contig.name, contig.length] for contig in contigs],
        "plan": plan_to_json(plan),
        "genome_blacklist": blacklist.provenance() if blacklist else None,
        "candidates": candidates,
        "records": len(records),
    }
    Path(f"{out_prefix}.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary
