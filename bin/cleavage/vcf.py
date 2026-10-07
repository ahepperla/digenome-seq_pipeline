"""Known indels from the optional, indexed variant VCF."""

from __future__ import annotations

from pathlib import Path

import pysam


def open_vcf(path: str, contig_lengths: dict[str, int], contigs: list[str]) -> pysam.VariantFile:
    """Open a bgzipped, indexed VCF that declares every analyzed contig.

    A missing contig or a different declared length means a different
    reference build or naming (chr1 versus 1); that stops the run rather than
    silently losing annotations.
    """
    if not any(Path(f"{path}{suffix}").is_file() for suffix in (".tbi", ".csi")):
        raise ValueError(f"Variant VCF must have a .tbi or .csi index: {path}")
    vcf = pysam.VariantFile(str(path))
    declared = vcf.header.contigs
    missing = [contig for contig in contigs if contig not in declared]
    if missing:
        preview = ", ".join(missing[:10])
        if len(missing) > 10:
            preview += f", ... ({len(missing) - 10} more)"
        raise ValueError(
            f"Variant VCF is missing analyzed BAM contig(s): {preview}. "
            "Check reference builds and contig naming."
        )
    mismatched = [
        f"{contig} (BAM {contig_lengths[contig]}, VCF {declared[contig].length})"
        for contig in contigs
        if declared[contig].length is not None
        and declared[contig].length != contig_lengths[contig]
    ]
    if mismatched:
        raise ValueError(
            "Variant VCF contig length does not match the BAM: " + ", ".join(mismatched[:10])
        )
    return vcf


def known_indels(
    vcf: pysam.VariantFile | None, contig: str, positions: list[int], window: int
) -> list[str]:
    """Indel records within `window` bases of any of `positions`, as sorted,
    de-duplicated "contig:pos:REF>ALT[,ALT]" strings (1-based VCF positions)."""
    if vcf is None:
        return []
    found = set()
    for position in positions:
        for record in vcf.fetch(contig, max(0, position - window), position + window + 1):
            indel_alts = [alt for alt in record.alts or () if alt and len(alt) != len(record.ref)]
            if indel_alts:
                found.add(f"{contig}:{record.pos}:{record.ref}>{','.join(indel_alts)}")
    return sorted(found)
