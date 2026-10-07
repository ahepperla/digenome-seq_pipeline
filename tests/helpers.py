"""Build small synthetic BAM, VCF, and BED inputs for the tests.

Reads are 50 bp, paired-flagged, and default to MAPQ 60 with NM 0. A forward
read's 5' endpoint is its start; a reverse read's is its last aligned base.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from cleavage.call import call_chunk  # noqa: E402
from cleavage.finalize import finalize  # noqa: E402
from cleavage.settings import CallerSettings, settings_from_dict  # noqa: E402

# Permissive settings used by the chunking tests.
CHUNK_SETTINGS = {
    "keep_multimappers": True,
    "long_reads": False,
    "ndigenome_min_count": 5,
    "ndigenome_min_fraction": 0.20,
    "ndigenome_min_mapq": 0,
    "ndigenome_opposite_window": 2,
    "ndigenome_ambiguous_min_count": 2,
    "ndigenome_ambiguous_min_fraction": 0.05,
    "digenome_overhang": 0,
    "digenome_pair_window": 2,
    "digenome_min_mapq": 0,
    "digenome_forward_cutoff": 4,
    "digenome_reverse_cutoff": 4,
    "digenome_depth_cutoff": 4,
    "digenome_fraction_cutoff": 0.20,
    "digenome_pair_score_cutoff": 0.20,
    "cleavage_artifact_window": 10,
    "cleavage_max_softclip_fraction": 0.20,
    "cleavage_max_indel_fraction": 0.20,
    "cleavage_min_support_mean_mapq": 0.0,
    "cleavage_control_min_depth": 1,
    "cleavage_control_max_fraction": 0.05,
    "cleavage_control_min_fold": 3.0,
    "cleavage_control_max_q": 0.05,
}

# Settings used by the single-site caller tests.
SITE_SETTINGS = {
    **CHUNK_SETTINGS,
    "keep_multimappers": False,
    "ndigenome_min_mapq": 1,
    "ndigenome_opposite_window": 5,
    "digenome_min_mapq": 1,
    "digenome_pair_score_cutoff": 1.0,
    "cleavage_min_support_mean_mapq": 10.0,
}


def make_settings(analysis: str, base: dict | None = None, **overrides) -> CallerSettings:
    """Caller settings for tests: SITE_SETTINGS unless `base` is given."""
    return settings_from_dict({**(base or SITE_SETTINGS), "analysis": analysis, **overrides})


READ_LENGTH = 50
QUERY_OPERATIONS = (0, 1, 4, 7, 8)  # M, I, S, =, X consume query bases


def make_read(
    name: str,
    start: int,
    *,
    contig: int = 0,
    reverse: bool = False,
    mapq: int = 60,
    flag: int = 0,
    cigar: list[tuple[int, int]] | None = None,
    nm: int = 0,
) -> pysam.AlignedSegment:
    """Return one aligned read; `flag` adds bits such as 256 or 1024."""
    cigar = cigar or [(0, READ_LENGTH)]
    query_length = sum(
        length for operation, length in cigar if operation in QUERY_OPERATIONS
    )
    read = pysam.AlignedSegment()
    read.query_name = name
    read.query_sequence = "A" * query_length
    read.flag = 1 | (16 if reverse else 0) | flag
    read.reference_id = contig
    read.reference_start = start
    read.mapping_quality = mapq
    read.cigartuples = cigar
    read.next_reference_id = contig
    read.next_reference_start = start + 100
    read.template_length = 150
    read.query_qualities = pysam.qualitystring_to_array("I" * query_length)
    read.set_tag("NM", nm)
    return read


def reverse_read_ending_at(
    name: str, endpoint: int, **options
) -> pysam.AlignedSegment:
    """Return a 50M reverse read whose 5' endpoint is `endpoint`."""
    return make_read(name, endpoint - READ_LENGTH + 1, reverse=True, **options)


def forward_background(
    prefix: str, position: int, count: int, **options
) -> list[pysam.AlignedSegment]:
    """Forward reads that cover `position` without ending at it."""
    return [
        make_read(f"{prefix}_fwd_bg_{position}_{index}", position - 20 - index, **options)
        for index in range(count)
    ]


def reverse_background(
    prefix: str, position: int, count: int, **options
) -> list[pysam.AlignedSegment]:
    """Reverse reads that cover `position` without ending at it."""
    return [
        make_read(
            f"{prefix}_rev_bg_{position}_{index}",
            position - 48 + index,
            reverse=True,
            **options,
        )
        for index in range(count)
    ]


def write_bam(
    path: Path,
    contigs: list[tuple[str, int]],
    reads: list[pysam.AlignedSegment],
) -> Path:
    """Write a coordinate-sorted, indexed BAM."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": length} for name, length in contigs],
    }
    ordered = sorted(
        reads,
        key=lambda read: (
            read.reference_id,
            read.reference_start,
            read.flag,
            read.query_name,
        ),
    )
    with pysam.AlignmentFile(str(path), "wb", header=header) as output:
        for read in ordered:
            output.write(read)
    pysam.index(str(path))
    return path


def write_vcf(
    path: Path,
    contigs: list[tuple[str, int]],
    records: list[str],
) -> Path:
    """Write a bgzipped, tabix-indexed VCF. `path` must end in .vcf.gz."""
    plain = path.with_suffix("")
    lines = ["##fileformat=VCFv4.2"]
    lines += [f"##contig=<ID={name},length={length}>" for name, length in contigs]
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    lines += records
    plain.write_text("\n".join(lines) + "\n")
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset="vcf", force=True)
    plain.unlink()
    return path


def run_caller(
    settings,
    inputs: dict[str, Path | None],
    out_dir: Path,
    chunks: int,
    sample: str = "Sample",
    control_sample: str = "Control",
) -> tuple[str, dict]:
    """Run every chunk and the finalizer the way main.nf does. Returns the
    output prefix and the QC report."""
    out_dir.mkdir(parents=True, exist_ok=True)
    control = inputs.get("control_bam")
    summaries = []
    for index in range(chunks):
        chunk_prefix = out_dir / f"chunk_{index:03d}"
        call_chunk(
            settings,
            bam_path=str(inputs["bam"]),
            sample=sample,
            index=index,
            count=chunks,
            out_prefix=str(chunk_prefix),
            control_path=str(control) if control else None,
            control_sample=control_sample if control else "",
            vcf_path=str(inputs["vcf"]) if inputs.get("vcf") else None,
            blacklist_path=str(inputs["blacklist"]) if inputs.get("blacklist") else None,
        )
        summaries.append(f"{chunk_prefix}.json")
    prefix = str(out_dir / sample)
    qc = finalize(settings, summaries, sample, control_sample if control else "", prefix)
    return prefix, qc
