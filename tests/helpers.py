"""Build small synthetic BAM, VCF, and BED inputs for the tests.

Reads are 50 bp, paired-flagged, and default to MAPQ 60 with NM 0. A forward
read's 5' endpoint is its start; a reverse read's is its last aligned base.
"""

from __future__ import annotations

from pathlib import Path

import pysam

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
