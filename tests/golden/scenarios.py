"""Synthetic inputs for the golden-output tests.

Each scenario is a treated BAM plus optional control, VCF, and blacklist, and
caller settings named exactly like the pipeline parameters. The scenarios cover
the calling rules in docs/cleavage_algorithm.md: endpoint counting filters,
every nDigenome class, Digenome pairing priority and conflict chains,
overhangs, blacklists, controls, known indels, artifacts, and MAPQ-0
multimappers.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helpers import (  # noqa: E402
    CHUNK_SETTINGS,
    SITE_SETTINGS,
    forward_background,
    make_read,
    reverse_background,
    reverse_read_ending_at,
    write_bam,
    write_vcf,
)

THREE_CONTIGS = [("chr1", 2000), ("chr2", 2000), ("chr3", 2000)]
ONE_CONTIG = [("chr1", 2000)]

@dataclass
class Scenario:
    name: str
    contigs: list[tuple[str, int]]
    settings: dict[str, object]
    treated: list
    control: list | None = None
    vcf_records: list[str] = field(default_factory=list)
    blacklist: str = ""

    def write_inputs(self, directory: Path) -> dict[str, Path | None]:
        """Write this scenario's files and return their paths."""
        directory.mkdir(parents=True, exist_ok=True)
        paths: dict[str, Path | None] = {
            "bam": write_bam(directory / "treated.bam", self.contigs, self.treated),
            "control_bam": None,
            "vcf": None,
            "blacklist": None,
        }
        if self.control is not None:
            paths["control_bam"] = write_bam(
                directory / "control.bam", self.contigs, self.control
            )
        if self.vcf_records:
            paths["vcf"] = write_vcf(
                directory / "variants.vcf.gz", self.contigs, self.vcf_records
            )
        if self.blacklist:
            paths["blacklist"] = directory / "blacklist.bed"
            paths["blacklist"].write_text(self.blacklist)
        return paths


def paired_site(contig: int, position: int, count: int, prefix: str) -> list:
    """MAPQ-0 forward and reverse pileups ending at `position`, with
    background on both strands and one secondary alignment."""
    reads = [
        make_read(f"{prefix}_fwd_{index}", position, contig=contig, mapq=0)
        for index in range(count)
    ]
    reads += [
        reverse_read_ending_at(f"{prefix}_rev_{index}", position, contig=contig, mapq=0)
        for index in range(count)
    ]
    reads += forward_background(prefix, position, 12, contig=contig, mapq=0)
    reads += reverse_background(prefix, position, 12, contig=contig, mapq=0)
    reads.append(
        make_read(f"{prefix}_secondary", position, contig=contig, mapq=0, flag=256)
    )
    return reads


def control_background(contig: int, position: int, prefix: str) -> list:
    """MAPQ-0 control reads covering `position` on both strands."""
    return forward_background(prefix, position, 20, contig=contig, mapq=0) + (
        reverse_background(prefix, position, 20, contig=contig, mapq=0)
    )


def scattered_reads(contig: int, start: int, prefix: str) -> list:
    return [
        make_read(f"{prefix}_{index}", start + index * 3, contig=contig, mapq=0)
        for index in range(12)
    ]


def multi_contig() -> Scenario:
    """Four contigs with chr10 before chr3 in the header, controls, a known
    indel on chr2, and a blacklisted site on chr1."""
    contigs = [("chr1", 2000), ("chr2", 2000), ("chr10", 2000), ("chr3", 2000)]
    treated = paired_site(0, 500, 8, "chr1_masked")
    treated += paired_site(0, 1500, 7, "chr1")
    treated += paired_site(1, 700, 6, "chr2")
    treated += paired_site(2, 900, 8, "chr10")
    treated += [
        make_read(f"chr10_ssb_{index}", 1300, contig=2, mapq=0) for index in range(8)
    ]
    treated += forward_background("chr10_ssb", 1300, 12, contig=2, mapq=0)
    treated += scattered_reads(3, 100, "chr3_scattered")
    control = control_background(0, 500, "c_chr1_masked")
    control += control_background(0, 1500, "c_chr1")
    control += control_background(1, 700, "c_chr2")
    control += control_background(2, 900, "c_chr10")
    control += control_background(2, 1300, "c_chr10_ssb")
    control += scattered_reads(3, 200, "c_chr3_scattered")
    return Scenario(
        name="multi_contig",
        contigs=contigs,
        settings=CHUNK_SETTINGS,
        treated=treated,
        control=control,
        vcf_records=["chr2\t501\t.\tAA\tA\t60\tPASS\t."],
        blacklist="chr1\t490\t510\n",
    )


def conflict_chain() -> Scenario:
    """Competing Digenome pairs every 2 bp, wider than any fixed padding."""
    reads = []
    for position in range(980, 1021, 2):
        reads += [
            make_read(f"chain_fwd_{position}_{index}", position, mapq=0)
            for index in range(5)
        ]
        reads += [
            reverse_read_ending_at(f"chain_rev_{position}_{index}", position, mapq=0)
            for index in range(5)
        ]
    return Scenario(
        name="conflict_chain",
        contigs=THREE_CONTIGS,
        settings=CHUNK_SETTINGS,
        treated=reads,
        vcf_records=["chr2\t501\t.\tAA\tA\t60\tPASS\t."],
    )


def boundary_overhang() -> Scenario:
    """A pair straddling the middle of chr1 with a 2 bp overhang."""
    reads = [make_read(f"bnd_fwd_{index}", 1001, mapq=0) for index in range(8)]
    reads += [
        reverse_read_ending_at(f"bnd_rev_{index}", 999, mapq=0) for index in range(8)
    ]
    reads += forward_background("bnd_f", 1001, 12, mapq=0)
    reads += reverse_background("bnd_r", 999, 12, mapq=0)
    return Scenario(
        name="boundary_overhang",
        contigs=THREE_CONTIGS,
        settings={**CHUNK_SETTINGS, "digenome_overhang": 2},
        treated=reads,
        vcf_records=["chr2\t501\t.\tAA\tA\t60\tPASS\t."],
    )


def single_sites() -> Scenario:
    """One contig with one site per behavior, a partial control, a known
    indel, and a blacklisted reverse endpoint."""
    reads = []
    # 100: forward-only pileup -> SSB.
    reads += [make_read(f"ssb_{index}", 100) for index in range(8)]
    reads += forward_background("ssb", 100, 12)
    # 300: strong reverse partner one base away -> POSSIBLE_DSB / Digenome pair.
    reads += [make_read(f"dsb_fwd_{index}", 300) for index in range(8)]
    reads += [reverse_read_ending_at(f"dsb_rev_{index}", 301) for index in range(8)]
    reads += forward_background("dsb", 300, 8)
    # 500: weak reverse partner -> AMBIGUOUS.
    reads += [make_read(f"amb_fwd_{index}", 500) for index in range(8)]
    reads += [reverse_read_ending_at(f"amb_rev_{index}", 501) for index in range(2)]
    reads += forward_background("amb", 500, 12)
    # 700: 5' soft clips -> HIGH_5P_SOFTCLIP.
    reads += [
        make_read(f"clip_{index}", 700, cigar=[(4, 5), (0, 45)]) for index in range(8)
    ]
    reads += forward_background("clip", 700, 10)
    # 900: insertion near the endpoint plus a VCF indel -> NEARBY_INDEL, KNOWN_INDEL.
    reads += [
        make_read(f"indel_{index}", 900, cigar=[(0, 5), (1, 1), (0, 44)])
        for index in range(8)
    ]
    reads += forward_background("indel", 900, 8)
    # 1100: only reads that must not count.
    for label, flag, mapq in (
        ("dup", 1024, 60),
        ("sec", 256, 60),
        ("sup", 2048, 60),
        ("lowq", 0, 0),
    ):
        reads += [
            make_read(f"{label}_{index}", 1100, flag=flag, mapq=mapq)
            for index in range(8)
        ]
    reads += forward_background("excluded", 1100, 10)
    # 1300: scattered starts, no pileup.
    reads += [make_read(f"random_{index}", 1300 + index * 3) for index in range(30)]
    # 1500: reverse endpoint 1498 is blacklisted.
    reads += [make_read(f"mask_fwd_{index}", 1500) for index in range(8)]
    reads += [reverse_read_ending_at(f"mask_rev_{index}", 1498) for index in range(8)]
    # 1650: exact pair that the control also shows.
    reads += [make_read(f"ctl_fwd_{index}", 1650) for index in range(10)]
    reads += [reverse_read_ending_at(f"ctl_rev_{index}", 1650) for index in range(10)]
    # 1850: forward-only pileup with no background.
    reads += [make_read(f"lone_{index}", 1850) for index in range(8)]

    control = forward_background("c_ssb", 100, 30)
    control += forward_background("c_dsb", 300, 30)
    control += [make_read(f"c_ctl_fwd_{index}", 1650) for index in range(10)]
    control += [
        reverse_read_ending_at(f"c_ctl_rev_{index}", 1650) for index in range(10)
    ]
    return Scenario(
        name="single_sites",
        contigs=ONE_CONTIG,
        settings=SITE_SETTINGS,
        treated=reads,
        control=control,
        vcf_records=["chr1\t904\t.\tAA\tA\t60\tPASS\t."],
        blacklist="chr1\t1498\t1499\n",
    )


def multimappers() -> Scenario:
    """MAPQ-0 primaries count once; secondaries are diagnostic only."""
    reads = [make_read(f"mm_{index}", 100, mapq=0) for index in range(8)]
    reads += [make_read(f"mm_sec_{index}", 100, mapq=0, flag=256) for index in range(8)]
    reads += forward_background("mm", 100, 12)
    reads += [make_read(f"mm_pair_fwd_{index}", 500, mapq=0) for index in range(8)]
    reads += [
        reverse_read_ending_at(f"mm_pair_rev_{index}", 500, mapq=0) for index in range(8)
    ]
    reads += [
        make_read(f"mm_pair_sec_{index}", 500, mapq=0, flag=256) for index in range(8)
    ]
    return Scenario(
        name="multimappers",
        contigs=ONE_CONTIG,
        settings={
            **SITE_SETTINGS,
            "keep_multimappers": True,
            "ndigenome_min_mapq": 0,
            "digenome_min_mapq": 0,
            "cleavage_min_support_mean_mapq": 0.0,
        },
        treated=reads,
    )


def digenome_pairs() -> Scenario:
    """Pairing priority, a low-depth pair, and a pair with an indel."""
    reads = []
    # 500: depth 8 fails the depth cutoff of 10 -> manual review.
    reads += [make_read(f"low_fwd_{index}", 500) for index in range(8)]
    reads += [reverse_read_ending_at(f"low_rev_{index}", 500) for index in range(8)]
    # 800: forward reads carry an insertion; VCF indel nearby.
    reads += [
        make_read(f"ind_fwd_{index}", 800, cigar=[(0, 5), (1, 1), (0, 44)])
        for index in range(8)
    ]
    reads += [reverse_read_ending_at(f"ind_rev_{index}", 800) for index in range(8)]
    # 1000: two reverse candidates; only the 998 pair passes every threshold.
    reads += [make_read(f"pri_fwd_{index}", 1000) for index in range(8)]
    reads += forward_background("pri", 1000, 12)
    reads += [reverse_read_ending_at(f"pri_exact_{index}", 1000) for index in range(8)]
    reads += [reverse_read_ending_at(f"pri_near_{index}", 998) for index in range(8)]
    return Scenario(
        name="digenome_pairs",
        contigs=ONE_CONTIG,
        settings={
            **SITE_SETTINGS,
            "digenome_depth_cutoff": 10,
            "digenome_pair_score_cutoff": 0.5,
        },
        treated=reads,
        vcf_records=["chr1\t805\t.\tAA\tA\t60\tPASS\t."],
    )


def all_scenarios() -> list[Scenario]:
    return [
        multi_contig(),
        conflict_chain(),
        boundary_overhang(),
        single_sites(),
        multimappers(),
        digenome_pairs(),
    ]
