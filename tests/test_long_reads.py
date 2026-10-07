"""Long-read calling: both aligned ends of each read count (decision
"Long-read input counts both ends of each read" in docs/decisions.md).

A synthetic 10 kb contig with reads of about 2 kb. A cut between 4999 and 5000
leaves molecules whose left end is 5000 and molecules whose right end is 4999,
and reads of either orientation can sequence them. The far ends of the reads
are staggered so that only the cut stacks.
"""

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
from helpers import make_read, make_settings, run_caller, write_bam  # noqa: E402

from cleavage.bam import ReadEnd, measure_site, read_ends  # noqa: E402
from cleavage.output import format_value  # noqa: E402
from cleavage.stats import fisher_exact_two_sided  # noqa: E402

CONTIG = [("chr1", 10000)]
CUT = 5000
FORWARD, REVERSE, BOTH = (False,), (True,), (False, True)
CLIP = (4, 30)


def long_read(name: str, start: int, end: int, reverse: bool) -> pysam.AlignedSegment:
    """A read aligned to [start, end) with no clips."""
    return make_read(name, start, reverse=reverse, cigar=[(0, end - start)])


def stagger(index: int, reverse: bool) -> int:
    """Distinct offsets for the far ends of reads built together."""
    return 10 * index + (5 if reverse else 0)


def cut_reads(prefix: str, per_orientation: int, orientations=BOTH) -> list:
    """Reads ending at the cut from each side: left ends at 5000, right ends at 4999."""
    reads = []
    for reverse in orientations:
        for index in range(per_orientation):
            length = 2000 + stagger(index, reverse)
            reads.append(long_read(f"{prefix}_right_{reverse}_{index}", CUT, CUT + length, reverse))
            reads.append(long_read(f"{prefix}_left_{reverse}_{index}", CUT - length, CUT, reverse))
    return reads


def spanning_reads(prefix: str, per_orientation: int, orientations=BOTH) -> list:
    """Uncut molecules covering both sides of the cut."""
    return [
        long_read(
            f"{prefix}_span_{reverse}_{index}",
            CUT - 1000 - stagger(index, reverse),
            CUT + 1000 + stagger(index, reverse),
            reverse,
        )
        for reverse in orientations
        for index in range(per_orientation)
    ]


def long_read_settings(analysis: str):
    """What main.nf passes with --long_reads: clip and indel limits at 1.0."""
    return make_settings(
        analysis,
        long_reads=True,
        cleavage_max_softclip_fraction=1.0,
        cleavage_max_indel_fraction=1.0,
    )


class LongReadEndTests(unittest.TestCase):
    def test_digenome_counts_left_ends_as_plus_and_right_ends_as_minus(self) -> None:
        settings = long_read_settings("digenome")
        for reverse in (False, True):
            self.assertEqual(
                read_ends(long_read("read", 1000, 3000, reverse), settings),
                [ReadEnd(1000, "+", True), ReadEnd(2999, "-", False)],
            )

    def test_ndigenome_keeps_both_ends_on_the_read_strand(self) -> None:
        settings = long_read_settings("ndigenome")
        self.assertEqual(
            read_ends(long_read("forward", 1000, 3000, False), settings),
            [ReadEnd(1000, "+", True), ReadEnd(2999, "+", False)],
        )
        self.assertEqual(
            read_ends(long_read("reverse", 1000, 3000, True), settings),
            [ReadEnd(1000, "-", True), ReadEnd(2999, "-", False)],
        )

    def test_one_base_alignment_has_one_end(self) -> None:
        read = long_read("tiny", 7000, 7001, False)
        for analysis in ("digenome", "ndigenome"):
            self.assertEqual(read_ends(read, long_read_settings(analysis)), [ReadEnd(7000, "+", True)])


class LongReadDigenomeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def call(self, treated: list, control: list | None = None, chunks: int = 1, name: str = "out") -> list[dict]:
        inputs = {
            "bam": write_bam(self.tmp / f"{name}_treated.bam", CONTIG, treated),
            "control_bam": write_bam(self.tmp / f"{name}_control.bam", CONTIG, control) if control else None,
            "vcf": None,
            "blacklist": None,
        }
        prefix, _qc = run_caller(long_read_settings("digenome"), inputs, self.tmp / name, chunks)
        with open(f"{prefix}.digenome.all.tsv", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_cut_ends_from_both_orientations_pair(self) -> None:
        rows = self.call(cut_reads("cut", 4) + spanning_reads("bg", 3))
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(
            (row["forward_position_0based"], row["reverse_position_0based"]),
            (str(CUT), str(CUT - 1)),
        )
        # 8 reads of either orientation end at each side of the cut; each
        # side's depth adds the 6 spanning reads, whatever their orientation.
        self.assertEqual((row["forward_endpoint_count"], row["forward_depth"]), ("8", "14"))
        self.assertEqual((row["reverse_endpoint_count"], row["reverse_depth"]), ("8", "14"))
        self.assertEqual(row["digenome_pair_score"], format_value((8 / 14) * (8 / 14) * (8 + 8) / 4))
        self.assertEqual(row["tier"], "high_confidence")

    def test_combined_depth_counts_each_spanning_read_once(self) -> None:
        rows = self.call(cut_reads("cut", 4) + spanning_reads("bg", 3), control=spanning_reads("ctl", 5))
        row = rows[0]
        # Reads covering 4999 or 5000: 8 + 8 cut reads and 6 spanning reads, not 14 + 14.
        self.assertEqual((row["combined_endpoint_count"], row["combined_depth"]), ("16", "22"))
        self.assertEqual(row["combined_fraction"], format_value(16 / 22))
        # The control's 10 spanning reads cover both endpoints and count once.
        self.assertEqual((row["control_forward_depth"], row["control_reverse_depth"]), ("10", "10"))
        self.assertEqual((row["control_combined_endpoint_count"], row["control_combined_depth"]), ("0", "10"))
        self.assertEqual(row["control_fisher_p"], format_value(fisher_exact_two_sided(16, 22 - 16, 0, 10)))

    def test_pair_across_a_chunk_boundary_matches_one_chunk(self) -> None:
        # Two chunks own [0, 5000) and [5000, 10000), so the forward endpoint
        # and its partner belong to different chunks.
        treated = cut_reads("cut", 4) + spanning_reads("bg", 3)
        control = spanning_reads("ctl", 5)
        one = self.call(treated, control, chunks=1, name="one")
        self.assertEqual(len(one), 1)
        for chunks in (2, 3):
            with self.subTest(chunks=chunks):
                self.assertEqual(self.call(treated, control, chunks=chunks, name=f"c{chunks}"), one)

    def test_clip_metric_uses_the_clip_at_the_endpoint(self) -> None:
        # Two reads per endpoint are clipped at the end that forms it and one
        # at its far end. The orientations make the 5' clip disagree.
        reads = cut_reads("cut", 4) + [
            make_read("plus_end_clip_0", CUT, reverse=True, cigar=[CLIP, (0, 2000)]),
            make_read("plus_end_clip_1", CUT, reverse=True, cigar=[CLIP, (0, 2001)]),
            make_read("plus_far_clip", CUT, cigar=[(0, 2002), CLIP]),
            make_read("minus_end_clip_0", CUT - 2000, cigar=[(0, 2000), CLIP]),
            make_read("minus_end_clip_1", CUT - 2001, cigar=[(0, 2001), CLIP]),
            make_read("minus_far_clip", CUT - 2002, reverse=True, cigar=[CLIP, (0, 2002)]),
        ]
        bam = write_bam(self.tmp / "clips.bam", CONTIG, reads)
        settings = long_read_settings("digenome")
        with pysam.AlignmentFile(str(bam), "rb") as handle:
            plus = measure_site(handle, "chr1", CUT, "+", settings)
            minus = measure_site(handle, "chr1", CUT - 1, "-", settings)
        self.assertEqual((plus.endpoint_count, plus.softclipped_count), (11, 2))
        self.assertEqual((minus.endpoint_count, minus.softclipped_count), (11, 2))


class LongReadNDigenomeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def call(self, reads: list) -> list[tuple[str, str, str]]:
        inputs = {"bam": write_bam(self.tmp / "treated.bam", CONTIG, reads), "control_bam": None, "vcf": None, "blacklist": None}
        prefix, _qc = run_caller(long_read_settings("ndigenome"), inputs, self.tmp / "out", chunks=2)
        with open(f"{prefix}.ndigenome.all.tsv", newline="") as handle:
            return [
                (row["position_0based"], row["strand"], row["signal_classification"])
                for row in csv.DictReader(handle, delimiter="\t")
            ]

    def nick_reads(self) -> list:
        """A nick in the top strand of a strand-keeping library: forward reads
        end at 4999 or start at 5000, and reverse reads span the site."""
        return (
            cut_reads("nick", 6, FORWARD)
            + spanning_reads("fwd", 4, FORWARD)
            + spanning_reads("rev", 6, REVERSE)
        )

    def test_nick_gives_two_same_strand_ssb_rows(self) -> None:
        self.assertEqual(
            self.call(self.nick_reads()),
            [(str(CUT - 1), "+", "SSB"), (str(CUT), "+", "SSB")],
        )

    def test_ends_on_both_strands_are_possible_dsb(self) -> None:
        reads = self.nick_reads() + cut_reads("dsb", 6, REVERSE)
        self.assertEqual(
            self.call(reads),
            [
                (str(CUT - 1), "+", "POSSIBLE_DSB"),
                (str(CUT - 1), "-", "POSSIBLE_DSB"),
                (str(CUT), "+", "POSSIBLE_DSB"),
                (str(CUT), "-", "POSSIBLE_DSB"),
            ],
        )


if __name__ == "__main__":
    unittest.main()
