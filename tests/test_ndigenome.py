"""Tests for nDigenome endpoint calling.

Tests the ndigenome module using the new cleavage package, with the same
read layouts and expectations from the old tests but adapted to the new
column names and structure.
"""

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

try:
    import pysam
except ImportError:
    pysam = None

sys.path.insert(0, str(Path(__file__).resolve().parent))
from helpers import (  # noqa: E402
    forward_background,
    make_read,
    make_settings,
    reverse_background,
    reverse_read_ending_at,
    run_caller,
    write_bam,
    write_vcf,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from cleavage.finalize import MULTIMAPPER_WARNING  # noqa: E402
from cleavage.bam import measure_site  # noqa: E402
from cleavage.ndigenome import strongest_opposite  # noqa: E402
from cleavage.output import format_value  # noqa: E402
from cleavage.stats import fisher_exact_two_sided, benjamini_hochberg  # noqa: E402


@unittest.skipIf(pysam is None, "pysam is not installed")
class NDigenomeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def read_rows(self, prefix: str) -> list[dict[str, str]]:
        path = self.tmp / f"{prefix}.ndigenome.all.tsv"
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def read_tier_rows(self, prefix: str, tier: str) -> list[dict[str, str]]:
        path = self.tmp / f"{prefix}.ndigenome.{tier}.tsv"
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_true_single_strand_signal_passes(self) -> None:
        """True single-strand signal passes: 1 row, signal SSB, filter_status PASS."""
        reads = [make_read(f"signal_{index}", 100) for index in range(8)]
        reads += forward_background("bg", 100, 12)
        bam = write_bam(self.tmp / "ssb.bam", [("chr1", 2000)], reads)
        settings = make_settings("ndigenome")
        prefix, qc = run_caller(settings, {"bam": bam, "control_bam": None, "vcf": None, "blacklist": None}, self.tmp, chunks=1, sample="ssb")
        rows = self.read_rows("ssb")
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["signal_classification"], "SSB")
        self.assertEqual(rows[0]["filter_status"], "PASS")
        self.assertEqual(rows[0]["tier"], "high_confidence")
        self.assertEqual(len(self.read_tier_rows("ssb", "high_confidence")), 1)
        self.assertEqual(self.read_tier_rows("ssb", "manual_review"), [])
        self.assertEqual(self.read_tier_rows("ssb", "artifact"), [])

    def test_opposite_strand_signal_is_possible_dsb(self) -> None:
        """Opposite-strand signal → POSSIBLE_DSB, FILTERED."""
        reads = [make_read(f"forward_{index}", 100) for index in range(8)]
        reads += [reverse_read_ending_at(f"reverse_{index}", 101) for index in range(8)]
        reads += forward_background("bg", 100, 8)
        bam = write_bam(self.tmp / "dsb.bam", [("chr1", 2000)], reads)
        settings = make_settings("ndigenome")
        prefix, qc = run_caller(settings, {"bam": bam, "control_bam": None, "vcf": None, "blacklist": None}, self.tmp, chunks=1, sample="dsb")
        rows = self.read_rows("dsb")
        self.assertEqual([(row["position_0based"], row["strand"]) for row in rows], [("100", "+"), ("101", "-")])
        for row in rows:
            self.assertEqual(row["signal_classification"], "POSSIBLE_DSB")
            self.assertEqual(row["filter_status"], "FILTERED")
            self.assertEqual(row["filter_reasons"], "POSSIBLE_DSB")

    def test_blacklisted_opposite_endpoint(self) -> None:
        """Blacklisted opposite endpoint (BED chr1 102-103) → 1 row at 100, SSB."""
        reads = [make_read(f"forward_{index}", 100) for index in range(8)]
        reads += [reverse_read_ending_at(f"reverse_{index}", 102) for index in range(8)]
        reads += forward_background("bg", 100, 8)
        bam = write_bam(self.tmp / "blacklisted_opposite.bam", [("chr1", 2000)], reads)
        blacklist = self.tmp / "blacklist.bed"
        blacklist.write_text("chr1\t102\t103\n")
        settings = make_settings("ndigenome")
        prefix, qc = run_caller(settings, {"bam": bam, "control_bam": None, "vcf": None, "blacklist": blacklist}, self.tmp, chunks=1, sample="blacklisted")
        rows = self.read_rows("blacklisted")
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["position_0based"], "100")
        self.assertEqual(rows[0]["signal_classification"], "SSB")
        self.assertEqual(rows[0]["opposite_position_0based"], "")
        self.assertEqual(rows[0]["opposite_count"], "0")

    def test_opposite_ranking_prefers_threshold_passing_fraction(self) -> None:
        """An opposite endpoint passing count and fraction outranks one with
        more reads whose fraction is below the threshold."""
        reads = [make_read(f"forward_{index}", 100) for index in range(8)]
        reads += [reverse_read_ending_at(f"reverse_100_{index}", 100) for index in range(10)]
        reads += [reverse_read_ending_at(f"reverse_101_{index}", 101) for index in range(11)]
        # Reverse reads covering 101 but not 100 dilute 101's fraction to 11/111.
        reads += [make_read(f"reverse_background_{index}", 101, reverse=True) for index in range(100)]
        bam = write_bam(self.tmp / "opposite_ranking.bam", [("chr1", 2000)], reads)
        settings = make_settings("ndigenome", ndigenome_min_count=10)
        with pysam.AlignmentFile(str(bam), "rb") as bam_handle:
            position, metrics = strongest_opposite(bam_handle, "chr1", 100, "+", settings, None)
            competitor = measure_site(bam_handle, "chr1", 101, "-", 10, 1)
        self.assertEqual(position, 100)
        self.assertEqual((metrics.endpoint_count, metrics.strand_depth), (10, 21))
        self.assertEqual((competitor.endpoint_count, competitor.strand_depth), (11, 111))

    def test_weak_opposite_signal_is_ambiguous(self) -> None:
        """Weak opposite → AMBIGUOUS, tier manual_review."""
        reads = [make_read(f"forward_{index}", 100) for index in range(8)]
        reads += [reverse_read_ending_at(f"weak_reverse_{index}", 101) for index in range(2)]
        reads += reverse_background("bg", 101, 48)
        bam = write_bam(self.tmp / "ambiguous.bam", [("chr1", 2000)], reads)
        settings = make_settings("ndigenome")
        prefix, qc = run_caller(settings, {"bam": bam, "control_bam": None, "vcf": None, "blacklist": None}, self.tmp, chunks=1, sample="ambiguous")
        rows = self.read_rows("ambiguous")
        self.assertEqual(rows[0]["signal_classification"], "AMBIGUOUS")
        self.assertEqual(len(self.read_tier_rows("ambiguous", "manual_review")), 1)
        self.assertEqual(self.read_tier_rows("ambiguous", "artifact"), [])

    def test_duplicate_secondary_supplementary_and_low_mapq_reads_excluded(self) -> None:
        """Duplicate, secondary, supplementary, and MAPQ-0 reads excluded."""
        reads = [make_read(f"duplicate_{index}", 100, flag=1024) for index in range(8)]
        reads += [make_read(f"secondary_{index}", 100, flag=256) for index in range(8)]
        reads += [make_read(f"supplementary_{index}", 100, flag=2048) for index in range(8)]
        reads += [make_read(f"low_mapq_{index}", 100, mapq=0) for index in range(8)]
        reads += forward_background("bg", 100, 10)
        bam = write_bam(self.tmp / "filtered.bam", [("chr1", 2000)], reads)
        settings = make_settings("ndigenome")
        prefix, qc = run_caller(settings, {"bam": bam, "control_bam": None, "vcf": None, "blacklist": None}, self.tmp, chunks=1, sample="filtered")
        self.assertEqual(qc["candidates_before_filters"], 0)

    def test_mapq_zero_primaries_counted_once_with_multimappers(self) -> None:
        """MAPQ-0 primaries counted once with keep_multimappers."""
        reads = [make_read(f"primary_{index}", 100, mapq=0) for index in range(8)]
        reads += [make_read(f"secondary_{index}", 100, flag=256, mapq=0) for index in range(8)]
        reads += forward_background("bg", 100, 12)
        bam = write_bam(self.tmp / "multimappers.bam", [("chr1", 2000)], reads)
        settings = make_settings("ndigenome", ndigenome_min_mapq=0, cleavage_min_support_mean_mapq=0.0, keep_multimappers=True)
        prefix, qc = run_caller(settings, {"bam": bam, "control_bam": None, "vcf": None, "blacklist": None}, self.tmp, chunks=1, sample="multimappers")
        rows = self.read_rows("multimappers")
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["endpoint_count"], "8")
        self.assertEqual(rows[0]["secondary_endpoint_count"], "8")
        self.assertEqual(rows[0]["filter_status"], "PASS")
        self.assertEqual(qc["warnings"], [MULTIMAPPER_WARNING])

    def test_random_endpoints_do_not_create_candidates(self) -> None:
        """Random endpoints → 0 rows."""
        reads = [make_read(f"random_{index}", 100 + index * 3) for index in range(30)]
        bam = write_bam(self.tmp / "random.bam", [("chr1", 2000)], reads)
        settings = make_settings("ndigenome")
        prefix, qc = run_caller(settings, {"bam": bam, "control_bam": None, "vcf": None, "blacklist": None}, self.tmp, chunks=1, sample="random")
        self.assertEqual(qc["candidates_before_filters"], 0)

    def test_five_prime_soft_clipping_creates_artifact_risk(self) -> None:
        """5' soft clipping → tier artifact, filter_reasons contains HIGH_5P_SOFTCLIP."""
        reads = [make_read(f"softclip_{index}", 200, cigar=[(4, 5), (0, 45)]) for index in range(8)]
        reads += forward_background("bg", 200, 10)
        bam = write_bam(self.tmp / "softclip.bam", [("chr1", 2000)], reads)
        settings = make_settings("ndigenome")
        prefix, qc = run_caller(settings, {"bam": bam, "control_bam": None, "vcf": None, "blacklist": None}, self.tmp, chunks=1, sample="softclip")
        row = self.read_rows("softclip")[0]
        self.assertEqual(row["tier"], "artifact")
        self.assertEqual(row["filter_reasons"], "HIGH_5P_SOFTCLIP")
        self.assertEqual(self.read_tier_rows("softclip", "manual_review"), [])
        self.assertEqual(len(self.read_tier_rows("softclip", "artifact")), 1)

    def test_nearby_indel_and_vcf_annotated(self) -> None:
        """Nearby insertion + VCF → filter_reasons exactly "NEARBY_INDEL;KNOWN_INDEL"."""
        reads = [make_read(f"indel_{index}", 300, cigar=[(0, 5), (1, 1), (0, 44)]) for index in range(8)]
        reads += forward_background("bg", 300, 8)
        bam = write_bam(self.tmp / "indel.bam", [("chr1", 2000)], reads)
        vcf = write_vcf(self.tmp / "variants.vcf.gz", [("chr1", 2000)], ["chr1\t304\t.\tAA\tA\t60\tPASS\t."])
        settings = make_settings("ndigenome")
        prefix, qc = run_caller(settings, {"bam": bam, "control_bam": None, "vcf": vcf, "blacklist": None}, self.tmp, chunks=1, sample="indel")
        row = self.read_rows("indel")[0]
        self.assertEqual(row["filter_reasons"], "NEARBY_INDEL;KNOWN_INDEL")
        self.assertEqual(row["known_indel_overlap"], "chr1:304:AA>A")

    def test_unindexed_vcf_raises_error(self) -> None:
        """Unindexed VCF → run_caller raises ValueError."""
        reads = [make_read(f"signal_{index}", 350) for index in range(8)]
        reads += forward_background("bg", 350, 8)
        bam = write_bam(self.tmp / "unindexed_vcf.bam", [("chr1", 2000)], reads)
        vcf = self.tmp / "unindexed.vcf.gz"
        vcf.write_text("not-empty\n")
        settings = make_settings("ndigenome")
        with self.assertRaisesRegex(ValueError, "must have a .tbi or .csi"):
            run_caller(settings, {"bam": bam, "control_bam": None, "vcf": vcf, "blacklist": None}, self.tmp, chunks=1, sample="bad")

    def test_matched_control_statistics_and_q_values(self) -> None:
        """Matched control: control_status MATCHED_CONTROL, control_fisher_q < 0.05."""
        treated_reads = [make_read(f"treated_signal_{index}", 400) for index in range(10)]
        treated_reads += forward_background("bg", 400, 20)
        control_reads = forward_background("control", 400, 30)
        treated = write_bam(self.tmp / "treated.bam", [("chr1", 2000)], treated_reads)
        control = write_bam(self.tmp / "control.bam", [("chr1", 2000)], control_reads)
        settings = make_settings("ndigenome")
        prefix, qc = run_caller(settings, {"bam": treated, "control_bam": control, "vcf": None, "blacklist": None}, self.tmp, chunks=1, sample="controlled", control_sample="Control")
        row = self.read_rows("controlled")[0]
        self.assertEqual(row["control_status"], "MATCHED_CONTROL")
        # One controlled row, so q equals p: Fisher on 10 of 30 treated vs 0 of 30 control.
        expected_p = fisher_exact_two_sided(10, 20, 0, 30)
        self.assertEqual(row["control_fisher_p"], format_value(expected_p))
        self.assertEqual(row["control_fisher_q"], format_value(benjamini_hochberg([expected_p])[0]))
        self.assertEqual(row["filter_status"], "PASS")

    def test_zero_control_depth_insufficient(self) -> None:
        """Zero control depth → control_status INSUFFICIENT_CONTROL_COVERAGE."""
        treated_reads = [make_read(f"treated_signal_{index}", 450) for index in range(10)]
        treated_reads += forward_background("bg", 450, 20)
        control_reads = [make_read("control_elsewhere", 100)]
        treated = write_bam(self.tmp / "zero_depth_treated.bam", [("chr1", 2000)], treated_reads)
        control = write_bam(self.tmp / "zero_depth_control.bam", [("chr1", 2000)], control_reads)
        settings = make_settings("ndigenome")
        prefix, qc = run_caller(settings, {"bam": treated, "control_bam": control, "vcf": None, "blacklist": None}, self.tmp, chunks=1, sample="zero_control_depth", control_sample="Control")
        row = self.read_rows("zero_control_depth")[0]
        self.assertEqual(row["control_status"], "INSUFFICIENT_CONTROL_COVERAGE")
        self.assertEqual(row["control_depth"], "0")
        self.assertEqual(row["control_fold_enrichment"], "")
        self.assertEqual(row["control_fisher_q"], "")
        self.assertEqual(row["filter_reasons"], "INSUFFICIENT_CONTROL_COVERAGE")

    def test_vcf_contig_mismatch_raises_error(self) -> None:
        """VCF declaring contig "1" instead of "chr1" → ValueError."""
        reads = [make_read(f"signal_{index}", 475) for index in range(8)]
        reads += forward_background("bg", 475, 8)
        bam = write_bam(self.tmp / "vcf_contig_mismatch.bam", [("chr1", 2000)], reads)
        vcf = write_vcf(self.tmp / "mismatched.vcf.gz", [("1", 2000)], ["1\t100\t.\tAA\tA\t60\tPASS\t."])
        settings = make_settings("ndigenome")
        with self.assertRaisesRegex(ValueError, "missing analyzed BAM contig"):
            run_caller(settings, {"bam": bam, "control_bam": None, "vcf": vcf, "blacklist": None}, self.tmp, chunks=1, sample="bad")

    def test_single_end_alignments_raise_error(self) -> None:
        """Single-end alignments → ValueError."""
        reads = [make_read(f"signal_{index}", 100) for index in range(8)]
        reads += forward_background("bg", 100, 12)
        # Remove paired bit from all reads
        for read in reads:
            read.flag &= ~1
        bam = write_bam(self.tmp / "single_end.bam", [("chr1", 2000)], reads)
        settings = make_settings("ndigenome")
        with self.assertRaisesRegex(ValueError, "requires paired-end alignments for nDigenome"):
            run_caller(settings, {"bam": bam, "control_bam": None, "vcf": None, "blacklist": None}, self.tmp, chunks=1, sample="bad")

    def test_focal_count_exactly_threshold(self) -> None:
        """Focal count exactly 5 (== ndigenome_min_count) → 1 row, SSB, PASS."""
        reads = [make_read(f"signal_{index}", 100) for index in range(5)]
        reads += forward_background("bg", 100, 20)
        bam = write_bam(self.tmp / "exact_threshold.bam", [("chr1", 2000)], reads)
        settings = make_settings("ndigenome")
        prefix, qc = run_caller(settings, {"bam": bam, "control_bam": None, "vcf": None, "blacklist": None}, self.tmp, chunks=1, sample="exact_threshold")
        rows = self.read_rows("exact_threshold")
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["signal_classification"], "SSB")
        self.assertEqual(rows[0]["filter_status"], "PASS")
        self.assertEqual(rows[0]["endpoint_count"], "5")

    def test_focal_count_below_threshold(self) -> None:
        """Focal count 4 → 0 rows."""
        reads = [make_read(f"signal_{index}", 100) for index in range(4)]
        reads += forward_background("bg", 100, 20)
        bam = write_bam(self.tmp / "below_threshold.bam", [("chr1", 2000)], reads)
        settings = make_settings("ndigenome")
        prefix, qc = run_caller(settings, {"bam": bam, "control_bam": None, "vcf": None, "blacklist": None}, self.tmp, chunks=1, sample="below_threshold")
        rows = self.read_rows("below_threshold")
        self.assertEqual(len(rows), 0)

    def test_opposite_endpoint_exactly_passing_fraction(self) -> None:
        """Opposite endpoint with exactly 5 reads and fraction exactly 0.20 → POSSIBLE_DSB."""
        reads = [make_read(f"forward_{index}", 100) for index in range(8)]
        reads += [reverse_read_ending_at(f"reverse_{index}", 101) for index in range(5)]
        reads += reverse_background("bg", 101, 20)
        bam = write_bam(self.tmp / "opposite_exact.bam", [("chr1", 2000)], reads)
        settings = make_settings("ndigenome")
        prefix, qc = run_caller(settings, {"bam": bam, "control_bam": None, "vcf": None, "blacklist": None}, self.tmp, chunks=1, sample="opposite_exact")
        rows = self.read_rows("opposite_exact")
        # Find the row at position 100 (the focal forward position)
        forward_rows = [r for r in rows if r["position_0based"] == "100" and r["strand"] == "+"]
        self.assertEqual(len(forward_rows), 1)
        self.assertEqual(forward_rows[0]["signal_classification"], "POSSIBLE_DSB")
        self.assertEqual(forward_rows[0]["opposite_count"], "5")
        self.assertAlmostEqual(float(forward_rows[0]["opposite_fraction"]), 0.2, places=5)

    def test_opposite_endpoint_ambiguous_threshold(self) -> None:
        """Opposite endpoint with exactly 2 reads (== ambiguous_min_count) and fraction below 0.05 → AMBIGUOUS."""
        reads = [make_read(f"forward_{index}", 100) for index in range(8)]
        reads += [reverse_read_ending_at(f"reverse_{index}", 101) for index in range(2)]
        reads += reverse_background("bg", 101, 50)
        bam = write_bam(self.tmp / "opposite_ambiguous.bam", [("chr1", 2000)], reads)
        settings = make_settings("ndigenome")
        prefix, qc = run_caller(settings, {"bam": bam, "control_bam": None, "vcf": None, "blacklist": None}, self.tmp, chunks=1, sample="opposite_ambiguous")
        rows = self.read_rows("opposite_ambiguous")
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["signal_classification"], "AMBIGUOUS")

    def test_serial_vs_chunked_equivalence(self) -> None:
        """Serial vs chunked: same all.tsv, tier TSVs, BED, and _mqc.tsv with chunks=1 and chunks=3."""
        reads = [make_read(f"fwd_{i}", 100 + i * 10) for i in range(8)]
        reads += [make_read(f"fwd2_{i}", 500 + i * 10) for i in range(6)]
        reads += forward_background("bg", 100, 10)
        reads += forward_background("bg2", 500, 8)
        bam = write_bam(self.tmp / "multisite.bam", [("chr1", 2000)], reads)

        settings = make_settings("ndigenome")

        # Serial run
        prefix_serial, qc_serial = run_caller(settings, {"bam": bam, "control_bam": None, "vcf": None, "blacklist": None}, self.tmp / "serial", chunks=1, sample="test")
        with open(f"{prefix_serial}.ndigenome.all.tsv") as f:
            serial_all = f.read()

        # Chunked run
        prefix_chunked, qc_chunked = run_caller(settings, {"bam": bam, "control_bam": None, "vcf": None, "blacklist": None}, self.tmp / "chunked", chunks=3, sample="test")
        with open(f"{prefix_chunked}.ndigenome.all.tsv") as f:
            chunked_all = f.read()

        self.assertEqual(serial_all, chunked_all)


if __name__ == "__main__":
    unittest.main()
