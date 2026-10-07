"""Tests for Digenome mode: paired endpoint calling and filtering."""

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
from helpers import (
    make_read,
    reverse_read_ending_at,
    forward_background,
    reverse_background,
    write_bam,
    write_vcf,
    run_caller,
    make_settings,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from cleavage.bam import SiteMetrics
from cleavage.digenome import (
    pair_score,
    caller_filter_reasons,
    PairCandidate,
    select_pairs,
)


@unittest.skipIf(pysam is None, "pysam is not installed")
class DigenomeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def read_rows(self, prefix: str) -> list[dict[str, str]]:
        """Read all rows from the digenome all.tsv output."""
        path = Path(f"{prefix}.digenome.all.tsv")
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_paired_endpoints_pass(self) -> None:
        """Forward 8 at 500, reverse 8 ending at 500 → 1 row, high_confidence."""
        reads = [make_read(f"fwd_{i}", 500) for i in range(8)]
        reads += [reverse_read_ending_at(f"rev_{i}", 500) for i in range(8)]
        bam = write_bam(self.tmp / "treated.bam", [("chr1", 2000)], reads)

        prefix, qc = run_caller(
            make_settings("digenome"),
            {"bam": bam, "control_bam": None, "vcf": None, "blacklist": None},
            self.tmp / "out",
            chunks=1,
        )
        rows = self.read_rows(prefix)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["forward_position_0based"], "500")
        self.assertEqual(rows[0]["forward_position_1based"], "501")
        self.assertEqual(rows[0]["reverse_position_0based"], "500")
        self.assertEqual(rows[0]["reverse_position_1based"], "501")
        self.assertEqual(rows[0]["tier"], "high_confidence")
        self.assertEqual(rows[0]["filter_status"], "PASS")
        self.assertEqual(qc["tiers"]["high_confidence"], 1)

        # Verify pair score is computed correctly: f_frac * r_frac * (f_count + r_count) / 4
        pair_score_val = float(rows[0]["digenome_pair_score"])
        f_count = int(rows[0]["forward_endpoint_count"])
        r_count = int(rows[0]["reverse_endpoint_count"])
        f_frac = float(rows[0]["forward_fraction"])
        r_frac = float(rows[0]["reverse_fraction"])
        expected = f_frac * r_frac * (f_count + r_count) / 4.0
        self.assertAlmostEqual(pair_score_val, expected)

        # RGEN score should also be filled
        self.assertNotEqual(rows[0]["rgen_digenome_score"], "")

    def test_pair_score_cutoff_filters(self) -> None:
        """Pair-score cutoff 100.0 → LOW_DIGENOME_PAIR_SCORE filter."""
        reads = [make_read(f"fwd_{i}", 520) for i in range(8)]
        reads += [make_read(f"rev_{i}", 471, reverse=True) for i in range(8)]
        bam = write_bam(self.tmp / "treated.bam", [("chr1", 2000)], reads)

        settings = make_settings("digenome", digenome_pair_score_cutoff=100.0)
        prefix, _ = run_caller(
            settings,
            {"bam": bam, "control_bam": None, "vcf": None, "blacklist": None},
            self.tmp / "out",
            chunks=1,
        )
        rows = self.read_rows(prefix)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["filter_reasons"], "LOW_DIGENOME_PAIR_SCORE")

    def test_one_endpoint_blacklisted(self) -> None:
        """Baseline without blacklist has reverse at 498; with blacklist, 0 rows."""
        reads = [make_read(f"fwd_{i}", 500) for i in range(8)]
        reads += [make_read(f"rev_{i}", 449, reverse=True) for i in range(8)]
        bam = write_bam(self.tmp / "treated.bam", [("chr1", 2000)], reads)

        # Baseline without blacklist
        prefix_baseline, _ = run_caller(
            make_settings("digenome"),
            {"bam": bam, "control_bam": None, "vcf": None, "blacklist": None},
            self.tmp / "baseline",
            chunks=1,
        )
        rows_baseline = self.read_rows(prefix_baseline)
        self.assertEqual(len(rows_baseline), 1)
        self.assertEqual(rows_baseline[0]["reverse_position_0based"], "498")

        # With blacklist covering reverse endpoint
        blacklist = self.tmp / "blacklist.bed"
        blacklist.write_text("chr1\t498\t499\n")
        prefix_masked, qc_masked = run_caller(
            make_settings("digenome"),
            {
                "bam": bam,
                "control_bam": None,
                "vcf": None,
                "blacklist": blacklist,
            },
            self.tmp / "masked",
            chunks=1,
        )
        rows_masked = self.read_rows(prefix_masked)
        self.assertEqual(len(rows_masked), 0)
        self.assertEqual(qc_masked["rows"], 0)

    def test_one_strand_only_no_pairs(self) -> None:
        """Only forward strand → 0 rows."""
        reads = [make_read(f"fwd_{i}", 550) for i in range(8)]
        bam = write_bam(self.tmp / "treated.bam", [("chr1", 2000)], reads)

        prefix, qc = run_caller(
            make_settings("digenome"),
            {"bam": bam, "control_bam": None, "vcf": None, "blacklist": None},
            self.tmp / "out",
            chunks=1,
        )
        rows = self.read_rows(prefix)
        self.assertEqual(len(rows), 0)
        self.assertEqual(qc["rows"], 0)

    def test_overhang_pairs_expected_coordinates(self) -> None:
        """Overhang 4, window 0 → forward 600, reverse 596."""
        reads = [make_read(f"fwd_{i}", 600) for i in range(8)]
        reads += [make_read(f"rev_{i}", 547, reverse=True) for i in range(8)]
        bam = write_bam(self.tmp / "treated.bam", [("chr1", 2000)], reads)

        settings = make_settings("digenome", digenome_overhang=4, digenome_pair_window=0)
        prefix, _ = run_caller(
            settings,
            {"bam": bam, "control_bam": None, "vcf": None, "blacklist": None},
            self.tmp / "out",
            chunks=1,
        )
        rows = self.read_rows(prefix)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["forward_position_0based"], "600")
        self.assertEqual(rows[0]["reverse_position_0based"], "596")

    def test_pair_score_formula(self) -> None:
        """pair_score(forward_count=8, forward_frac=0.5, reverse_count=12, reverse_frac=0.25) == 0.625."""
        forward = SiteMetrics(endpoint_count=8, endpoint_fraction=0.5)
        reverse = SiteMetrics(endpoint_count=12, endpoint_fraction=0.25)
        score = pair_score(forward, reverse)
        self.assertAlmostEqual(score, 0.625)

    def test_pairing_prefers_threshold_passing_pair(self) -> None:
        """When a threshold-passing pair exists, it is preferred over failing ones."""
        reads = [make_read(f"fwd_{i}", 1000) for i in range(8)]
        reads += forward_background("fbg", 1000, 12)
        reads += [make_read(f"rev_exact_{i}", 951, reverse=True) for i in range(8)]
        reads += [make_read(f"rev_nearby_{i}", 949, reverse=True) for i in range(8)]
        bam = write_bam(self.tmp / "treated.bam", [("chr1", 2000)], reads)

        settings = make_settings(
            "digenome",
            digenome_depth_cutoff=10,
            digenome_pair_score_cutoff=0.5,
        )
        prefix, _ = run_caller(
            settings,
            {"bam": bam, "control_bam": None, "vcf": None, "blacklist": None},
            self.tmp / "out",
            chunks=1,
        )
        rows = self.read_rows(prefix)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["reverse_position_0based"], "998")
        self.assertEqual(rows[0]["filter_status"], "PASS")

    def test_select_pairs_maximizes_nonconflicting_pairs(self) -> None:
        """Two unit tests: max pairs from (100,100), (100,101), (101,100) → (100,101), (101,100)."""
        metrics = SiteMetrics()
        candidates = [
            PairCandidate("chr1", 100, 100, metrics, metrics, 10.0, []),
            PairCandidate("chr1", 100, 101, metrics, metrics, 9.0, []),
            PairCandidate("chr1", 101, 100, metrics, metrics, 9.0, []),
        ]

        selected = select_pairs(candidates)

        self.assertEqual(
            [(c.forward_position, c.reverse_position) for c in selected],
            [(100, 101), (101, 100)],
        )

    def test_select_pairs_handles_exact_float_costs(self) -> None:
        """Exact float cost handling with binary fractions."""
        metrics = SiteMetrics()
        candidates = [
            PairCandidate("chr1", 1065, 1066, metrics, metrics, 0.2599142252422858, []),
            PairCandidate("chr1", 1067, 1066, metrics, metrics, 0.49406008082185116, []),
            PairCandidate("chr1", 1067, 1068, metrics, metrics, 2.863272311212815, []),
            PairCandidate("chr1", 1067, 1069, metrics, metrics, 2.863272311212815, []),
        ]

        selected = select_pairs(candidates)

        self.assertEqual(
            [(c.forward_position, c.reverse_position) for c in selected],
            [(1065, 1066), (1067, 1068)],
        )

    def test_select_pairs_order_independent(self) -> None:
        """Result does not depend on input order."""
        metrics = SiteMetrics()
        candidates = [
            PairCandidate("chr1", 1065, 1066, metrics, metrics, 0.2599142252422858, []),
            PairCandidate("chr1", 1067, 1066, metrics, metrics, 0.49406008082185116, []),
            PairCandidate("chr1", 1067, 1068, metrics, metrics, 2.863272311212815, []),
            PairCandidate("chr1", 1067, 1069, metrics, metrics, 2.863272311212815, []),
        ]

        # Reverse the order
        selected1 = select_pairs(candidates)
        selected2 = select_pairs(list(reversed(candidates)))

        result1 = [(c.forward_position, c.reverse_position) for c in selected1]
        result2 = [(c.forward_position, c.reverse_position) for c in selected2]
        self.assertEqual(result1, result2)

    def test_control_evidence_filters_pair(self) -> None:
        """Matched control with all three failure reasons."""
        treated_reads = [make_read(f"t_fwd_{i}", 700) for i in range(10)]
        treated_reads += [make_read(f"t_rev_{i}", 651, reverse=True) for i in range(10)]
        control_reads = [make_read(f"c_fwd_{i}", 700) for i in range(10)]
        control_reads += [make_read(f"c_rev_{i}", 651, reverse=True) for i in range(10)]

        treated_bam = write_bam(self.tmp / "treated.bam", [("chr1", 2000)], treated_reads)
        control_bam = write_bam(self.tmp / "control.bam", [("chr1", 2000)], control_reads)

        prefix, _ = run_caller(
            make_settings("digenome"),
            {
                "bam": treated_bam,
                "control_bam": control_bam,
                "vcf": None,
                "blacklist": None,
            },
            self.tmp / "out",
            chunks=1,
            sample="Sample",
            control_sample="Control",
        )
        rows = self.read_rows(prefix)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["control_status"], "MATCHED_CONTROL")
        self.assertEqual(
            rows[0]["digenome_pair_score"],
            rows[0]["control_digenome_pair_score"],
        )
        self.assertEqual(
            rows[0]["rgen_digenome_score"],
            rows[0]["control_rgen_digenome_score"],
        )
        self.assertEqual(rows[0]["tier"], "artifact")
        self.assertEqual(rows[0]["filter_reasons"], "HIGH_CONTROL_FRACTION;LOW_CONTROL_FOLD;CONTROL_Q_FAIL")

    def test_indel_and_vcf_artifacts_are_shared(self) -> None:
        """Both endpoints share indel and known-indel artifact reasons."""
        # Forward reads with insertion at position 805
        forward_cigar = [(0, 5), (1, 1), (0, 44)]
        reads = [make_read(f"fwd_indel_{i}", 800, cigar=forward_cigar) for i in range(8)]
        reads += [make_read(f"rev_{i}", 751, reverse=True) for i in range(8)]
        bam = write_bam(self.tmp / "treated.bam", [("chr1", 2000)], reads)

        # VCF with indel at position 805 (1-based in VCF, so POS=805, REF=AA, ALT=A is a deletion)
        vcf = write_vcf(
            self.tmp / "variants.vcf.gz",
            [("chr1", 2000)],
            ["chr1\t805\t.\tAA\tA\t60\tPASS\t."],
        )

        prefix, _ = run_caller(
            make_settings("digenome"),
            {"bam": bam, "control_bam": None, "vcf": vcf, "blacklist": None},
            self.tmp / "out",
            chunks=1,
        )
        rows = self.read_rows(prefix)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["tier"], "artifact")
        self.assertEqual(rows[0]["filter_reasons"], "NEARBY_INDEL;KNOWN_INDEL")

    def test_known_indel_overlap_union(self) -> None:
        """Known indel within artifact window of both endpoints appears exactly once."""
        # Forward endpoint 800 and reverse endpoint 798: both artifact windows
        # (10 bp) include the VCF indel at 0-based 804.
        reads = [make_read(f"fwd_{i}", 800) for i in range(8)]
        reads += [reverse_read_ending_at(f"rev_{i}", 798) for i in range(8)]
        reads += forward_background("fbg", 800, 12)
        reads += reverse_background("rbg", 798, 12)
        bam = write_bam(self.tmp / "treated.bam", [("chr1", 2000)], reads)

        # VCF indel at position 805 (within artifact window of both 800 and 751)
        vcf = write_vcf(
            self.tmp / "variants.vcf.gz",
            [("chr1", 2000)],
            ["chr1\t805\t.\tAA\tA\t60\tPASS\t."],
        )

        prefix, _ = run_caller(
            make_settings("digenome"),
            {"bam": bam, "control_bam": None, "vcf": vcf, "blacklist": None},
            self.tmp / "out",
            chunks=1,
        )
        rows = self.read_rows(prefix)

        self.assertEqual(len(rows), 1)
        # known_indel_overlap should contain the single overlap (chr1:804 or similar)
        self.assertEqual(rows[0]["reverse_position_0based"], "798")
        self.assertEqual(rows[0]["known_indel_overlap"], "chr1:805:AA>A")

    def test_strict_cutoff_low_forward_count(self) -> None:
        """Forward count exactly equal to cutoff → LOW_FORWARD_COUNT filter."""
        forward = SiteMetrics(endpoint_count=4, endpoint_fraction=0.5, strand_depth=8)
        reverse = SiteMetrics(endpoint_count=8, endpoint_fraction=0.5, strand_depth=16)
        score = pair_score(forward, reverse)

        settings = make_settings("digenome", digenome_forward_cutoff=4, digenome_pair_score_cutoff=0.0)
        reasons = caller_filter_reasons(forward, reverse, score, settings)

        self.assertEqual(reasons, ["LOW_FORWARD_COUNT"])

    def test_strict_cutoff_forward_count_just_above(self) -> None:
        """Forward count one above cutoff → no LOW_FORWARD_COUNT."""
        forward = SiteMetrics(endpoint_count=5, endpoint_fraction=0.5, strand_depth=10)
        reverse = SiteMetrics(endpoint_count=8, endpoint_fraction=0.5, strand_depth=16)
        score = pair_score(forward, reverse)

        settings = make_settings("digenome", digenome_forward_cutoff=4, digenome_pair_score_cutoff=0.0)
        reasons = caller_filter_reasons(forward, reverse, score, settings)

        self.assertEqual(reasons, [])

    def test_strict_cutoff_low_forward_depth(self) -> None:
        """Strand depth exactly equal to cutoff → LOW_FORWARD_DEPTH."""
        forward = SiteMetrics(endpoint_count=5, endpoint_fraction=0.5, strand_depth=4)
        reverse = SiteMetrics(endpoint_count=8, endpoint_fraction=0.5, strand_depth=16)
        score = pair_score(forward, reverse)

        settings = make_settings("digenome", digenome_depth_cutoff=4, digenome_pair_score_cutoff=0.0)
        reasons = caller_filter_reasons(forward, reverse, score, settings)

        self.assertEqual(reasons, ["LOW_FORWARD_DEPTH"])

    def test_duplicate_endpoint_records_error(self) -> None:
        """Duplicate endpoint records on same contig raise ValueError."""
        from cleavage.digenome import pair_rows

        records = [
            {
                "contig": "chr1",
                "position_0based": 100,
                "strand": "+",
                "metrics": {"endpoint_count": 8, "endpoint_fraction": 0.5, "strand_depth": 16},
                "control_metrics": None,
                "known_indels": [],
                "rgen_score": 0.0,
                "control_rgen_score": None,
            },
            {
                "contig": "chr1",
                "position_0based": 100,
                "strand": "+",
                "metrics": {"endpoint_count": 5, "endpoint_fraction": 0.5, "strand_depth": 10},
                "control_metrics": None,
                "known_indels": [],
                "rgen_score": 0.0,
                "control_rgen_score": None,
            },
        ]

        with self.assertRaisesRegex(ValueError, "Duplicate Digenome endpoint records"):
            pair_rows(records, make_settings("digenome"), "Sample", "")


if __name__ == "__main__":
    unittest.main()
