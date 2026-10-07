"""Tests for BAM reading: endpoints, site metrics, RGEN scoring."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from helpers import make_read, reverse_read_ending_at, forward_background, write_bam  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from cleavage.bam import (  # noqa: E402
    endpoint_position,
    five_prime_clip_length,
    count_endpoints,
    measure_site,
    SiteMetrics,
    combine,
    rgen_counts,
    rgen_score,
    RgenCounts,
)


class EndpointCoordinatesTests(unittest.TestCase):
    def test_forward_read_endpoint_is_start_position(self) -> None:
        read = make_read("forward", 100)
        self.assertEqual(endpoint_position(read), 100)

    def test_reverse_read_endpoint_is_end_minus_one(self) -> None:
        read = reverse_read_ending_at("reverse", 100)
        self.assertEqual(endpoint_position(read), 100)

    def test_complex_cigar_forward_endpoint(self) -> None:
        read = make_read("complex_fwd", 100, cigar=[(4, 5), (0, 45)])
        self.assertEqual(endpoint_position(read), 100)
        self.assertEqual(five_prime_clip_length(read), 5)

    def test_complex_cigar_reverse_endpoint(self) -> None:
        read = make_read("complex_rev", 100, reverse=True, cigar=[(5, 2), (0, 20), (2, 3), (0, 20), (4, 5)])
        self.assertEqual(endpoint_position(read), 142)
        self.assertEqual(five_prime_clip_length(read), 5)


class CountEndpointsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_count_endpoints_counts_forward_and_reverse_separately(self) -> None:
        reads = [
            make_read(f"forward_{i}", 100, mapq=1) for i in range(3)
        ]
        reads += [
            reverse_read_ending_at(f"reverse_{i}", 300, mapq=1) for i in range(3)
        ]
        bam_path = write_bam(self.tmp / "count_test.bam", [("chr1", 2000)], reads)
        with __import__("pysam").AlignmentFile(str(bam_path), "rb") as bam:
            counts = count_endpoints(bam, "chr1", 0, 2000, min_mapq=1, min_count=3)
        self.assertEqual(counts, {(100, "+"): 3, (300, "-"): 3})

    def test_count_endpoints_filters_below_min_count(self) -> None:
        reads = [
            make_read(f"forward_{i}", 100, mapq=1) for i in range(2)
        ]
        reads += [
            reverse_read_ending_at(f"reverse_{i}", 300, mapq=1) for i in range(3)
        ]
        bam_path = write_bam(self.tmp / "min_count_test.bam", [("chr1", 2000)], reads)
        with __import__("pysam").AlignmentFile(str(bam_path), "rb") as bam:
            counts = count_endpoints(bam, "chr1", 0, 2000, min_mapq=1, min_count=3)
        self.assertEqual(counts, {(300, "-"): 3})

    def test_count_endpoints_filters_by_flags(self) -> None:
        reads = [make_read("primary", 100, mapq=1)]
        reads += [make_read(f"duplicate_{i}", 700, mapq=1, flag=1024) for i in range(3)]
        reads += [make_read(f"secondary_{i}", 700, mapq=1, flag=256) for i in range(3)]
        reads += [make_read(f"supplementary_{i}", 700, mapq=1, flag=2048) for i in range(3)]
        reads += [make_read(f"qcfail_{i}", 700, mapq=1, flag=512) for i in range(3)]
        bam_path = write_bam(self.tmp / "flags_test.bam", [("chr1", 2000)], reads)
        with __import__("pysam").AlignmentFile(str(bam_path), "rb") as bam:
            counts = count_endpoints(bam, "chr1", 0, 2000, min_mapq=1, min_count=1)
        self.assertEqual(counts, {(100, "+"): 1})

    def test_count_endpoints_filters_by_mapq(self) -> None:
        reads = [make_read("high_mapq", 100, mapq=30)]
        reads += [make_read(f"low_mapq_{i}", 500, mapq=0) for i in range(3)]
        bam_path = write_bam(self.tmp / "mapq_test.bam", [("chr1", 2000)], reads)
        with __import__("pysam").AlignmentFile(str(bam_path), "rb") as bam:
            counts = count_endpoints(bam, "chr1", 0, 2000, min_mapq=1, min_count=1)
        self.assertEqual(counts, {(100, "+"): 1})

    def test_count_endpoints_respects_range(self) -> None:
        reads = [
            make_read(f"at_100_{i}", 100, mapq=1) for i in range(3)
        ]
        reads += [
            reverse_read_ending_at(f"at_300_{i}", 300, mapq=1) for i in range(3)
        ]
        bam_path = write_bam(self.tmp / "range_test.bam", [("chr1", 2000)], reads)
        with __import__("pysam").AlignmentFile(str(bam_path), "rb") as bam:
            counts = count_endpoints(bam, "chr1", 200, 400, min_mapq=1, min_count=1)
        self.assertEqual(counts, {(300, "-"): 3})
        self.assertNotIn((100, "+"), counts)


class MeasureSiteTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_deletion_ending_at_window_start_is_not_nearby(self) -> None:
        read = make_read("deletion_read", 280, cigar=[(0, 5), (2, 5), (0, 30)])
        bam_path = write_bam(self.tmp / "deletion_test.bam", [("chr1", 2000)], [read])
        with __import__("pysam").AlignmentFile(str(bam_path), "rb") as bam:
            metrics = measure_site(bam, "chr1", 300, "+", window=10, min_mapq=1)
        self.assertEqual(metrics.indel_read_count, 0)
        self.assertIsNone(metrics.indel_position)

    def test_measure_site_counts_endpoints_and_depth(self) -> None:
        reads = [
            make_read(f"endpoint_{i}", 100, mapq=60) for i in range(6)
        ]
        reads += forward_background("bg", 100, 4, mapq=60)
        bam_path = write_bam(self.tmp / "measure_test.bam", [("chr1", 2000)], reads)
        with __import__("pysam").AlignmentFile(str(bam_path), "rb") as bam:
            metrics = measure_site(bam, "chr1", 100, "+", window=10, min_mapq=1)
        self.assertEqual(metrics.endpoint_count, 6)
        self.assertEqual(metrics.strand_depth, 10)
        self.assertAlmostEqual(metrics.endpoint_fraction, 0.6, places=5)

    def test_measure_site_counts_secondary_endpoints_separately(self) -> None:
        reads = [
            make_read(f"primary_{i}", 100, mapq=60) for i in range(6)
        ]
        reads += [
            make_read(f"secondary_{i}", 100, mapq=60, flag=256) for i in range(2)
        ]
        reads += forward_background("bg", 100, 4, mapq=60)
        bam_path = write_bam(self.tmp / "secondary_test.bam", [("chr1", 2000)], reads)
        with __import__("pysam").AlignmentFile(str(bam_path), "rb") as bam:
            metrics = measure_site(bam, "chr1", 100, "+", window=10, min_mapq=1)
        self.assertEqual(metrics.endpoint_count, 6)
        self.assertEqual(metrics.secondary_endpoint_count, 2)


class CombineMetricsTests(unittest.TestCase):
    def test_combine_sums_endpoint_counts(self) -> None:
        fwd = SiteMetrics(endpoint_count=6, strand_depth=10)
        rev = SiteMetrics(endpoint_count=4, strand_depth=8)
        combined = combine(fwd, rev)
        self.assertEqual(combined.endpoint_count, 10)
        self.assertEqual(combined.strand_depth, 18)

    def test_combine_tie_rule_prefers_forward_indel(self) -> None:
        fwd = SiteMetrics(
            endpoint_count=1,
            strand_depth=2,
            indel_fraction=0.5,
            indel_read_count=1,
            indel_position=100,
            indel_type="DEL",
            indel_length=3,
        )
        rev = SiteMetrics(
            endpoint_count=1,
            strand_depth=2,
            indel_fraction=0.5,
            indel_read_count=1,
            indel_position=200,
            indel_type="INS",
            indel_length=2,
        )
        combined = combine(fwd, rev)
        self.assertEqual(combined.indel_position, 100)
        self.assertEqual(combined.indel_type, "DEL")
        self.assertEqual(combined.indel_length, 3)

    def test_combine_swapping_arguments_changes_tie_result(self) -> None:
        fwd = SiteMetrics(
            endpoint_count=1,
            strand_depth=2,
            indel_fraction=0.5,
            indel_read_count=1,
            indel_position=100,
            indel_type="DEL",
            indel_length=3,
        )
        rev = SiteMetrics(
            endpoint_count=1,
            strand_depth=2,
            indel_fraction=0.5,
            indel_read_count=1,
            indel_position=200,
            indel_type="INS",
            indel_length=2,
        )
        combined1 = combine(fwd, rev)
        combined2 = combine(rev, fwd)
        self.assertEqual(combined1.indel_position, 100)
        self.assertEqual(combined2.indel_position, 200)


class RgenCountsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_rgen_counts_includes_supplementary_not_secondary(self) -> None:
        reads = [
            make_read("primary", 900, mapq=1),
            make_read("supplementary", 900, mapq=1, flag=2048),
            make_read("secondary", 900, mapq=1, flag=256),
            make_read("duplicate", 900, mapq=1, flag=1024),
        ]
        bam_path = write_bam(self.tmp / "rgen_test.bam", [("chr1", 2000)], reads)
        with __import__("pysam").AlignmentFile(str(bam_path), "rb") as bam:
            counts = rgen_counts(bam, "chr1", 900, min_mapq=1)
        self.assertEqual(counts.forward_endpoints, 2)
        self.assertEqual(counts.reverse_endpoints, 0)
        self.assertEqual(counts.depth, 2)


class RgenScoreTests(unittest.TestCase):
    def test_rgen_score_with_float32_math(self) -> None:
        table = {
            95: RgenCounts(1, 1, 20),
            96: RgenCounts(1, 1, 20),
            97: RgenCounts(1, 1, 20),
            98: RgenCounts(1, 1, 20),
            99: RgenCounts(1, 8, 16),
            100: RgenCounts(10, 1, 20),
            101: RgenCounts(1, 1, 20),
            102: RgenCounts(1, 1, 20),
            103: RgenCounts(1, 1, 20),
            104: RgenCounts(1, 1, 20),
            105: RgenCounts(1, 1, 20),
        }
        counts_at = lambda position: table.get(position)
        score = rgen_score(counts_at, 100, 0)
        self.assertAlmostEqual(score, 6.300000190734863, places=8)

    def test_rgen_score_returns_zero_when_blacklist_blocks_anchor(self) -> None:
        table = {
            position: RgenCounts(1, 1, 20)
            for position in range(95, 106)
        }
        table[99] = RgenCounts(1, 8, 16)
        table[100] = RgenCounts(10, 1, 20)
        counts_at = lambda position: table.get(position) if position != 99 else None
        score = rgen_score(counts_at, 100, 0)
        self.assertEqual(score, 0.0)


if __name__ == "__main__":
    unittest.main()
