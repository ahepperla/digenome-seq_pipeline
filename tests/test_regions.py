"""Tests for blacklist and chunk planning: load_blacklist, plan_chunks, check_plan, etc."""

from __future__ import annotations

import gzip
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from helpers import write_bam, make_read  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from cleavage.regions import (  # noqa: E402
    Blacklist,
    load_blacklist,
    plan_chunks,
    check_plan,
    write_plan,
    plan_to_json,
    plan_from_json,
    MappedContig,
    OwnedInterval,
)


class LoadBlacklistTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tempdir.name)
        self.contig_lengths = {"chr1": 1000, "chr2": 500}

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_overlapping_regions_merge_and_contain_works(self) -> None:
        path = self.tmp / "blacklist.bed"
        path.write_text(
            "# comment\n"
            "chr1\t100\t200\n"
            "chr1\t150\t250\n"
            "chr1\t250\t275\n"
            "chr2\t10\t20\toptional-name\n"
        )
        blacklist = load_blacklist(path, self.contig_lengths)
        self.assertEqual(
            blacklist.intervals,
            {"chr1": [(100, 275)], "chr2": [(10, 20)]},
        )
        self.assertEqual(blacklist.interval_count, 2)
        self.assertEqual(blacklist.excluded_bases, 185)
        self.assertFalse(blacklist.contains("chr1", 99))
        self.assertTrue(blacklist.contains("chr1", 100))
        self.assertTrue(blacklist.contains("chr1", 274))
        self.assertFalse(blacklist.contains("chr1", 275))

    def test_subtract_excludes_blacklisted_regions(self) -> None:
        path = self.tmp / "blacklist.bed"
        path.write_text(
            "chr1\t100\t200\n"
            "chr1\t150\t250\n"
            "chr1\t250\t275\n"
        )
        blacklist = load_blacklist(path, self.contig_lengths)
        result = blacklist.subtract("chr1", 50, 300)
        self.assertEqual(result, [(50, 100), (275, 300)])

    def test_provenance_returns_file_info_and_counts(self) -> None:
        path = self.tmp / "blacklist.bed"
        path.write_text(
            "chr1\t100\t200\n"
            "chr1\t150\t250\n"
            "chr1\t250\t275\n"
            "chr2\t10\t20\n"
        )
        blacklist = load_blacklist(path, self.contig_lengths)
        prov = blacklist.provenance()
        self.assertEqual(prov["file"], "blacklist.bed")
        self.assertIn("sha256", prov)
        self.assertEqual(len(prov["sha256"]), 64)
        self.assertEqual(prov["intervals"], 2)
        self.assertEqual(prov["excluded_bases"], 185)

    def test_bgzip_bed_is_supported(self) -> None:
        path = self.tmp / "blacklist.bed.gz"
        with gzip.open(path, "wt") as handle:
            handle.write("chr1\t10\t20\n")
        blacklist = load_blacklist(path, self.contig_lengths)
        self.assertTrue(blacklist.contains("chr1", 15))

    def test_unknown_contig_raises_error(self) -> None:
        path = self.tmp / "unknown.bed"
        path.write_text("chrMissing\t0\t10\n")
        with self.assertRaisesRegex(ValueError, "absent from the BAM"):
            load_blacklist(path, self.contig_lengths)

    def test_coordinates_outside_contig_raise_error(self) -> None:
        path = self.tmp / "outside.bed"
        path.write_text("chr1\t900\t1100\n")
        with self.assertRaisesRegex(ValueError, "outside chr1"):
            load_blacklist(path, self.contig_lengths)

    def test_empty_blacklist_raises_error(self) -> None:
        path = self.tmp / "empty.bed"
        path.write_text("# no intervals\n")
        with self.assertRaisesRegex(ValueError, "contains no BED intervals"):
            load_blacklist(path, self.contig_lengths)

    def test_fewer_than_three_columns_raises_error(self) -> None:
        path = self.tmp / "short.bed"
        path.write_text("chr1\t100\n")
        with self.assertRaisesRegex(ValueError, "fewer than three BED columns"):
            load_blacklist(path, self.contig_lengths)

    def test_non_integer_coordinates_raise_error(self) -> None:
        path = self.tmp / "invalid.bed"
        path.write_text("chr1\tabc\t200\n")
        with self.assertRaisesRegex(ValueError, "invalid coordinates"):
            load_blacklist(path, self.contig_lengths)


class PlanChunksTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_planner_balances_mapped_contigs_into_chunks(self) -> None:
        contigs = [
            MappedContig("chr1", 2000, 8),
            MappedContig("chr2", 2000, 5),
            MappedContig("chr3", 2000, 3),
        ]
        plan = plan_chunks(contigs, 2, None)
        self.assertEqual(len(plan), 2)
        self.assertEqual(len(plan[0]), 1)
        self.assertEqual(plan[0][0].contig, "chr1")
        self.assertEqual(plan[0][0].start, 0)
        self.assertEqual(plan[0][0].end, 2000)
        self.assertEqual(plan[0][0].callable_bases, 2000)
        self.assertEqual(len(plan[1]), 2)
        self.assertTrue(any(i.contig == "chr2" for i in plan[1]))
        self.assertTrue(any(i.contig == "chr3" for i in plan[1]))

    def test_planner_splits_large_contigs_across_many_chunks(self) -> None:
        contigs = [
            MappedContig("chr1", 2000, 8),
            MappedContig("chr2", 2000, 5),
            MappedContig("chr3", 2000, 3),
        ]
        plan = plan_chunks(contigs, 20, None)
        self.assertEqual(len(plan), 20)
        non_empty = [slot for slot in plan if slot]
        self.assertGreater(len(non_empty), 1)
        self.assertTrue(any(i.start > 0 or i.end < 2000 for slot in plan for i in slot))
        check_plan(plan, {"chr1": 2000, "chr2": 2000, "chr3": 2000})

    def test_planner_with_blacklist_reduces_callable_bases(self) -> None:
        path = self.tmp / "blacklist.bed"
        path.write_text("chr1\t0\t1800\n")
        blacklist = load_blacklist(path, {"chr1": 2000, "chr2": 2000})
        contigs = [
            MappedContig("chr1", 2000, 100),
            MappedContig("chr2", 2000, 100),
        ]
        plan = plan_chunks(contigs, 2, blacklist)
        self.assertEqual(len(plan), 2)
        callable_sum = sum(
            i.callable_bases for slot in plan for i in slot
        )
        self.assertEqual(callable_sum, 2200)

    def test_planner_handles_fully_blacklisted_mapped_contig(self) -> None:
        path = self.tmp / "fully_blacklist.bed"
        path.write_text("chr1\t0\t2000\n")
        blacklist = load_blacklist(path, {"chr1": 2000})
        contigs = [MappedContig("chr1", 2000, 10)]
        plan = plan_chunks(contigs, 8, blacklist)
        self.assertEqual(len(plan), 8)
        found_interval = None
        for slot in plan:
            if slot:
                self.assertEqual(len(slot), 1)
                found_interval = slot[0]
        self.assertIsNotNone(found_interval)
        self.assertEqual(found_interval.start, 0)
        self.assertEqual(found_interval.end, 2000)
        self.assertEqual(found_interval.callable_bases, 0)

    def test_planner_produces_deterministic_output(self) -> None:
        contigs = [
            MappedContig("chr1", 2000, 8),
            MappedContig("chr2", 2000, 5),
        ]
        plan1 = plan_chunks(contigs, 3, None)
        plan2 = plan_chunks(contigs, 3, None)
        self.assertEqual(plan1, plan2)


class CheckPlanTests(unittest.TestCase):
    def test_check_plan_accepts_valid_plan(self) -> None:
        plan = [
            [OwnedInterval("chr1", 0, 1000, 1000)],
            [OwnedInterval("chr1", 1000, 2000, 1000)],
        ]
        contig_lengths = {"chr1": 2000}
        check_plan(plan, contig_lengths)

    def test_check_plan_detects_gaps(self) -> None:
        plan = [
            [OwnedInterval("chr1", 0, 500, 500)],
            [OwnedInterval("chr1", 1500, 2000, 500)],
        ]
        contig_lengths = {"chr1": 2000}
        with self.assertRaisesRegex(ValueError, "gap on chr1"):
            check_plan(plan, contig_lengths)

    def test_check_plan_detects_overlaps(self) -> None:
        plan = [
            [OwnedInterval("chr1", 0, 1200, 1200)],
            [OwnedInterval("chr1", 1000, 2000, 1000)],
        ]
        contig_lengths = {"chr1": 2000}
        with self.assertRaisesRegex(ValueError, "overlaps on chr1"):
            check_plan(plan, contig_lengths)

    def test_check_plan_detects_unexpected_contig(self) -> None:
        plan = [
            [OwnedInterval("chrX", 0, 1000, 1000)],
        ]
        contig_lengths = {"chr1": 2000}
        with self.assertRaisesRegex(ValueError, "unexpected contig"):
            check_plan(plan, contig_lengths)


class WritePlanTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_write_plan_produces_correct_tsv_format(self) -> None:
        path = self.tmp / "plan.tsv"
        plan = [
            [
                OwnedInterval("chr1", 0, 1000, 900),
                OwnedInterval("chr2", 0, 500, 500),
            ],
            [
                OwnedInterval("chr1", 1000, 2000, 1000),
            ],
        ]
        write_plan(plan, str(path))
        content = path.read_text()
        lines = content.strip().split("\n")
        self.assertEqual(lines[0], "chunk\tcontig\tstart\tend\tcallable_bases\texcluded_bases")
        self.assertEqual(lines[1], "0\tchr1\t0\t1000\t900\t100")
        self.assertEqual(lines[2], "0\tchr2\t0\t500\t500\t0")
        self.assertEqual(lines[3], "1\tchr1\t1000\t2000\t1000\t0")

    def test_plan_roundtrip_through_json(self) -> None:
        plan = [
            [
                OwnedInterval("chr1", 0, 500, 500),
                OwnedInterval("chr2", 100, 600, 500),
            ],
            [
                OwnedInterval("chr1", 500, 2000, 1500),
            ],
        ]
        json_data = plan_to_json(plan)
        restored = plan_from_json(json_data)
        self.assertEqual(restored, plan)


if __name__ == "__main__":
    unittest.main()
