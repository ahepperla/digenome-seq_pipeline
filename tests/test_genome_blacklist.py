from __future__ import annotations

import gzip
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
from genome_blacklist import load_genome_blacklist  # noqa: E402


class GenomeBlacklistTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tempdir.name)
        self.lengths = {"chr1": 1000, "chr2": 500}

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_overlaps_are_merged_and_regions_are_subtracted(self) -> None:
        path = self.tmp / "blacklist.bed"
        path.write_text(
            "# comment\n"
            "chr1\t100\t200\n"
            "chr1\t150\t250\n"
            "chr1\t250\t275\n"
            "chr2\t10\t20\toptional-name\n"
        )
        blacklist = load_genome_blacklist(str(path), self.lengths)
        self.assertEqual(
            blacklist.intervals,
            {
                "chr1": ((100, 275),),
                "chr2": ((10, 20),),
            },
        )
        self.assertEqual(blacklist.interval_count, 2)
        self.assertEqual(blacklist.excluded_bases, 185)
        self.assertFalse(blacklist.contains("chr1", 99))
        self.assertTrue(blacklist.contains("chr1", 100))
        self.assertTrue(blacklist.contains("chr1", 274))
        self.assertFalse(blacklist.contains("chr1", 275))
        self.assertEqual(
            blacklist.subtract("chr1", 50, 300),
            [(50, 100), (275, 300)],
        )
        self.assertEqual(
            blacklist.overlap_bases("chr1", 50, 300),
            175,
        )

    def test_bgzip_style_bed_is_supported(self) -> None:
        path = self.tmp / "blacklist.bed.gz"
        with gzip.open(path, "wt") as handle:
            handle.write("chr1\t10\t20\n")
        blacklist = load_genome_blacklist(str(path), self.lengths)
        self.assertTrue(blacklist.contains("chr1", 15))

    def test_invalid_reference_coordinates_are_rejected(self) -> None:
        unknown = self.tmp / "unknown.bed"
        unknown.write_text("chrMissing\t0\t10\n")
        with self.assertRaisesRegex(ValueError, "absent from the BAM"):
            load_genome_blacklist(str(unknown), self.lengths)

        outside = self.tmp / "outside.bed"
        outside.write_text("chr1\t900\t1100\n")
        with self.assertRaisesRegex(ValueError, "outside chr1"):
            load_genome_blacklist(str(outside), self.lengths)

    def test_empty_blacklist_is_rejected_when_supplied(self) -> None:
        path = self.tmp / "empty.bed"
        path.write_text("# no intervals\n")
        with self.assertRaisesRegex(ValueError, "contains no BED intervals"):
            load_genome_blacklist(str(path), self.lengths)


if __name__ == "__main__":
    unittest.main()
