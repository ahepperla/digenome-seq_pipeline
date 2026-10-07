"""Tests for output filtering, formatting, and file writing."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from helpers import make_settings  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from cleavage.output import (  # noqa: E402
    DIGENOME_COLUMNS,
    NDIGENOME_COLUMNS,
    apply_filters,
    bed_line,
    format_value,
)
from cleavage.stats import MATCHED_CONTROL  # noqa: E402


class OutputFilterTests(unittest.TestCase):
    def test_apply_filters_reports_each_failed_control_criterion(self) -> None:
        """apply_filters reports each failed control criterion on its own."""
        settings = make_settings("ndigenome")

        cases = {
            "HIGH_CONTROL_FRACTION": {
                "control_fraction": 0.10,
                "control_fold_enrichment": 10.0,
                "control_fisher_q": 0.01,
            },
            "LOW_CONTROL_FOLD": {
                "control_fraction": 0.01,
                "control_fold_enrichment": 2.0,
                "control_fisher_q": 0.01,
            },
            "CONTROL_Q_FAIL": {
                "control_fraction": 0.01,
                "control_fold_enrichment": 10.0,
                "control_fisher_q": 0.10,
            },
        }

        for expected_reason, control_values in cases.items():
            with self.subTest(expected_reason=expected_reason):
                row = {
                    "softclip_fraction": 0.0,
                    "indel_fraction": 0.0,
                    "known_indel_overlap": "",
                    "support_mean_mapq": 60.0,
                    "control_status": MATCHED_CONTROL,
                    "signal_classification": "SSB",
                    **control_values,
                }
                apply_filters(row, settings)
                self.assertEqual(row["filter_reasons"], expected_reason)
                self.assertEqual(row["filter_status"], "FILTERED")
                self.assertEqual(row["tier"], "artifact")

    def test_softclip_fraction_boundary(self) -> None:
        """softclip_fraction == 0.20 → HIGH_5P_SOFTCLIP (fails at >=)."""
        settings = make_settings("ndigenome")
        row = {
            "softclip_fraction": 0.20,
            "indel_fraction": 0.0,
            "known_indel_overlap": "",
            "support_mean_mapq": 60.0,
            "control_status": "UNCONTROLLED",
            "signal_classification": "SSB",
        }
        apply_filters(row, settings)
        self.assertEqual(row["filter_reasons"], "HIGH_5P_SOFTCLIP")
        self.assertEqual(row["filter_status"], "FILTERED")
        self.assertEqual(row["tier"], "artifact")

    def test_indel_fraction_boundary(self) -> None:
        """indel_fraction == 0.20 → NEARBY_INDEL."""
        settings = make_settings("ndigenome")
        row = {
            "softclip_fraction": 0.0,
            "indel_fraction": 0.20,
            "known_indel_overlap": "",
            "support_mean_mapq": 60.0,
            "control_status": "UNCONTROLLED",
            "signal_classification": "SSB",
        }
        apply_filters(row, settings)
        self.assertEqual(row["filter_reasons"], "NEARBY_INDEL")
        self.assertEqual(row["filter_status"], "FILTERED")

    def test_support_mean_mapq_boundary(self) -> None:
        """support_mean_mapq == 10.0 → no LOW_SUPPORT_MAPQ (fails only below)."""
        settings = make_settings("ndigenome")
        row = {
            "softclip_fraction": 0.0,
            "indel_fraction": 0.0,
            "known_indel_overlap": "",
            "support_mean_mapq": 10.0,
            "control_status": "UNCONTROLLED",
            "signal_classification": "SSB",
        }
        apply_filters(row, settings)
        self.assertNotIn("LOW_SUPPORT_MAPQ", row["filter_reasons"])
        self.assertEqual(row["filter_status"], "PASS")

    def test_control_boundaries_pass(self) -> None:
        """control_fraction == 0.05, fold == 3.0, q == 0.05 → PASS."""
        settings = make_settings("ndigenome")
        row = {
            "softclip_fraction": 0.0,
            "indel_fraction": 0.0,
            "known_indel_overlap": "",
            "support_mean_mapq": 60.0,
            "control_status": MATCHED_CONTROL,
            "signal_classification": "SSB",
            "control_fraction": 0.05,
            "control_fold_enrichment": 3.0,
            "control_fisher_q": 0.05,
        }
        apply_filters(row, settings)
        self.assertEqual(row["filter_status"], "PASS")
        self.assertEqual(row["tier"], "high_confidence")

    def test_digenome_uses_combined_control_fraction(self) -> None:
        """A Digenome row uses control_combined_fraction for HIGH_CONTROL_FRACTION."""
        settings = make_settings("digenome")
        row = {
            "softclip_fraction": 0.0,
            "indel_fraction": 0.0,
            "known_indel_overlap": "",
            "support_mean_mapq": 60.0,
            "control_status": MATCHED_CONTROL,
            "control_combined_fraction": 0.10,
            "control_fold_enrichment": 10.0,
            "control_fisher_q": 0.01,
        }
        apply_filters(row, settings)
        self.assertEqual(row["filter_reasons"], "HIGH_CONTROL_FRACTION")

    def test_digenome_caller_filter_reasons_come_first(self) -> None:
        """Digenome: caller_filter_reasons come first, then artifact reasons."""
        settings = make_settings("digenome")
        row = {
            "softclip_fraction": 0.30,
            "indel_fraction": 0.0,
            "known_indel_overlap": "",
            "support_mean_mapq": 60.0,
            "control_status": "UNCONTROLLED",
            "caller_filter_reasons": ["LOW_DIGENOME_PAIR_SCORE"],
        }
        apply_filters(row, settings)
        self.assertEqual(row["filter_reasons"], "LOW_DIGENOME_PAIR_SCORE;HIGH_5P_SOFTCLIP")
        self.assertEqual(row["tier"], "artifact")

    def test_ndigenome_non_ssb_appends_signal_class(self) -> None:
        """nDigenome non-SSB rows get their signal class appended as a reason."""
        settings = make_settings("ndigenome")
        row = {
            "softclip_fraction": 0.0,
            "indel_fraction": 0.0,
            "known_indel_overlap": "",
            "support_mean_mapq": 60.0,
            "control_status": "UNCONTROLLED",
            "signal_classification": "POSSIBLE_DSB",
        }
        apply_filters(row, settings)
        self.assertEqual(row["filter_reasons"], "POSSIBLE_DSB")
        self.assertEqual(row["tier"], "manual_review")

    def test_ndigenome_ssb_no_signal_class_appended(self) -> None:
        """nDigenome SSB rows do not get signal class appended."""
        settings = make_settings("ndigenome")
        row = {
            "softclip_fraction": 0.0,
            "indel_fraction": 0.0,
            "known_indel_overlap": "",
            "support_mean_mapq": 60.0,
            "control_status": "UNCONTROLLED",
            "signal_classification": "SSB",
        }
        apply_filters(row, settings)
        self.assertEqual(row["filter_reasons"], "")
        self.assertEqual(row["filter_status"], "PASS")


class FormatValueTests(unittest.TestCase):
    def test_format_value_none(self) -> None:
        """format_value(None) → ""."""
        self.assertEqual(format_value(None), "")

    def test_format_value_float(self) -> None:
        """format_value(0.1) → "0.1"."""
        self.assertEqual(format_value(0.1), "0.1")

    def test_format_value_fraction(self) -> None:
        """format_value(1/3) → "0.33333333"."""
        result = format_value(1 / 3)
        self.assertEqual(result, "0.33333333")

    def test_format_value_int(self) -> None:
        """format_value(5) → "5"."""
        self.assertEqual(format_value(5), "5")

    def test_format_value_string(self) -> None:
        """format_value("x") → "x"."""
        self.assertEqual(format_value("x"), "x")


class BedLineTests(unittest.TestCase):
    def test_bed_line_digenome(self) -> None:
        """Digenome row at forward position 500 with combined_fraction 0.5 → BED line."""
        settings = make_settings("digenome")
        row = {
            "sample": "Sample",
            "forward_position_0based": 500,
            "combined_fraction": 0.5,
            "contig": "chr1",
        }
        line = bed_line(row, settings)
        self.assertEqual(line, "chr1\t500\t501\tSample|both|DSB\t500\t.\n")

    def test_bed_line_ndigenome(self) -> None:
        """nDigenome "-" row at 100 with endpoint_fraction 0.9996 → score 1000."""
        settings = make_settings("ndigenome")
        row = {
            "sample": "Sample",
            "position_0based": 100,
            "strand": "-",
            "endpoint_fraction": 0.9996,
            "contig": "chr1",
            "signal_classification": "SSB",
        }
        line = bed_line(row, settings)
        self.assertEqual(line, "chr1\t100\t101\tSample|-|SSB\t1000\t-\n")


class ColumnTests(unittest.TestCase):
    def test_digenome_columns_count(self) -> None:
        """DIGENOME_COLUMNS has 50 unique names."""
        self.assertEqual(len(DIGENOME_COLUMNS), 50)
        self.assertEqual(len(set(DIGENOME_COLUMNS)), 50)

    def test_ndigenome_columns_count(self) -> None:
        """NDIGENOME_COLUMNS has 39."""
        self.assertEqual(len(NDIGENOME_COLUMNS), 39)
        self.assertEqual(len(set(NDIGENOME_COLUMNS)), 39)

    def test_no_classification_column(self) -> None:
        """Neither DIGENOME_COLUMNS nor NDIGENOME_COLUMNS contains "classification"."""
        self.assertNotIn("classification", DIGENOME_COLUMNS)
        self.assertNotIn("classification", NDIGENOME_COLUMNS)

    def test_no_old_digenome_score_columns(self) -> None:
        """No "digenome_score" or "control_digenome_score" columns."""
        self.assertNotIn("digenome_score", DIGENOME_COLUMNS)
        self.assertNotIn("control_digenome_score", DIGENOME_COLUMNS)
        self.assertNotIn("digenome_score", NDIGENOME_COLUMNS)
        self.assertNotIn("control_digenome_score", NDIGENOME_COLUMNS)

    def test_ndigenome_no_rgen_or_forward_columns(self) -> None:
        """NDIGENOME_COLUMNS contains no "rgen" or "forward_" columns."""
        for col in NDIGENOME_COLUMNS:
            self.assertNotIn("rgen", col.lower())
            self.assertFalse(col.startswith("forward_"))


if __name__ == "__main__":
    unittest.main()
