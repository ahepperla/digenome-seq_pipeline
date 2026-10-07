"""Tests for statistics: Fisher exact test, Benjamini-Hochberg, control comparison."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from cleavage.stats import (  # noqa: E402
    fisher_exact_two_sided,
    benjamini_hochberg,
    compare_to_control,
    MATCHED_CONTROL,
    INSUFFICIENT_CONTROL,
)


class FisherExactTests(unittest.TestCase):
    def test_fisher_exact_2x2_table(self) -> None:
        p_value = fisher_exact_two_sided(3, 1, 1, 3)
        self.assertAlmostEqual(p_value, 17.0 / 35.0, places=12)

    def test_fisher_exact_all_zeros_returns_one(self) -> None:
        p_value = fisher_exact_two_sided(0, 0, 0, 0)
        self.assertEqual(p_value, 1.0)

    def test_fisher_exact_is_symmetric_under_row_swap(self) -> None:
        p1 = fisher_exact_two_sided(1, 3, 3, 1)
        p2 = fisher_exact_two_sided(3, 1, 1, 3)
        self.assertAlmostEqual(p1, p2, places=12)


class BenjaminiHochbergTests(unittest.TestCase):
    def test_benjamini_hochberg_example(self) -> None:
        p_values = [0.01, 0.04, 0.03, 0.20]
        q_values = benjamini_hochberg(p_values)
        self.assertAlmostEqual(q_values[0], 0.04, places=15)
        self.assertAlmostEqual(q_values[1], 0.04 * 4 / 3, places=15)
        self.assertAlmostEqual(q_values[2], 0.04 * 4 / 3, places=15)
        self.assertAlmostEqual(q_values[3], 0.2, places=15)

    def test_benjamini_hochberg_tied_p_values(self) -> None:
        p_values = [0.02, 0.02, 0.5]
        q_values = benjamini_hochberg(p_values)
        self.assertEqual(q_values[0], q_values[1])
        self.assertAlmostEqual(q_values[0], 0.02 * 3 / 2, places=15)

    def test_benjamini_hochberg_empty_list(self) -> None:
        q_values = benjamini_hochberg([])
        self.assertEqual(q_values, [])

    def test_benjamini_hochberg_monotonicity(self) -> None:
        p_values = [0.001, 0.01, 0.05, 0.1, 0.2]
        q_values = benjamini_hochberg(p_values)
        for i in range(len(q_values) - 1):
            self.assertLessEqual(q_values[i], q_values[i + 1])


class CompareToControlTests(unittest.TestCase):
    def test_insufficient_control_depth_returns_status_and_none_values(self) -> None:
        status, fold, p = compare_to_control(10, 20, 0, 0, min_control_depth=5)
        self.assertEqual(status, INSUFFICIENT_CONTROL)
        self.assertIsNone(fold)
        self.assertIsNone(p)

    def test_control_depth_equal_to_minimum_returns_matched_status(self) -> None:
        status, fold, p = compare_to_control(10, 20, 5, 20, min_control_depth=20)
        self.assertEqual(status, MATCHED_CONTROL)
        self.assertIsNotNone(fold)
        self.assertIsNotNone(p)

    def test_fold_enrichment_formula(self) -> None:
        status, fold, p = compare_to_control(10, 20, 0, 30, min_control_depth=1)
        self.assertEqual(status, MATCHED_CONTROL)
        treated_rate = (10 + 0.5) / (20 + 1.0)
        control_rate = (0 + 0.5) / (30 + 1.0)
        expected_fold = treated_rate / control_rate
        self.assertAlmostEqual(fold, expected_fold, places=15)


if __name__ == "__main__":
    unittest.main()
