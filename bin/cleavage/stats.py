"""Matched-control statistics: fold enrichment, Fisher's exact test, and
Benjamini-Hochberg q-values. Implemented here to avoid a SciPy dependency."""

from __future__ import annotations

import math

UNCONTROLLED = "UNCONTROLLED"
MATCHED_CONTROL = "MATCHED_CONTROL"
INSUFFICIENT_CONTROL = "INSUFFICIENT_CONTROL_COVERAGE"


def log_choose(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float:
    """Two-sided p-value for the 2x2 table [[a, b], [c, d]]."""
    row_one, row_two, column_one = a + b, c + d, a + c
    total = row_one + row_two
    if total == 0:
        return 1.0

    def probability(cell: int) -> float:
        return math.exp(
            log_choose(row_one, cell)
            + log_choose(row_two, column_one - cell)
            - log_choose(total, column_one)
        )

    observed = probability(a)
    tables = [probability(cell) for cell in range(max(0, column_one - row_two), min(row_one, column_one) + 1)]
    # The relative tolerance keeps tables tied with the observed one.
    return min(1.0, sum(p for p in tables if p <= observed * (1.0 + 1e-12)))


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    """Step-up adjusted q-values, in the order of `p_values`."""
    total = len(p_values)
    ranked = sorted(range(total), key=lambda index: p_values[index])
    q_values = [1.0] * total
    running = 1.0
    for rank in range(total, 0, -1):
        index = ranked[rank - 1]
        running = min(running, p_values[index] * total / rank)
        q_values[index] = min(1.0, running)
    return q_values


def compare_to_control(
    treated_count: int,
    treated_depth: int,
    control_count: int,
    control_depth: int,
    min_control_depth: int,
) -> tuple[str, float | None, float | None]:
    """Return (control status, fold enrichment, Fisher p-value).

    Below `min_control_depth` the comparison would be misleading, so fold and
    p are left blank and the row is marked for filtering.
    """
    if control_depth < min_control_depth:
        return INSUFFICIENT_CONTROL, None, None
    p_value = fisher_exact_two_sided(
        treated_count,
        max(0, treated_depth - treated_count),
        control_count,
        max(0, control_depth - control_count),
    )
    # 0.5 pseudocount: a site absent from the control still gets a finite fold.
    treated_rate = (treated_count + 0.5) / (treated_depth + 1.0)
    control_rate = (control_count + 0.5) / (control_depth + 1.0)
    return MATCHED_CONTROL, treated_rate / control_rate, p_value
