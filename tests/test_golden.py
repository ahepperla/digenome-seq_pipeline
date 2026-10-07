"""The refactored caller must reproduce the pre-refactor caller's calls.

tests/golden/legacy/ holds the outputs of the old serial caller (commit
b9455e2) for every scenario in tests/golden/scenarios.py. Each old row is
mapped to the new column layout (decision "Streamlined, mode-specific output
columns" in docs/decisions.md) and rows are reordered by BAM header contig
order. Every value must then match exactly, at 1, 2, and 4 chunks.
"""

from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "golden"))
from helpers import run_caller  # noqa: E402
from scenarios import all_scenarios  # noqa: E402

from cleavage.output import DIGENOME_COLUMNS, NDIGENOME_COLUMNS  # noqa: E402
from cleavage.settings import settings_from_dict  # noqa: E402

LEGACY = HERE / "golden" / "legacy"
CHUNK_COUNTS = (1, 2, 4)
# New Digenome column -> old column, where the name changed.
RENAMED = {
    "combined_endpoint_count": "endpoint_count",
    "combined_depth": "strand_depth",
    "combined_fraction": "endpoint_fraction",
    "control_combined_endpoint_count": "control_endpoint_count",
    "control_combined_depth": "control_depth",
    "control_combined_fraction": "control_fraction",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def legacy_tier(row: dict[str, str]) -> str:
    if row["filter_status"] == "PASS":
        return "high_confidence"
    return "artifact" if row["classification"] == "ARTIFACT_RISK" else "manual_review"


def expected_rows(path: Path, analysis: str, contig_order: dict[str, int]) -> list[list[str]]:
    """Old rows in the new layout and order."""
    columns = DIGENOME_COLUMNS if analysis == "digenome" else NDIGENOME_COLUMNS
    position = "forward_position_0based" if analysis == "digenome" else "position_0based"
    rows = read_tsv(path)
    rows.sort(key=lambda row: (contig_order[row["contig"]], int(row[position]), row["strand"]))
    mapped = []
    for row in rows:
        row["tier"] = legacy_tier(row)
        mapped.append([row[RENAMED.get(column, column)] for column in columns])
    return mapped


def expected_bed(path: Path, contig_order: dict[str, int]) -> list[str]:
    lines = path.read_text().splitlines()
    return sorted(lines, key=lambda line: (contig_order[line.split("\t")[0]], int(line.split("\t")[1])))


@unittest.skipUnless(LEGACY.is_dir(), "golden outputs are missing")
class GoldenOutputTests(unittest.TestCase):
    def test_outputs_match_the_pre_refactor_caller(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary)
            for scenario in all_scenarios():
                inputs = scenario.write_inputs(workspace / scenario.name)
                contig_order = {name: index for index, (name, _length) in enumerate(scenario.contigs)}
                for analysis in ("digenome", "ndigenome"):
                    settings = settings_from_dict({**scenario.settings, "analysis": analysis})
                    for chunks in CHUNK_COUNTS:
                        with self.subTest(scenario=scenario.name, analysis=analysis, chunks=chunks):
                            out_dir = workspace / f"{scenario.name}_{analysis}_{chunks}"
                            prefix, qc = run_caller(settings, inputs, out_dir, chunks)
                            self.check_outputs(scenario.name, analysis, prefix, qc, contig_order)

    def check_outputs(self, name, analysis, prefix, qc, contig_order) -> None:
        legacy_stem = LEGACY / f"{name}.{analysis}"
        for tier in ("all", "high_confidence", "manual_review", "artifact"):
            with Path(f"{prefix}.{analysis}.{tier}.tsv").open(newline="") as handle:
                observed = list(csv.reader(handle, delimiter="\t"))[1:]
            expected = expected_rows(Path(f"{legacy_stem}.{tier}.tsv"), analysis, contig_order)
            self.assertEqual(observed, expected, tier)
        self.assertEqual(
            Path(f"{prefix}.{analysis}.bed").read_text().splitlines(),
            expected_bed(Path(f"{legacy_stem}.bed"), contig_order),
        )
        counts = json.loads(Path(f"{legacy_stem}.counts.json").read_text())
        self.assertEqual(qc["candidates_before_filters"], counts["candidate_endpoints_or_pairs_before_filters"])
        self.assertEqual(qc["rows"], counts["reported_candidates"])
        self.assertEqual(qc["tiers"]["high_confidence"], counts["high_confidence_calls"])
        self.assertEqual(qc["tiers"]["manual_review"], counts["manual_review_candidates"])
        self.assertEqual(qc["tiers"]["artifact"], counts["artifact_candidates"])
        if analysis == "ndigenome":
            self.assertEqual(qc["signal_classifications"], counts["signal_classifications"])


if __name__ == "__main__":
    unittest.main()
