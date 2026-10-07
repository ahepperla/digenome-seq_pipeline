"""Every output of every golden scenario must match tests/golden/expected/
exactly, whatever the chunk count.

The expected files were written by tests/golden/build_expected.py after the
rebuilt caller was shown to reproduce every value of the pre-refactor caller
(commit b9455e2; see docs/decisions.md). Changing them needs the project lead.
"""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "golden"))
from build_expected import OUTPUT_DIR, OUTPUT_FILES, comparable_qc  # noqa: E402
from helpers import make_settings, run_caller  # noqa: E402
from scenarios import all_scenarios  # noqa: E402

CHUNK_COUNTS = (1, 2, 4)


class GoldenOutputTests(unittest.TestCase):
    def test_outputs_match_the_golden_files(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary)
            for scenario in all_scenarios():
                inputs = scenario.write_inputs(workspace / scenario.name)
                for analysis in ("digenome", "ndigenome"):
                    settings = make_settings(analysis, base=scenario.settings)
                    expected = OUTPUT_DIR / f"{scenario.name}.{analysis}"
                    for chunks in CHUNK_COUNTS:
                        with self.subTest(scenario=scenario.name, analysis=analysis, chunks=chunks):
                            out_dir = workspace / f"{scenario.name}_{analysis}_{chunks}"
                            prefix, qc = run_caller(settings, inputs, out_dir, chunks)
                            for suffix in OUTPUT_FILES:
                                self.assertEqual(
                                    Path(f"{prefix}.{analysis}.{suffix}").read_text(),
                                    Path(f"{expected}.{suffix}").read_text(),
                                    suffix,
                                )
                            self.assertEqual(
                                Path(f"{prefix}.{analysis}_mqc.tsv").read_text(),
                                Path(f"{expected}_mqc.tsv").read_text(),
                            )
                            self.assertEqual(
                                comparable_qc(qc),
                                json.loads(Path(f"{expected}.qc.json").read_text()),
                            )


if __name__ == "__main__":
    unittest.main()
