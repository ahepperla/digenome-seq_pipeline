#!/usr/bin/env python3
"""Write the golden outputs in tests/golden/expected/ from the current caller.

These files are the regression anchor for tests/test_golden.py. Rewriting them
changes what counts as correct, so do it only with the project lead's
decision (docs/decisions.md).

Run from the repository root:  python3 tests/golden/build_expected.py
"""

from __future__ import annotations

import json
import shutil
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUTPUT_DIR = HERE / "expected"
OUTPUT_FILES = ["all.tsv", "high_confidence.tsv", "manual_review.tsv", "artifact.tsv", "bed"]

sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
from helpers import make_settings, run_caller  # noqa: E402
from scenarios import all_scenarios  # noqa: E402


def comparable_qc(qc: dict) -> dict:
    """The QC report without the chunk count, the one field that depends on chunking."""
    return {key: value for key, value in qc.items() if key != "chunks"}


def main() -> None:
    if OUTPUT_DIR.exists():
        shutil.rmtree(OUTPUT_DIR)
    OUTPUT_DIR.mkdir()
    with tempfile.TemporaryDirectory() as temporary:
        workspace = Path(temporary)
        for scenario in all_scenarios():
            inputs = scenario.write_inputs(workspace / scenario.name)
            for analysis in ("digenome", "ndigenome"):
                settings = make_settings(analysis, base=scenario.settings)
                prefix, qc = run_caller(settings, inputs, workspace / f"{scenario.name}_{analysis}", chunks=1)
                stem = f"{scenario.name}.{analysis}"
                for suffix in OUTPUT_FILES:
                    shutil.copy(f"{prefix}.{analysis}.{suffix}", OUTPUT_DIR / f"{stem}.{suffix}")
                shutil.copy(f"{prefix}.{analysis}_mqc.tsv", OUTPUT_DIR / f"{stem}_mqc.tsv")
                (OUTPUT_DIR / f"{stem}.qc.json").write_text(json.dumps(comparable_qc(qc), indent=2) + "\n")
                print(f"{stem}: {qc['rows']} row(s)")


if __name__ == "__main__":
    main()
