#!/usr/bin/env python3
"""Write golden outputs from the pre-refactor caller (commit b9455e2).

The old `bin/` is extracted from git into a temporary directory, so this
script keeps working after the old files are deleted. It runs the old serial
caller on every scenario in both modes and writes the five call files, the
MultiQC table, and the QC counts to tests/golden/legacy/.

Run from the repository root:  python3 tests/golden/build_legacy.py
"""

from __future__ import annotations

import argparse
import importlib
import io
import json
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LEGACY_COMMIT = "b9455e2"
OUTPUT_DIR = HERE / "legacy"
CALL_FILES = ["all.tsv", "high_confidence.tsv", "manual_review.tsv", "artifact.tsv", "bed"]
QC_COUNTS = [
    "candidate_endpoints_or_pairs_before_filters",
    "reported_candidates",
    "high_confidence_calls",
    "manual_review_candidates",
    "artifact_candidates",
    "classifications",
    "signal_classifications",
]
# Old caller option names for the pipeline's cleavage_* parameters.
LEGACY_NAMES = {
    "cleavage_artifact_window": "artifact_window",
    "cleavage_max_softclip_fraction": "max_softclip_fraction",
    "cleavage_max_indel_fraction": "max_indel_fraction",
    "cleavage_min_support_mean_mapq": "min_support_mean_mapq",
    "cleavage_control_min_depth": "control_min_depth",
    "cleavage_control_max_fraction": "control_max_fraction",
    "cleavage_control_min_fold": "control_min_fold",
    "cleavage_control_max_q": "control_max_q",
}

sys.path.insert(0, str(HERE))
from scenarios import all_scenarios  # noqa: E402


def load_legacy_caller(directory: Path):
    """Extract bin/ from LEGACY_COMMIT and import its call_cleavage module."""
    archive = subprocess.run(
        ["git", "-C", str(ROOT), "archive", LEGACY_COMMIT, "bin"],
        check=True,
        capture_output=True,
    ).stdout
    with tarfile.open(fileobj=io.BytesIO(archive)) as tar:
        tar.extractall(directory, filter="data")
    sys.path.insert(0, str(directory / "bin"))
    return importlib.import_module("call_cleavage")


def legacy_args(settings, paths, analysis: str, prefix: Path) -> argparse.Namespace:
    options = {LEGACY_NAMES.get(name, name): value for name, value in settings.items()}
    control = paths["control_bam"]
    return argparse.Namespace(
        analysis=analysis,
        bam=str(paths["bam"]),
        sample="Sample",
        output_prefix=str(prefix),
        control_bam=str(control) if control else None,
        control_sample="Control" if control else "",
        variant_vcf=str(paths["vcf"]) if paths["vcf"] else None,
        genome_blacklist=str(paths["blacklist"]) if paths["blacklist"] else None,
        intervals_file=None,
        chunk_id=None,
        raw_output=None,
        summary_output=None,
        **options,
    )


def main() -> None:
    if OUTPUT_DIR.exists():
        shutil.rmtree(OUTPUT_DIR)
    OUTPUT_DIR.mkdir()
    with tempfile.TemporaryDirectory() as temporary:
        workspace = Path(temporary)
        cleavage = load_legacy_caller(workspace / "legacy")
        for scenario in all_scenarios():
            paths = scenario.write_inputs(workspace / scenario.name)
            for analysis in ("digenome", "ndigenome"):
                prefix = workspace / scenario.name / "out"
                qc = cleavage.run_serial_calling(
                    legacy_args(scenario.settings, paths, analysis, prefix)
                )
                stem = f"{scenario.name}.{analysis}"
                for suffix in CALL_FILES:
                    shutil.copy(f"{prefix}.{analysis}.{suffix}", OUTPUT_DIR / f"{stem}.{suffix}")
                shutil.copy(f"{prefix}.{analysis}_mqc.tsv", OUTPUT_DIR / f"{stem}_mqc.tsv")
                counts = {key: qc[key] for key in QC_COUNTS}
                (OUTPUT_DIR / f"{stem}.counts.json").write_text(
                    json.dumps(counts, indent=2, sort_keys=True) + "\n"
                )
                print(f"{stem}: {qc['reported_candidates']} row(s)")


if __name__ == "__main__":
    main()
