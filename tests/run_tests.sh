#!/usr/bin/env bash
# The full local check suite. Run from anywhere, with a python3 that has
# pysam (the main checkout's .venv). Nextflow lint and stub runs are part of
# the unittest suite and are skipped when `nextflow` is not installed.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
bash -n bin/prepare_bwamem2_index.sh containers/build_cleavage.sh tests/run_tests.sh
python3 -m unittest discover -s tests -p 'test_*.py' "$@"
