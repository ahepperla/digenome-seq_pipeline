from __future__ import annotations

import hashlib
import os
import subprocess
import tempfile
import time
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INDEX_HELPER = ROOT / "bin" / "prepare_bwamem2_index.sh"
SUFFIXES = [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]


class IndexCacheTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tempdir.name)
        self.bin_dir = self.tmp / "bin"
        self.bin_dir.mkdir()
        self.count_file = self.tmp / "build_count.txt"
        fake_bwa = self.bin_dir / "bwa-mem2"
        fake_bwa.write_text(
            """#!/usr/bin/env bash
set -euo pipefail
if [[ "$1" == "version" ]]; then
  echo "2.3-test"
  for ((i = 0; i < 10000; i++)); do
    echo "additional version output $i"
  done
  exit 0
fi
if [[ "$1" == "index" && "$2" == "-p" ]]; then
  prefix=$3
  for suffix in .0123 .amb .ann .bwt.2bit.64 .pac; do
    printf 'index\\n' > "${prefix}${suffix}"
  done
  count=0
  [[ ! -s "$FAKE_BWA_COUNT" ]] || count=$(cat "$FAKE_BWA_COUNT")
  printf '%s\\n' "$((count + 1))" > "$FAKE_BWA_COUNT"
  exit 0
fi
exit 2
"""
        )
        fake_bwa.chmod(0o755)
        self.fasta = self.tmp / "tiny.fa"
        self.fasta.write_text(">chr1\nACGTACGTACGT\n")
        self.cache = self.tmp / "cache"
        self.ready = self.tmp / "ready.tsv"

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def run_helper(self, stale_seconds: int = 172800) -> subprocess.CompletedProcess[str]:
        environment = os.environ.copy()
        environment["PATH"] = f"{self.bin_dir}:{environment['PATH']}"
        environment["FAKE_BWA_COUNT"] = str(self.count_file)
        return subprocess.run(
            [
                str(INDEX_HELPER),
                "--genome",
                "tiny",
                "--fasta",
                str(self.fasta),
                "--cache-dir",
                str(self.cache),
                "--ready",
                str(self.ready),
                "--lock-timeout-seconds",
                "1",
                "--stale-lock-seconds",
                str(stale_seconds),
            ],
            text=True,
            capture_output=True,
            env=environment,
            check=False,
        )

    def ready_values(self) -> dict[str, str]:
        return dict(
            line.split("\t", 1)
            for line in self.ready.read_text().splitlines()
            if "\t" in line
        )

    def test_complete_index_is_reused(self) -> None:
        first = self.run_helper()
        self.assertEqual(first.returncode, 0, first.stderr)
        second = self.run_helper()
        self.assertEqual(second.returncode, 0, second.stderr)
        self.assertEqual(self.count_file.read_text().strip(), "1")
        prefix = self.ready_values()["index_prefix"]
        self.assertEqual(self.ready_values()["bwa_mem2_version"], "2.3-test")
        for suffix in SUFFIXES:
            self.assertTrue(Path(prefix + suffix).stat().st_size > 0)

    def test_changed_fasta_gets_a_new_content_addressed_index(self) -> None:
        self.assertEqual(self.run_helper().returncode, 0)
        first_prefix = self.ready_values()["index_prefix"]
        self.fasta.write_text(">chr1\nACGTACGTACGTAAAA\n")
        self.assertEqual(self.run_helper().returncode, 0)
        second_prefix = self.ready_values()["index_prefix"]
        self.assertNotEqual(first_prefix, second_prefix)
        self.assertEqual(self.count_file.read_text().strip(), "2")

    def test_partial_index_is_quarantined_and_rebuilt(self) -> None:
        self.assertEqual(self.run_helper().returncode, 0)
        prefix = self.ready_values()["index_prefix"]
        Path(prefix + ".pac").unlink()
        result = self.run_helper()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.count_file.read_text().strip(), "2")
        quarantined = list(Path(prefix).parent.parent.glob("bwamem2.incomplete.*"))
        self.assertEqual(len(quarantined), 1)

    def test_stale_local_dead_pid_lock_is_recovered(self) -> None:
        digest = hashlib.sha256(self.fasta.read_bytes()).hexdigest()
        lock = self.cache / "tiny" / f".{digest}.build.lock"
        lock.mkdir(parents=True)
        (lock / "owner.tsv").write_text(
            f"hostname\t{os.uname().nodename}\n"
            "pid\t99999999\n"
            f"created_epoch\t{int(time.time()) - 100}\n"
        )
        result = self.run_helper(stale_seconds=1)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(lock.exists())


if __name__ == "__main__":
    unittest.main()
