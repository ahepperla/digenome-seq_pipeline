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
  simd=${FAKE_BWA_SIMD:-avx2}
  version=${FAKE_BWA_VERSION:-2.3-test}
  echo "Looking to launch executable \\"/opt/bwa-mem2-${version}_x64-linux/bwa-mem2.${simd}\\", simd = .${simd}"
  echo "Launching executable \\"/opt/bwa-mem2-${version}_x64-linux/bwa-mem2.${simd}\\""
  echo "${version}"
  for ((i = 0; i < 10000; i++)); do
    echo "additional version output $i"
  done
  exit 0
fi
if [[ "$1" == "index" && "$2" == "-p" ]]; then
  prefix=$3
  if [[ -n "${FAKE_BWA_STARTED:-}" ]]; then
    : > "$FAKE_BWA_STARTED"
  fi
  if [[ -n "${FAKE_BWA_RELEASE:-}" ]]; then
    while [[ ! -e "$FAKE_BWA_RELEASE" ]]; do
      sleep 0.05
    done
  fi
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

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def helper_environment(self, simd: str = "avx2", version: str = "2.3-test") -> dict[str, str]:
        environment = os.environ.copy()
        environment["PATH"] = f"{self.bin_dir}:{environment['PATH']}"
        environment["FAKE_BWA_COUNT"] = str(self.count_file)
        environment["FAKE_BWA_SIMD"] = simd
        environment["FAKE_BWA_VERSION"] = version
        return environment

    def helper_command(self, stale_seconds: int = 172800) -> list[str]:
        return [
            str(INDEX_HELPER),
            "--genome", "tiny",
            "--fasta", str(self.fasta),
            "--cache-dir", str(self.cache),
            "--lock-timeout-seconds", "1",
            "--stale-lock-seconds", str(stale_seconds),
        ]

    def run_helper(
        self,
        stale_seconds: int = 172800,
        simd: str = "avx2",
        version: str = "2.3-test",
        umask: int | None = None,
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            self.helper_command(stale_seconds),
            text=True,
            capture_output=True,
            env=self.helper_environment(simd, version),
            preexec_fn=(lambda: os.umask(umask)) if umask is not None else None,
            check=False,
        )

    def prefix(self, **options) -> Path:
        """Run the helper, require success, and return the printed prefix."""
        result = self.run_helper(**options)
        self.assertEqual(result.returncode, 0, result.stderr)
        return Path(result.stdout.strip())

    def builds(self) -> int:
        return int(self.count_file.read_text())

    def test_stdout_is_only_the_index_prefix(self) -> None:
        result = self.run_helper()
        self.assertEqual(result.returncode, 0, result.stderr)
        digest = hashlib.sha256(self.fasta.read_bytes()).hexdigest()
        expected = self.cache.resolve() / "tiny" / digest / "2.3-test" / "bwamem2" / "tiny"
        self.assertEqual(result.stdout, f"{expected}\n")

    def test_complete_index_is_reused(self) -> None:
        first = self.prefix()
        manifest = first.parent / "index.complete.tsv"
        manifest_mtime = manifest.stat().st_mtime_ns
        self.assertEqual(self.prefix(), first)
        self.assertEqual(self.builds(), 1)
        self.assertEqual(manifest.stat().st_mtime_ns, manifest_mtime)
        self.assertIn("bwa_mem2_version\t2.3-test\n", manifest.read_text())
        for suffix in SUFFIXES:
            self.assertEqual(Path(f"{first}{suffix}").read_text(), "index\n")

    def test_shared_cache_permissions_override_restrictive_umask(self) -> None:
        prefix = self.prefix(umask=0o077)
        final_directory = prefix.parent
        version_directory = final_directory.parent
        fingerprint_directory = version_directory.parent
        genome_directory = fingerprint_directory.parent
        expected_directory_modes = {
            self.cache: 0o1777,
            genome_directory: 0o1777,
            fingerprint_directory: 0o1777,
            version_directory: 0o1777,
            final_directory: 0o755,
        }
        for path, expected_mode in expected_directory_modes.items():
            with self.subTest(path=path):
                self.assertEqual(path.stat().st_mode & 0o7777, expected_mode)
        for path in final_directory.iterdir():
            with self.subTest(path=path):
                self.assertTrue(path.is_file())
                self.assertEqual(path.stat().st_mode & 0o7777, 0o644)
        self.assertEqual(self.prefix(version="2.4-test", umask=0o077).parent.parent.name, "2.4-test")

    def test_setgid_modes_satisfy_production_mode_check(self) -> None:
        helper_text = INDEX_HELPER.read_text()
        start = helper_text.index("mode_satisfies_requirement() {")
        end = helper_text.index("\n}\n", start) + 2
        function_text = helper_text[start:end]
        for current, required, expected in (
            ("3777", "1777", True),
            ("2755", "755", True),
            ("1777", "1777", True),
            ("0777", "1777", False),
            ("1775", "1777", False),
            ("5777", "1777", False),
        ):
            with self.subTest(current=current, required=required):
                result = subprocess.run(
                    ["bash", "-c", f"{function_text}\nmode_satisfies_requirement {current} {required}"],
                    check=False,
                )
                self.assertEqual(result.returncode == 0, expected)

    def test_build_lock_is_readable_with_restrictive_umask(self) -> None:
        started = self.tmp / "build.started"
        release = self.tmp / "build.release"
        environment = self.helper_environment()
        environment["FAKE_BWA_STARTED"] = str(started)
        environment["FAKE_BWA_RELEASE"] = str(release)
        process = subprocess.Popen(
            self.helper_command(),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=environment,
            preexec_fn=lambda: os.umask(0o077),
        )
        digest = hashlib.sha256(self.fasta.read_bytes()).hexdigest()
        lock = self.cache / "tiny" / f".{digest}.2.3-test.build.lock"
        owner = lock / "owner.tsv"
        try:
            deadline = time.time() + 5
            while time.time() < deadline and not started.exists():
                time.sleep(0.05)
            self.assertTrue(started.exists())
            self.assertTrue(owner.is_file())
            self.assertEqual(lock.stat().st_mode & 0o7777, 0o755)
            self.assertEqual(owner.stat().st_mode & 0o7777, 0o644)
        finally:
            release.touch()
            stdout, stderr = process.communicate(timeout=5)
        self.assertEqual(process.returncode, 0, stderr or stdout)

    def test_cpu_dispatch_message_does_not_invalidate_index(self) -> None:
        first = self.prefix(simd="avx2")
        self.assertEqual(self.prefix(simd="avx512"), first)
        self.assertEqual(self.builds(), 1)
        self.assertEqual(first.parent.parent.name, "2.3-test")

    def test_changed_fasta_gets_a_new_content_addressed_index(self) -> None:
        first = self.prefix()
        self.fasta.write_text(">chr1\nACGTACGTACGTAAAA\n")
        second = self.prefix()
        self.assertNotEqual(first, second)
        self.assertEqual(self.builds(), 2)

    def test_bwa_versions_use_isolated_cache_directories(self) -> None:
        first = self.prefix(version="2.3-test")
        second = self.prefix(version="2.4-test")
        self.assertEqual(first.parent.parent.name, "2.3-test")
        self.assertEqual(second.parent.parent.name, "2.4-test")
        self.assertEqual(first.parent.parent.parent, second.parent.parent.parent)
        self.assertEqual(self.builds(), 2)
        self.assertEqual(self.prefix(version="2.3-test"), first)
        self.assertEqual(self.builds(), 2)

    def test_partial_index_is_quarantined_and_rebuilt(self) -> None:
        prefix = self.prefix()
        Path(f"{prefix}.pac").unlink()
        self.assertEqual(self.prefix(), prefix)
        self.assertEqual(self.builds(), 2)
        quarantined = list(prefix.parent.parent.glob("bwamem2.incomplete.*"))
        self.assertEqual(len(quarantined), 1)
        self.assertFalse((quarantined[0] / "tiny.pac").exists())

    def test_stale_local_dead_pid_lock_is_recovered(self) -> None:
        digest = hashlib.sha256(self.fasta.read_bytes()).hexdigest()
        lock = self.cache / "tiny" / f".{digest}.2.3-test.build.lock"
        lock.mkdir(parents=True)
        (lock / "owner.tsv").write_text(
            f"hostname\t{os.uname().nodename}\n"
            "pid\t99999999\n"
            f"created_epoch\t{int(time.time()) - 100}\n"
        )
        self.prefix(stale_seconds=1)
        self.assertFalse(lock.exists())

    def test_live_lock_from_another_host_is_not_removed(self) -> None:
        digest = hashlib.sha256(self.fasta.read_bytes()).hexdigest()
        lock = self.cache / "tiny" / f".{digest}.2.3-test.build.lock"
        lock.mkdir(parents=True)
        (lock / "owner.tsv").write_text(
            "hostname\tanother-node\n"
            "pid\t1\n"
            f"created_epoch\t{int(time.time()) - 100}\n"
        )
        result = self.run_helper(stale_seconds=1)
        self.assertEqual(result.returncode, 1)
        self.assertIn("stale or unverifiable index lock detected", result.stderr)
        self.assertTrue(lock.exists())


if __name__ == "__main__":
    unittest.main()
