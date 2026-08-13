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
        self.ready = self.tmp / "ready.tsv"

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def helper_environment(
        self,
        simd: str = "avx2",
        version: str = "2.3-test",
    ) -> dict[str, str]:
        environment = os.environ.copy()
        environment["PATH"] = f"{self.bin_dir}:{environment['PATH']}"
        environment["FAKE_BWA_COUNT"] = str(self.count_file)
        environment["FAKE_BWA_SIMD"] = simd
        environment["FAKE_BWA_VERSION"] = version
        return environment

    def helper_command(
        self,
        stale_seconds: int = 172800,
    ) -> list[str]:
        return [
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
            preexec_fn=(
                (lambda: os.umask(umask))
                if umask is not None
                else None
            ),
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
        first_ready_mtime = self.ready.stat().st_mtime_ns
        second = self.run_helper()
        self.assertEqual(second.returncode, 0, second.stderr)
        self.assertEqual(self.count_file.read_text().strip(), "1")
        self.assertEqual(self.ready.stat().st_mtime_ns, first_ready_mtime)
        prefix = self.ready_values()["index_prefix"]
        self.assertEqual(self.ready_values()["bwa_mem2_version"], "2.3-test")
        for suffix in SUFFIXES:
            self.assertTrue(Path(prefix + suffix).stat().st_size > 0)

    def test_shared_cache_permissions_override_restrictive_umask(
        self,
    ) -> None:
        result = self.run_helper(umask=0o077)
        self.assertEqual(result.returncode, 0, result.stderr)

        prefix = Path(self.ready_values()["index_prefix"])
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
                self.assertEqual(
                    path.stat().st_mode & 0o7777,
                    expected_mode,
                )

        for path in final_directory.iterdir():
            with self.subTest(path=path):
                self.assertTrue(path.is_file())
                self.assertEqual(path.stat().st_mode & 0o7777, 0o644)

        second_version = self.run_helper(
            version="2.4-test",
            umask=0o077,
        )
        self.assertEqual(
            second_version.returncode,
            0,
            second_version.stderr,
        )
        self.assertIn(
            "/2.4-test/bwamem2/",
            self.ready_values()["index_prefix"],
        )

    def test_setgid_modes_satisfy_production_mode_check(self) -> None:
        helper_text = INDEX_HELPER.read_text()
        start = helper_text.index("mode_satisfies_requirement() {")
        end = helper_text.index(
            "\n}\n\nensure_directory_mode()",
            start,
        ) + 2
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
                    [
                        "bash",
                        "-c",
                        (
                            f"{function_text}\n"
                            "mode_satisfies_requirement "
                            f"{current} {required}"
                        ),
                    ],
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
        first = self.run_helper(simd="avx2")
        self.assertEqual(first.returncode, 0, first.stderr)
        second = self.run_helper(simd="avx512")
        self.assertEqual(second.returncode, 0, second.stderr)
        self.assertEqual(self.count_file.read_text().strip(), "1")
        self.assertEqual(self.ready_values()["bwa_mem2_version"], "2.3-test")

    def test_legacy_dispatch_manifest_is_migrated_without_rebuild(self) -> None:
        first = self.run_helper(simd="avx2")
        self.assertEqual(first.returncode, 0, first.stderr)

        prefix = Path(self.ready_values()["index_prefix"])
        manifest = prefix.parent / "index.complete.tsv"
        manifest_text = manifest.read_text()
        manifest.write_text(
            manifest_text.replace(
                "schema_version\t2",
                "schema_version\t1",
            ).replace(
                "bwa_mem2_version\t2.3-test",
                "bwa_mem2_version\tLooking to launch executable "
                '"/opt/bwa-mem2-2.3-test_x64-linux/bwa-mem2.avx2", '
                "simd = .avx2",
            )
        )

        second = self.run_helper(simd="avx512")
        self.assertEqual(second.returncode, 0, second.stderr)
        self.assertEqual(self.count_file.read_text().strip(), "1")
        self.assertIn("schema_version\t2", manifest.read_text())
        self.assertIn("bwa_mem2_version\t2.3-test", manifest.read_text())
        self.assertIn("version_metadata_migrated_utc\t", manifest.read_text())

    def test_changed_fasta_gets_a_new_content_addressed_index(self) -> None:
        self.assertEqual(self.run_helper().returncode, 0)
        first_prefix = self.ready_values()["index_prefix"]
        self.fasta.write_text(">chr1\nACGTACGTACGTAAAA\n")
        self.assertEqual(self.run_helper().returncode, 0)
        second_prefix = self.ready_values()["index_prefix"]
        self.assertNotEqual(first_prefix, second_prefix)
        self.assertEqual(self.count_file.read_text().strip(), "2")

    def test_bwa_versions_use_isolated_cache_directories(self) -> None:
        self.assertEqual(self.run_helper(version="2.3-test").returncode, 0)
        first_prefix = self.ready_values()["index_prefix"]
        self.assertEqual(self.run_helper(version="2.4-test").returncode, 0)
        second_prefix = self.ready_values()["index_prefix"]

        self.assertNotEqual(first_prefix, second_prefix)
        self.assertIn("/2.3-test/bwamem2/", first_prefix)
        self.assertIn("/2.4-test/bwamem2/", second_prefix)
        self.assertEqual(self.count_file.read_text().strip(), "2")
        for prefix in (first_prefix, second_prefix):
            self.assertTrue(
                all(Path(prefix + suffix).is_file() for suffix in SUFFIXES)
            )

        self.assertEqual(self.run_helper(version="2.3-test").returncode, 0)
        self.assertEqual(self.ready_values()["index_prefix"], first_prefix)
        self.assertEqual(self.count_file.read_text().strip(), "2")

    def test_unversioned_cache_is_migrated_without_rebuild(self) -> None:
        first = self.run_helper()
        self.assertEqual(first.returncode, 0, first.stderr)
        versioned_prefix = Path(self.ready_values()["index_prefix"])
        versioned_directory = versioned_prefix.parent
        fingerprint_directory = versioned_directory.parent.parent
        legacy_directory = fingerprint_directory / "bwamem2"
        versioned_directory.rename(legacy_directory)
        versioned_directory.parent.rmdir()

        second = self.run_helper()

        self.assertEqual(second.returncode, 0, second.stderr)
        self.assertIn(
            "Migrated unversioned bwa-mem2 index cache",
            second.stdout,
        )
        migrated_prefix = Path(self.ready_values()["index_prefix"])
        self.assertEqual(
            migrated_prefix,
            fingerprint_directory / "2.3-test" / "bwamem2" / "tiny",
        )
        self.assertEqual(self.count_file.read_text().strip(), "1")
        self.assertEqual(
            (legacy_directory / "tiny.0123").stat().st_ino,
            Path(f"{migrated_prefix}.0123").stat().st_ino,
        )

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
        lock = self.cache / "tiny" / f".{digest}.2.3-test.build.lock"
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
