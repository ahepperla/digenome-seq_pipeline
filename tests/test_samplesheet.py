from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VALIDATOR = ROOT / "bin" / "validate_samplesheet.py"


class SamplesheetTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def fastq(self, name: str, directory: Path | None = None) -> Path:
        target_dir = directory or self.tmp
        target_dir.mkdir(parents=True, exist_ok=True)
        path = target_dir / name
        path.touch()
        return path

    def run_validator(
        self, rows: list[dict[str, str]], analysis: str = "digenome"
    ) -> subprocess.CompletedProcess[str]:
        input_csv = self.tmp / "input.csv"
        output_csv = self.tmp / "output.csv"
        fieldnames = list(rows[0])
        with input_csv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        return subprocess.run(
            [
                sys.executable,
                str(VALIDATOR),
                str(input_csv),
                str(output_csv),
                analysis,
            ],
            text=True,
            capture_output=True,
            check=False,
        )

    def read_output(self) -> list[dict[str, str]]:
        with (self.tmp / "output.csv").open(newline="") as handle:
            return list(csv.DictReader(handle))

    def test_explicit_lane_column_remains_valid(self) -> None:
        r1 = self.fastq("sample_R1.fastq.gz")
        r2 = self.fastq("sample_R2.fastq.gz")
        result = self.run_validator(
            [
                {
                    "sample": "SampleA",
                    "lane": "L001",
                    "fastq_1": str(r1),
                    "fastq_2": str(r2),
                }
            ]
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        row = self.read_output()[0]
        self.assertEqual(row["lane"], "L001")
        self.assertEqual(row["control"], "")
        self.assertEqual(row["is_control"], "false")

    def test_lane_column_is_optional_and_rows_are_auto_numbered(self) -> None:
        result = self.run_validator(
            [
                {
                    "sample": "SampleA",
                    "fastq_1": str(self.fastq("SampleA_L001_R1.fastq.gz")),
                    "fastq_2": str(self.fastq("SampleA_L001_R2.fastq.gz")),
                },
                {
                    "sample": "SampleA",
                    "fastq_1": str(self.fastq("SampleA_L002_R1.fastq.gz")),
                    "fastq_2": str(self.fastq("SampleA_L002_R2.fastq.gz")),
                },
            ]
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(
            [row["lane"] for row in self.read_output()],
            ["L001", "L002"],
        )

    def test_blank_explicit_lane_is_rejected(self) -> None:
        result = self.run_validator(
            [
                {
                    "sample": "SampleA",
                    "lane": "",
                    "fastq_1": str(self.fastq("SampleA_R1.fastq.gz")),
                    "fastq_2": str(self.fastq("SampleA_R2.fastq.gz")),
                }
            ]
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("omit the lane column", result.stderr)

    def test_fastq_paths_with_spaces_are_preserved(self) -> None:
        directory = self.tmp / "directory with spaces"
        r1 = self.fastq("sample R1.fastq.gz", directory)
        r2 = self.fastq("sample R2.fastq.gz", directory)
        result = self.run_validator(
            [
                {
                    "sample": "SampleA",
                    "lane": "L001",
                    "fastq_1": str(r1),
                    "fastq_2": str(r2),
                }
            ]
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        row = self.read_output()[0]
        self.assertEqual(row["fastq_1"], str(r1.resolve()))

    def test_blank_control_is_allowed(self) -> None:
        result = self.run_validator(
            [
                {
                    "sample": "Uncontrolled",
                    "lane": "L001",
                    "fastq_1": str(self.fastq("Uncontrolled_R1.fastq.gz")),
                    "fastq_2": str(self.fastq("Uncontrolled_R2.fastq.gz")),
                    "control": "",
                }
            ],
            "ndigenome",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        row = self.read_output()[0]
        self.assertEqual(row["control"], "")
        self.assertEqual(row["is_control"], "false")

    def test_control_and_variant_metadata_are_normalized(self) -> None:
        vcf = self.tmp / "donor.vcf.gz"
        vcf.touch()
        Path(f"{vcf}.tbi").touch()
        rows = [
            {
                "sample": "Treated",
                "lane": "L001",
                "fastq_1": str(self.fastq("Treated_R1.fastq.gz")),
                "fastq_2": str(self.fastq("Treated_R2.fastq.gz")),
                "control": "Control",
                "variant_vcf": str(vcf),
            },
            {
                "sample": "Control",
                "lane": "L001",
                "fastq_1": str(self.fastq("Control_R1.fastq.gz")),
                "fastq_2": str(self.fastq("Control_R2.fastq.gz")),
                "control": "",
                "variant_vcf": "",
            },
        ]
        result = self.run_validator(rows, "ndigenome")
        self.assertEqual(result.returncode, 0, result.stderr)
        output = {row["sample"]: row for row in self.read_output()}
        self.assertEqual(output["Treated"]["control"], "Control")
        self.assertEqual(output["Treated"]["is_control"], "false")
        self.assertEqual(output["Control"]["control"], "")
        self.assertEqual(output["Control"]["is_control"], "true")
        self.assertTrue(output["Treated"]["variant_index"].endswith(".tbi"))

    def test_ndigenome_rejects_single_end_data(self) -> None:
        result = self.run_validator(
            [
                {
                    "sample": "SampleA",
                    "lane": "L001",
                    "fastq_1": str(self.fastq("sample.fastq.gz")),
                    "fastq_2": "",
                }
            ],
            "ndigenome",
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("requires paired-end", result.stderr)

    def test_duplicate_staged_basenames_are_rejected(self) -> None:
        first = self.fastq("reads_R1.fastq.gz", self.tmp / "lane1")
        second = self.fastq("reads_R1.fastq.gz", self.tmp / "lane2")
        result = self.run_validator(
            [
                {
                    "sample": "SampleA",
                    "lane": "L001",
                    "fastq_1": str(first),
                    "fastq_2": "",
                },
                {
                    "sample": "SampleA",
                    "lane": "L002",
                    "fastq_1": str(second),
                    "fastq_2": "",
                },
            ]
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("same staged basename", result.stderr)

    def test_one_control_can_be_shared_by_multiple_samples(self) -> None:
        rows = []
        for sample, control in (
            ("TreatedA", "Control"),
            ("TreatedB", "Control"),
            ("Control", ""),
        ):
            rows.append(
                {
                    "sample": sample,
                    "lane": "L001",
                    "fastq_1": str(self.fastq(f"{sample}_R1.fastq.gz")),
                    "fastq_2": str(self.fastq(f"{sample}_R2.fastq.gz")),
                    "control": control,
                    "variant_vcf": "",
                }
            )
        result = self.run_validator(rows, "ndigenome")
        self.assertEqual(result.returncode, 0, result.stderr)
        output = {row["sample"]: row for row in self.read_output()}
        self.assertEqual(output["TreatedA"]["control"], "Control")
        self.assertEqual(output["TreatedB"]["control"], "Control")
        self.assertEqual(output["Control"]["is_control"], "true")

    def test_nonexistent_named_control_is_rejected(self) -> None:
        result = self.run_validator(
            [
                {
                    "sample": "Treated",
                    "lane": "L001",
                    "fastq_1": str(self.fastq("Treated_R1.fastq.gz")),
                    "fastq_2": str(self.fastq("Treated_R2.fastq.gz")),
                    "control": "MissingControl",
                }
            ],
            "ndigenome",
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no sample named 'MissingControl'", result.stderr)

    def test_self_control_is_rejected(self) -> None:
        result = self.run_validator(
            [
                {
                    "sample": "SampleA",
                    "lane": "L001",
                    "fastq_1": str(self.fastq("SampleA_R1.fastq.gz")),
                    "fastq_2": str(self.fastq("SampleA_R2.fastq.gz")),
                    "control": "SampleA",
                }
            ],
            "ndigenome",
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("cannot use itself", result.stderr)

    def test_control_chain_is_rejected(self) -> None:
        rows = []
        for sample, control in (
            ("Treated", "ControlA"),
            ("ControlA", "ControlB"),
            ("ControlB", ""),
        ):
            rows.append(
                {
                    "sample": sample,
                    "lane": "L001",
                    "fastq_1": str(self.fastq(f"{sample}_R1.fastq.gz")),
                    "fastq_2": str(self.fastq(f"{sample}_R2.fastq.gz")),
                    "control": control,
                }
            )
        result = self.run_validator(rows, "ndigenome")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Control rows must leave the control column blank", result.stderr)

    def test_lane_metadata_must_be_consistent(self) -> None:
        rows = []
        for lane, control in (("L001", "ControlA"), ("L002", "ControlB")):
            rows.append(
                {
                    "sample": "SampleA",
                    "lane": lane,
                    "fastq_1": str(self.fastq(f"{lane}_R1.fastq.gz")),
                    "fastq_2": str(self.fastq(f"{lane}_R2.fastq.gz")),
                    "control": control,
                    "variant_vcf": "",
                }
            )
        result = self.run_validator(rows, "ndigenome")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("inconsistent control", result.stderr)

    def test_unknown_columns_are_rejected(self) -> None:
        result = self.run_validator(
            [
                {
                    "sample": "SampleA",
                    "lane": "L001",
                    "fastq_1": str(self.fastq("SampleA_R1.fastq.gz")),
                    "fastq_2": str(self.fastq("SampleA_R2.fastq.gz")),
                    "contol": "",
                }
            ]
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("unknown column(s): contol", result.stderr)

    def test_identical_r1_and_r2_are_rejected(self) -> None:
        shared = self.fastq("SampleA.fastq.gz")
        result = self.run_validator(
            [
                {
                    "sample": "SampleA",
                    "fastq_1": str(shared),
                    "fastq_2": str(shared),
                }
            ]
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("fastq_2 reuses FASTQ", result.stderr)

    def test_fastq_reuse_across_rows_is_rejected(self) -> None:
        shared = self.fastq("shared_R1.fastq.gz")
        result = self.run_validator(
            [
                {
                    "sample": "SampleA",
                    "fastq_1": str(shared),
                    "fastq_2": str(self.fastq("SampleA_R2.fastq.gz")),
                },
                {
                    "sample": "SampleB",
                    "fastq_1": str(shared),
                    "fastq_2": str(self.fastq("SampleB_R2.fastq.gz")),
                },
            ]
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("already used as fastq_1 for sample 'SampleA'", result.stderr)


if __name__ == "__main__":
    unittest.main()
