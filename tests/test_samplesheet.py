from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

from cleavage.samplesheet import validate_samplesheet, write_samples_json  # noqa: E402


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

    def test_valid_pe_sheet_produces_exact_record(self) -> None:
        """Test valid paired-end sheet produces expected record."""
        r1 = self.fastq("sample_R1.fastq.gz")
        r2 = self.fastq("sample_R2.fastq.gz")
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2\n"
            f"SampleA,{r1},{r2}\n"
        )
        records = validate_samplesheet(input_csv, "digenome")
        self.assertEqual(len(records), 1)
        record = records[0]
        self.assertEqual(record["sample"], "SampleA")
        self.assertFalse(record["single_end"])
        self.assertEqual(record["fastq_1"], [str(r1.resolve())])
        self.assertEqual(record["fastq_2"], [str(r2.resolve())])
        self.assertEqual(record["control"], "")
        self.assertEqual(record["variant_vcf"], "")
        self.assertEqual(record["variant_index"], "")
        self.assertFalse(record["is_control"])

    def test_fastq_paths_with_spaces_are_preserved(self) -> None:
        """Test that paths with spaces are correctly resolved."""
        directory = self.tmp / "directory with spaces"
        r1 = self.fastq("sample R1.fastq.gz", directory)
        r2 = self.fastq("sample R2.fastq.gz", directory)
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2\n"
            f"SampleA,{r1},{r2}\n"
        )
        records = validate_samplesheet(input_csv, "digenome")
        self.assertEqual(records[0]["fastq_1"], [str(r1.resolve())])
        self.assertEqual(records[0]["fastq_2"], [str(r2.resolve())])

    def test_blank_control_is_allowed(self) -> None:
        """Test that blank control is allowed."""
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2,control\n"
            f"Uncontrolled,{self.fastq('U_R1.fastq.gz')},{self.fastq('U_R2.fastq.gz')},\n"
        )
        records = validate_samplesheet(input_csv, "ndigenome")
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["control"], "")
        self.assertFalse(records[0]["is_control"])

    def test_control_and_variant_metadata_are_normalized(self) -> None:
        """Test VCF and index are normalized and is_control is set correctly."""
        vcf = self.tmp / "donor.vcf.gz"
        vcf.touch()
        Path(f"{vcf}.tbi").touch()
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2,control,variant_vcf\n"
            f"Treated,{self.fastq('T_R1.fastq.gz')},{self.fastq('T_R2.fastq.gz')},Control,{vcf}\n"
            f"Control,{self.fastq('C_R1.fastq.gz')},{self.fastq('C_R2.fastq.gz')},\n"
        )
        records = validate_samplesheet(input_csv, "ndigenome")
        records_by_sample = {r["sample"]: r for r in records}
        self.assertEqual(records_by_sample["Treated"]["control"], "Control")
        self.assertFalse(records_by_sample["Treated"]["is_control"])
        self.assertEqual(records_by_sample["Control"]["control"], "")
        self.assertTrue(records_by_sample["Control"]["is_control"])
        self.assertTrue(
            records_by_sample["Treated"]["variant_index"].endswith(".tbi")
        )

    def test_variant_index_csi_is_accepted(self) -> None:
        """Test that .csi index is accepted as well as .tbi."""
        vcf = self.tmp / "donor.vcf.gz"
        vcf.touch()
        Path(f"{vcf}.csi").touch()
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2,variant_vcf\n"
            f"Sample,{self.fastq('S_R1.fastq.gz')},{self.fastq('S_R2.fastq.gz')},{vcf}\n"
        )
        records = validate_samplesheet(input_csv, "digenome")
        self.assertTrue(
            records[0]["variant_index"].endswith(".csi")
        )

    def test_ndigenome_rejects_single_end_data(self) -> None:
        """Test that nDigenome rejects single-end data."""
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2\n"
            f"SampleA,{self.fastq('sample.fastq.gz')},\n"
        )
        with self.assertRaisesRegex(ValueError, "nDigenome requires paired-end"):
            validate_samplesheet(input_csv, "ndigenome")

    def test_digenome_accepts_single_end_data(self) -> None:
        """Test that digenome accepts single-end data with single_end: true."""
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2\n"
            f"SampleA,{self.fastq('sample.fastq.gz')},\n"
        )
        records = validate_samplesheet(input_csv, "digenome")
        self.assertEqual(len(records), 1)
        self.assertTrue(records[0]["single_end"])
        self.assertEqual(records[0]["fastq_2"], [])

    def test_two_lanes_of_one_sample_grouped_in_input_order(self) -> None:
        """Test that multiple lanes of one sample are grouped with input order preserved."""
        r1_l1 = self.fastq("SampleA_L001_R1.fastq.gz")
        r2_l1 = self.fastq("SampleA_L001_R2.fastq.gz")
        r1_l2 = self.fastq("SampleA_L002_R1.fastq.gz")
        r2_l2 = self.fastq("SampleA_L002_R2.fastq.gz")
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2\n"
            f"SampleA,{r1_l1},{r2_l1}\n"
            f"SampleA,{r1_l2},{r2_l2}\n"
        )
        records = validate_samplesheet(input_csv, "digenome")
        self.assertEqual(len(records), 1)
        self.assertEqual(len(records[0]["fastq_1"]), 2)
        self.assertEqual(len(records[0]["fastq_2"]), 2)
        self.assertEqual(records[0]["fastq_1"][0], str(r1_l1.resolve()))
        self.assertEqual(records[0]["fastq_1"][1], str(r1_l2.resolve()))

    def test_one_control_shared_by_multiple_samples(self) -> None:
        """Test that one control can be shared by multiple treated samples."""
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2,control\n"
            f"TreatedA,{self.fastq('TA_R1.fastq.gz')},{self.fastq('TA_R2.fastq.gz')},Control\n"
            f"TreatedB,{self.fastq('TB_R1.fastq.gz')},{self.fastq('TB_R2.fastq.gz')},Control\n"
            f"Control,{self.fastq('C_R1.fastq.gz')},{self.fastq('C_R2.fastq.gz')},\n"
        )
        records = validate_samplesheet(input_csv, "ndigenome")
        records_by_sample = {r["sample"]: r for r in records}
        self.assertEqual(records_by_sample["TreatedA"]["control"], "Control")
        self.assertEqual(records_by_sample["TreatedB"]["control"], "Control")
        self.assertTrue(records_by_sample["Control"]["is_control"])
        self.assertFalse(records_by_sample["TreatedA"]["is_control"])
        self.assertFalse(records_by_sample["TreatedB"]["is_control"])

    def test_nonexistent_named_control_is_rejected(self) -> None:
        """Test that referencing a non-existent control is rejected."""
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2,control\n"
            f"Treated,{self.fastq('T_R1.fastq.gz')},{self.fastq('T_R2.fastq.gz')},MissingControl\n"
        )
        with self.assertRaisesRegex(ValueError, "no sample named"):
            validate_samplesheet(input_csv, "ndigenome")

    def test_self_control_is_rejected(self) -> None:
        """Test that a sample using itself as control is rejected."""
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2,control\n"
            f"SampleA,{self.fastq('S_R1.fastq.gz')},{self.fastq('S_R2.fastq.gz')},SampleA\n"
        )
        with self.assertRaisesRegex(ValueError, "cannot use itself"):
            validate_samplesheet(input_csv, "ndigenome")

    def test_control_chain_is_rejected(self) -> None:
        """Test that control chains (control having its own control) are rejected."""
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2,control\n"
            f"Treated,{self.fastq('T_R1.fastq.gz')},{self.fastq('T_R2.fastq.gz')},ControlA\n"
            f"ControlA,{self.fastq('CA_R1.fastq.gz')},{self.fastq('CA_R2.fastq.gz')},ControlB\n"
            f"ControlB,{self.fastq('CB_R1.fastq.gz')},{self.fastq('CB_R2.fastq.gz')},\n"
        )
        with self.assertRaisesRegex(
            ValueError,
            "sample 'ControlA' is the control of 'Treated' but declares control 'ControlB'",
        ):
            validate_samplesheet(input_csv, "ndigenome")

    def test_inconsistent_metadata_across_lanes_rejected(self) -> None:
        """Test that inconsistent control or VCF across lanes is rejected."""
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2,control\n"
            f"SampleA,{self.fastq('L1_R1.fastq.gz')},{self.fastq('L1_R2.fastq.gz')},ControlA\n"
            f"SampleA,{self.fastq('L2_R1.fastq.gz')},{self.fastq('L2_R2.fastq.gz')},ControlB\n"
        )
        with self.assertRaisesRegex(ValueError, "inconsistent control"):
            validate_samplesheet(input_csv, "ndigenome")

    def test_unknown_columns_rejected_including_lane(self) -> None:
        """Test that unknown columns including 'lane' are rejected."""
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2,lane\n"
            f"SampleA,{self.fastq('S_R1.fastq.gz')},{self.fastq('S_R2.fastq.gz')},L001\n"
        )
        with self.assertRaisesRegex(ValueError, "unknown column"):
            validate_samplesheet(input_csv, "digenome")

    def test_unknown_column_typo_rejected(self) -> None:
        """Test that typos in column names are rejected."""
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2,contol\n"
            f"SampleA,{self.fastq('S_R1.fastq.gz')},{self.fastq('S_R2.fastq.gz')},\n"
        )
        with self.assertRaisesRegex(ValueError, "unknown column"):
            validate_samplesheet(input_csv, "digenome")

    def test_duplicate_column_names_rejected(self) -> None:
        """Test that duplicate column names are rejected."""
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2,fastq_1\n"
            f"SampleA,{self.fastq('S_R1.fastq.gz')},{self.fastq('S_R2.fastq.gz')},extra\n"
        )
        with self.assertRaisesRegex(ValueError, "duplicate column"):
            validate_samplesheet(input_csv, "digenome")

    def test_identical_r1_and_r2_rejected(self) -> None:
        """Test that identical fastq_1 and fastq_2 are rejected."""
        shared = self.fastq("SampleA.fastq.gz")
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2\n"
            f"SampleA,{shared},{shared}\n"
        )
        with self.assertRaisesRegex(ValueError, "fastq_2 reuses"):
            validate_samplesheet(input_csv, "digenome")

    def test_fastq_reuse_across_rows_rejected(self) -> None:
        """Test that a FASTQ file reused across different samples is rejected."""
        shared = self.fastq("shared_R1.fastq.gz")
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2\n"
            f"SampleA,{shared},{self.fastq('A_R2.fastq.gz')}\n"
            f"SampleB,{shared},{self.fastq('B_R2.fastq.gz')}\n"
        )
        with self.assertRaisesRegex(ValueError, r"line 3: fastq_1 reuses FASTQ .* already used on line 2 \(SampleA\) as fastq_1"):
            validate_samplesheet(input_csv, "digenome")

    def test_same_basename_different_directories_accepted(self) -> None:
        """Test that same FASTQ basename in different directories is now accepted."""
        dir1 = self.tmp / "lane1"
        dir2 = self.tmp / "lane2"
        r1_d1 = self.fastq("reads_R1.fastq.gz", dir1)
        r2_d1 = self.fastq("reads_R2.fastq.gz", dir1)
        r1_d2 = self.fastq("reads_R1.fastq.gz", dir2)
        r2_d2 = self.fastq("reads_R2.fastq.gz", dir2)
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2\n"
            f"SampleA,{r1_d1},{r2_d1}\n"
            f"SampleA,{r1_d2},{r2_d2}\n"
        )
        records = validate_samplesheet(input_csv, "digenome")
        self.assertEqual(len(records), 1)
        self.assertEqual(len(records[0]["fastq_1"]), 2)

    def test_missing_fastq_file_rejected(self) -> None:
        """Test that missing FASTQ file is rejected."""
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2\n"
            f"SampleA,/nonexistent/file.fastq.gz,{self.fastq('S_R2.fastq.gz')}\n"
        )
        with self.assertRaisesRegex(ValueError, "file does not exist"):
            validate_samplesheet(input_csv, "digenome")

    def test_wrong_fastq_suffix_rejected(self) -> None:
        """Test that FASTQ with wrong suffix is rejected."""
        fastq = self.fastq("sample.fastq")
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2\n"
            f"SampleA,{fastq},{self.fastq('S_R2.fastq.gz')}\n"
        )
        with self.assertRaisesRegex(ValueError, "must end with"):
            validate_samplesheet(input_csv, "digenome")

    def test_unindexed_vcf_rejected(self) -> None:
        """Test that VCF without index is rejected."""
        vcf = self.tmp / "donor.vcf.gz"
        vcf.touch()
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2,variant_vcf\n"
            f"Sample,{self.fastq('S_R1.fastq.gz')},{self.fastq('S_R2.fastq.gz')},{vcf}\n"
        )
        with self.assertRaisesRegex(ValueError, "not indexed"):
            validate_samplesheet(input_csv, "digenome")

    def test_wrong_vcf_suffix_rejected(self) -> None:
        """Test that VCF with wrong suffix is rejected."""
        vcf = self.tmp / "donor.vcf"
        vcf.touch()
        Path(f"{vcf}.tbi").touch()
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2,variant_vcf\n"
            f"Sample,{self.fastq('S_R1.fastq.gz')},{self.fastq('S_R2.fastq.gz')},{vcf}\n"
        )
        with self.assertRaisesRegex(ValueError, "bgzip-compressed"):
            validate_samplesheet(input_csv, "digenome")

    def test_multiple_problems_all_reported(self) -> None:
        """Test that all problems in one sheet are reported together."""
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2\n"
            f"SampleA,{self.fastq('S_R1.fastq.gz')},\n"
            f"SampleA,/missing.fastq.gz,{self.fastq('S_R2.fastq.gz')}\n"
        )
        try:
            validate_samplesheet(input_csv, "ndigenome")
            self.fail("Should have raised ValueError")
        except ValueError as e:
            error_msg = str(e)
            self.assertIn("nDigenome requires", error_msg)
            self.assertIn("file does not exist", error_msg)

    def test_write_samples_json_round_trips(self) -> None:
        """Test that write_samples_json produces valid JSON that round-trips."""
        r1 = self.fastq("sample_R1.fastq.gz")
        r2 = self.fastq("sample_R2.fastq.gz")
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2\n"
            f"SampleA,{r1},{r2}\n"
        )
        records = validate_samplesheet(input_csv, "digenome")
        output_json = self.tmp / "output.json"
        write_samples_json(records, output_json)

        # Verify JSON is valid and round-trips
        with output_json.open() as f:
            content = f.read()
        self.assertTrue(content.endswith("\n"))
        loaded = json.loads(content)
        self.assertEqual(loaded, records)

    def test_blank_sample_skipped(self) -> None:
        """Test that rows with blank sample are skipped."""
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2\n"
            f",{self.fastq('B1_R1.fastq.gz')},{self.fastq('B1_R2.fastq.gz')}\n"
            f"SampleA,{self.fastq('S_R1.fastq.gz')},{self.fastq('S_R2.fastq.gz')}\n"
        )
        with self.assertRaisesRegex(ValueError, "sample is blank"):
            validate_samplesheet(input_csv, "digenome")

    def test_blank_fastq_1_rejected(self) -> None:
        """Test that blank fastq_1 is rejected."""
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2\n"
            f"SampleA,,{self.fastq('S_R2.fastq.gz')}\n"
        )
        with self.assertRaisesRegex(ValueError, "fastq_1 is blank"):
            validate_samplesheet(input_csv, "digenome")

    def test_valid_long_read_sheet_with_treated_and_control(self) -> None:
        """Test valid long-read sheet with treated sample and control."""
        treated_bam = self.tmp / "treated.bam"
        treated_bam.touch()
        Path(f"{treated_bam}.bai").touch()
        control_bam = self.tmp / "control.bam"
        control_bam.touch()
        Path(f"{control_bam}.bai").touch()
        vcf = self.tmp / "donor.vcf.gz"
        vcf.touch()
        Path(f"{vcf}.tbi").touch()

        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,bam,control,variant_vcf\n"
            f"Treated,{treated_bam},Control,{vcf}\n"
            f"Control,{control_bam},,\n"
        )
        records = validate_samplesheet(input_csv, "digenome", long_reads=True)
        self.assertEqual(len(records), 2)

        treated = records[0]
        self.assertEqual(treated["sample"], "Control")
        self.assertEqual(treated["bam"], str(control_bam.resolve()))
        self.assertEqual(treated["bam_index"], str(control_bam.resolve()) + ".bai")
        self.assertEqual(treated["control"], "")
        self.assertEqual(treated["variant_vcf"], "")
        self.assertEqual(treated["variant_index"], "")
        self.assertTrue(treated["is_control"])

        control = records[1]
        self.assertEqual(control["sample"], "Treated")
        self.assertEqual(control["bam"], str(treated_bam.resolve()))
        self.assertEqual(control["bam_index"], str(treated_bam.resolve()) + ".bai")
        self.assertEqual(control["control"], "Control")
        self.assertEqual(control["variant_vcf"], str(vcf.resolve()))
        self.assertEqual(control["variant_index"], str(vcf.resolve()) + ".tbi")
        self.assertFalse(control["is_control"])

    def test_long_read_sheet_with_csi_index(self) -> None:
        """Test that .csi index is accepted for long-read BAM."""
        bam = self.tmp / "sample.bam"
        bam.touch()
        Path(f"{bam}.csi").touch()

        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,bam\n"
            f"Sample,{bam}\n"
        )
        records = validate_samplesheet(input_csv, "digenome", long_reads=True)
        self.assertEqual(len(records), 1)
        self.assertTrue(records[0]["bam_index"].endswith(".csi"))

    def test_long_read_missing_bam_index_rejected(self) -> None:
        """Test that BAM without .bai or .csi index is rejected."""
        bam = self.tmp / "sample.bam"
        bam.touch()

        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,bam\n"
            f"Sample,{bam}\n"
        )
        with self.assertRaisesRegex(ValueError, "not indexed with .bai or .csi"):
            validate_samplesheet(input_csv, "digenome", long_reads=True)

    def test_long_read_non_bam_suffix_rejected(self) -> None:
        """Test that non-.bam path is rejected."""
        file_path = self.tmp / "sample.txt"
        file_path.touch()

        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,bam\n"
            f"Sample,{file_path}\n"
        )
        with self.assertRaisesRegex(ValueError, "must end with .bam"):
            validate_samplesheet(input_csv, "digenome", long_reads=True)

    def test_long_read_missing_bam_file_rejected(self) -> None:
        """Test that missing BAM file is rejected."""
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,bam\n"
            "Sample,/nonexistent/sample.bam\n"
        )
        with self.assertRaisesRegex(ValueError, "file does not exist"):
            validate_samplesheet(input_csv, "digenome", long_reads=True)

    def test_long_read_blank_bam_rejected(self) -> None:
        """Test that blank bam is rejected."""
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,bam\n"
            "Sample,\n"
        )
        with self.assertRaisesRegex(ValueError, "bam is blank"):
            validate_samplesheet(input_csv, "digenome", long_reads=True)

    def test_long_read_bam_reuse_rejected(self) -> None:
        """Test that the same BAM on two rows is rejected."""
        bam = self.tmp / "shared.bam"
        bam.touch()
        Path(f"{bam}.bai").touch()

        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,bam\n"
            f"Sample1,{bam}\n"
            f"Sample2,{bam}\n"
        )
        with self.assertRaisesRegex(ValueError, "reuses BAM"):
            validate_samplesheet(input_csv, "digenome", long_reads=True)

    def test_long_read_one_sample_one_row_enforced(self) -> None:
        """Test that one sample can have only one row in long-read mode."""
        bam1 = self.tmp / "sample1.bam"
        bam1.touch()
        Path(f"{bam1}.bai").touch()
        bam2 = self.tmp / "sample2.bam"
        bam2.touch()
        Path(f"{bam2}.bai").touch()

        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,bam\n"
            f"Sample,{bam1}\n"
            f"Sample,{bam2}\n"
        )
        with self.assertRaisesRegex(ValueError, "has more than one row"):
            validate_samplesheet(input_csv, "digenome", long_reads=True)

    def test_short_read_rejects_fastq_columns_in_long_read_mode(self) -> None:
        """Test that fastq_1/fastq_2 columns are rejected with --long_reads."""
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,fastq_1,fastq_2\n"
            f"Sample,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz\n"
        )
        with self.assertRaisesRegex(
            ValueError,
            "With --long_reads the samplesheet lists aligned BAMs in a bam column, not fastq_1/fastq_2"
        ):
            validate_samplesheet(input_csv, "digenome", long_reads=True)

    def test_short_read_rejects_bam_column_without_long_reads(self) -> None:
        """Test that bam column is rejected without --long_reads."""
        bam = self.tmp / "sample.bam"
        bam.touch()
        Path(f"{bam}.bai").touch()

        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,bam\n"
            f"Sample,{bam}\n"
        )
        with self.assertRaisesRegex(ValueError, "The bam column needs --long_reads"):
            validate_samplesheet(input_csv, "digenome", long_reads=False)

    def test_long_read_shared_control_allowed(self) -> None:
        """Test that one control can be shared by multiple treated samples."""
        treated1_bam = self.tmp / "treated1.bam"
        treated1_bam.touch()
        Path(f"{treated1_bam}.bai").touch()
        treated2_bam = self.tmp / "treated2.bam"
        treated2_bam.touch()
        Path(f"{treated2_bam}.bai").touch()
        control_bam = self.tmp / "control.bam"
        control_bam.touch()
        Path(f"{control_bam}.bai").touch()

        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,bam,control\n"
            f"Control,{control_bam},\n"
            f"Treated1,{treated1_bam},Control\n"
            f"Treated2,{treated2_bam},Control\n"
        )
        records = validate_samplesheet(input_csv, "digenome", long_reads=True)
        records_by_sample = {r["sample"]: r for r in records}
        self.assertEqual(records_by_sample["Treated1"]["control"], "Control")
        self.assertEqual(records_by_sample["Treated2"]["control"], "Control")
        self.assertTrue(records_by_sample["Control"]["is_control"])
        self.assertFalse(records_by_sample["Treated1"]["is_control"])
        self.assertFalse(records_by_sample["Treated2"]["is_control"])

    def test_long_read_control_chain_rejected(self) -> None:
        """Test that control chains are rejected in long-read mode."""
        treated_bam = self.tmp / "treated.bam"
        treated_bam.touch()
        Path(f"{treated_bam}.bai").touch()
        control_a_bam = self.tmp / "control_a.bam"
        control_a_bam.touch()
        Path(f"{control_a_bam}.bai").touch()
        control_b_bam = self.tmp / "control_b.bam"
        control_b_bam.touch()
        Path(f"{control_b_bam}.bai").touch()

        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,bam,control\n"
            f"Treated,{treated_bam},ControlA\n"
            f"ControlA,{control_a_bam},ControlB\n"
            f"ControlB,{control_b_bam},\n"
        )
        with self.assertRaisesRegex(ValueError, "but declares control"):
            validate_samplesheet(input_csv, "digenome", long_reads=True)

    def test_long_read_ndigenome_allowed(self) -> None:
        """Test that ndigenome analysis is allowed with long reads."""
        bam = self.tmp / "sample.bam"
        bam.touch()
        Path(f"{bam}.bai").touch()

        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,bam\n"
            f"Sample,{bam}\n"
        )
        records = validate_samplesheet(input_csv, "ndigenome", long_reads=True)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["sample"], "Sample")

    def test_long_read_multiple_problems_reported(self) -> None:
        """Test that multiple problems are reported together."""
        input_csv = self.tmp / "input.csv"
        input_csv.write_text(
            "sample,bam\n"
            "Sample1,/missing1.bam\n"
            "Sample2,/missing2.bam\n"
        )
        try:
            validate_samplesheet(input_csv, "digenome", long_reads=True)
            self.fail("Should have raised ValueError")
        except ValueError as e:
            error_msg = str(e)
            self.assertIn("missing1.bam", error_msg)
            self.assertIn("missing2.bam", error_msg)


if __name__ == "__main__":
    unittest.main()
