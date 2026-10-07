"""Tests for finalization: merging chunks, filtering, and validation."""

from __future__ import annotations

import csv
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

try:
    import pysam
except ImportError:
    pysam = None

sys.path.insert(0, str(Path(__file__).resolve().parent))
from helpers import (
    make_read,
    reverse_read_ending_at,
    forward_background,
    reverse_background,
    write_bam,
    write_vcf,
    run_caller,
    make_settings,
    CHUNK_SETTINGS,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from cleavage.call import call_chunk
from cleavage.finalize import finalize
from cleavage.output import format_value
from cleavage.stats import benjamini_hochberg


@unittest.skipIf(pysam is None, "pysam is not installed")
class FinalizeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def endpoint_site_reads(
        self,
        contig_id: int,
        position: int,
        signal_count: int,
        prefix: str,
    ) -> list:
        """Reads for a single Digenome site: forward signal, reverse signal, background."""
        reads = [
            make_read(f"{prefix}_fwd_sig_{i}", position, contig=contig_id, mapq=0)
            for i in range(signal_count)
        ]
        reads += [
            make_read(f"{prefix}_rev_sig_{i}", position - 49, contig=contig_id, reverse=True, mapq=0)
            for i in range(signal_count)
        ]
        reads += forward_background(f"{prefix}_fwd_bg", position, 12, contig=contig_id, mapq=0)
        reads += [
            make_read(f"{prefix}_rev_bg_{i}", position - 48 + i, contig=contig_id, reverse=True, mapq=0)
            for i in range(12)
        ]
        return reads

    def control_site_reads(
        self,
        contig_id: int,
        position: int,
        prefix: str,
    ) -> list:
        """Control reads at a single position."""
        reads = [
            make_read(f"{prefix}_fwd_{i}", position - 20 - i, contig=contig_id, mapq=0)
            for i in range(20)
        ]
        reads += [
            make_read(f"{prefix}_rev_{i}", position - 48 + i, contig=contig_id, reverse=True, mapq=0)
            for i in range(20)
        ]
        return reads

    def test_serial_and_chunked_outputs_are_equivalent(self) -> None:
        """Multi-contig with control, VCF, blacklist: chunks=1 and chunks=2,3,4,6 are identical."""
        # Build multi-contig treated and control BAMs
        treated_reads = self.endpoint_site_reads(0, 500, 8, "chr1")
        treated_reads += self.endpoint_site_reads(1, 700, 6, "chr2")
        treated_reads += [
            make_read(f"chr3_bg_{i}", 100 + i * 3, contig=2, mapq=0)
            for i in range(12)
        ]
        control_reads = self.control_site_reads(0, 500, "chr1_ctrl")
        control_reads += self.control_site_reads(1, 700, "chr2_ctrl")
        control_reads += [
            make_read(f"chr3_ctrl_{i}", 200 + i * 3, contig=2, mapq=0)
            for i in range(12)
        ]

        treated_bam = write_bam(
            self.tmp / "treated.bam",
            [("chr1", 2000), ("chr2", 2000), ("chr3", 2000)],
            treated_reads,
        )
        control_bam = write_bam(
            self.tmp / "control.bam",
            [("chr1", 2000), ("chr2", 2000), ("chr3", 2000)],
            control_reads,
        )
        vcf = write_vcf(
            self.tmp / "variants.vcf.gz",
            [("chr1", 2000), ("chr2", 2000), ("chr3", 2000)],
            ["chr2\t701\t.\tAA\tA\t60\tPASS\t."],
        )
        blacklist = self.tmp / "blacklist.bed"
        blacklist.write_text("chr1\t490\t510\n")

        settings = make_settings("digenome", base=CHUNK_SETTINGS)

        # Serial run with chunks=1
        serial_prefix = self.tmp / "serial"
        serial_prefix, serial_qc = run_caller(
            settings,
            {
                "bam": treated_bam,
                "control_bam": control_bam,
                "vcf": vcf,
                "blacklist": blacklist,
            },
            self.tmp / "serial_out",
            chunks=1,
            sample="Sample",
            control_sample="Control",
        )

        # Chunked runs with various chunk counts
        for chunk_count in [2, 3, 4, 6]:
            with self.subTest(chunks=chunk_count):
                chunked_prefix = self.tmp / f"chunked_{chunk_count}"
                chunked_prefix, chunked_qc = run_caller(
                    settings,
                    {
                        "bam": treated_bam,
                        "control_bam": control_bam,
                        "vcf": vcf,
                        "blacklist": blacklist,
                    },
                    self.tmp / f"chunked_{chunk_count}_out",
                    chunks=chunk_count,
                    sample="Sample",
                    control_sample="Control",
                )

                # Check all output files are byte-identical
                for suffix in [
                    ".digenome.all.tsv",
                    ".digenome.high_confidence.tsv",
                    ".digenome.manual_review.tsv",
                    ".digenome.artifact.tsv",
                    ".digenome.bed",
                    ".digenome_mqc.tsv",
                ]:
                    with self.subTest(suffix=suffix):
                        serial_path = Path(f"{serial_prefix}{suffix}")
                        chunked_path = Path(f"{chunked_prefix}{suffix}")
                        if serial_path.exists() and chunked_path.exists():
                            self.assertEqual(
                                serial_path.read_bytes(),
                                chunked_path.read_bytes(),
                                f"Mismatch in {suffix}",
                            )

                # Check QC stats match
                self.assertEqual(serial_qc["rows"], chunked_qc["rows"])
                self.assertEqual(serial_qc["tiers"], chunked_qc["tiers"])
                self.assertEqual(
                    serial_qc["candidates_before_filters"],
                    chunked_qc["candidates_before_filters"],
                )

                # Verify q-values match benjamini_hochberg
                all_tsv_path = Path(f"{chunked_prefix}.digenome.all.tsv")
                with all_tsv_path.open(newline="") as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                p_values = [
                    float(row["control_fisher_p"])
                    for row in rows
                    if row["control_fisher_p"]
                ]
                if p_values:
                    q_values = [
                        float(row["control_fisher_q"])
                        for row in rows
                        if row["control_fisher_q"]
                    ]
                    expected_q = benjamini_hochberg(p_values)
                    for observed, wanted in zip(q_values, expected_q):
                        self.assertAlmostEqual(observed, wanted)

    def test_conflict_chain_across_chunk_boundary(self) -> None:
        """Chain 980..1020 step 2 with 5 reads each: chunks=2 == chunks=1."""
        positions = list(range(980, 1021, 2))
        reads = []
        for position in positions:
            reads.extend(
                make_read(f"chain_fwd_{position}_{i}", position, mapq=0)
                for i in range(5)
            )
            reads.extend(
                make_read(f"chain_rev_{position}_{i}", position - 49, reverse=True, mapq=0)
                for i in range(5)
            )

        bam = write_bam(self.tmp / "treated.bam", [("chr1", 2000)], reads)
        vcf = write_vcf(self.tmp / "variants.vcf.gz", [("chr1", 2000)], [])

        settings = make_settings(
            "digenome",
            base=CHUNK_SETTINGS,
            digenome_overhang=0,
            digenome_pair_window=2,
        )

        # Serial run
        serial_prefix, serial_qc = run_caller(
            settings,
            {"bam": bam, "control_bam": None, "vcf": vcf, "blacklist": None},
            self.tmp / "serial_out",
            chunks=1,
            sample="Sample",
        )

        # Chunked run with 2 chunks
        chunked_prefix, chunked_qc = run_caller(
            settings,
            {"bam": bam, "control_bam": None, "vcf": vcf, "blacklist": None},
            self.tmp / "chunked_out",
            chunks=2,
            sample="Sample",
        )

        # Outputs should be byte-identical
        self.assertEqual(
            Path(f"{serial_prefix}.digenome.all.tsv").read_bytes(),
            Path(f"{chunked_prefix}.digenome.all.tsv").read_bytes(),
        )

        # Candidate counts should match
        self.assertEqual(
            serial_qc["candidates_before_filters"],
            chunked_qc["candidates_before_filters"],
        )

        # Each of the 21 forward endpoints can pair with reverse endpoints at
        # -2, 0, and +2, except at the two ends of the chain.
        self.assertEqual(serial_qc["candidates_before_filters"], 21 * 3 - 2)

    def test_boundary_pair_with_overhang(self) -> None:
        """Forward 1001, reverse 999: chunks=2 == chunks=1, no duplicates."""
        forward_pos = 1001
        reverse_pos = 999
        reads = [
            make_read(f"fwd_sig_{i}", forward_pos, mapq=0)
            for i in range(8)
        ]
        reads += [
            make_read(f"rev_sig_{i}", reverse_pos - 49, reverse=True, mapq=0)
            for i in range(8)
        ]
        reads += forward_background("fwd_bg", forward_pos, 12, mapq=0)
        reads += [
            make_read(f"rev_bg_{i}", reverse_pos - 48 + i, reverse=True, mapq=0)
            for i in range(12)
        ]

        bam = write_bam(self.tmp / "treated.bam", [("chr1", 2000)], reads)
        vcf = write_vcf(self.tmp / "variants.vcf.gz", [("chr1", 2000)], [])

        settings = make_settings(
            "digenome",
            base=CHUNK_SETTINGS,
            digenome_overhang=2,
        )

        # Serial run
        serial_prefix, _ = run_caller(
            settings,
            {"bam": bam, "control_bam": None, "vcf": vcf, "blacklist": None},
            self.tmp / "serial_out",
            chunks=1,
            sample="Sample",
        )

        # Chunked run
        chunked_prefix, _ = run_caller(
            settings,
            {"bam": bam, "control_bam": None, "vcf": vcf, "blacklist": None},
            self.tmp / "chunked_out",
            chunks=2,
            sample="Sample",
        )

        # Outputs should be byte-identical
        self.assertEqual(
            Path(f"{serial_prefix}.digenome.all.tsv").read_bytes(),
            Path(f"{chunked_prefix}.digenome.all.tsv").read_bytes(),
        )

        # Check for duplicates
        all_tsv = Path(f"{chunked_prefix}.digenome.all.tsv")
        with all_tsv.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        keys = [
            (row["contig"], row["forward_position_0based"], row["reverse_position_0based"])
            for row in rows
        ]
        self.assertEqual(len(keys), len(set(keys)))

    def test_chunk_files_may_arrive_in_any_order(self) -> None:
        """Nextflow collects chunks in completion order; the finalizer must
        pair each summary with its own records regardless."""
        sys.path.insert(0, str(Path(__file__).resolve().parent / "golden"))
        from scenarios import multi_contig

        scenario = multi_contig()
        inputs = scenario.write_inputs(self.tmp / "inputs")
        for analysis in ("digenome", "ndigenome"):
            settings = make_settings(analysis, base=scenario.settings)
            summaries = []
            for index in range(3):
                prefix = self.tmp / f"{analysis}_chunk_{index}"
                call_chunk(
                    settings,
                    bam_path=str(inputs["bam"]),
                    sample="Sample",
                    index=index,
                    count=3,
                    out_prefix=str(prefix),
                    control_path=str(inputs["control_bam"]),
                    control_sample="Control",
                    vcf_path=str(inputs["vcf"]),
                    blacklist_path=str(inputs["blacklist"]),
                )
                summaries.append(f"{prefix}.json")
            in_order = str(self.tmp / f"{analysis}_in_order")
            reversed_order = str(self.tmp / f"{analysis}_reversed")
            finalize(settings, summaries, "Sample", "Control", in_order)
            finalize(settings, list(reversed(summaries)), "Sample", "Control", reversed_order)
            for suffix in ("all.tsv", "bed", "qc.json"):
                self.assertEqual(
                    Path(f"{in_order}.{analysis}.{suffix}").read_bytes(),
                    Path(f"{reversed_order}.{analysis}.{suffix}").read_bytes(),
                )

    def test_finalizer_missing_chunks(self) -> None:
        """Only chunk 1 of 2 supplied → error message contains expected text."""
        treated_reads = self.endpoint_site_reads(0, 500, 8, "chr1")
        treated_bam = write_bam(
            self.tmp / "treated.bam",
            [("chr1", 2000)],
            treated_reads,
        )

        settings = make_settings("digenome", base=CHUNK_SETTINGS)

        # Create chunks but only supply one
        summaries = []
        for index in range(2):
            chunk_prefix = self.tmp / f"chunk_{index:03d}"
            call_chunk(
                settings,
                bam_path=str(treated_bam),
                sample="Sample",
                index=index,
                count=2,
                out_prefix=str(chunk_prefix),
                control_path=None,
                control_sample="",
                vcf_path=None,
                blacklist_path=None,
            )
            if index == 0:  # Only add first chunk
                summaries.append(str(chunk_prefix) + ".json")

        # Try to finalize with incomplete summaries
        with self.assertRaisesRegex(ValueError, "Chunk .* expects"):
            finalize(settings, summaries, "Sample", "", str(self.tmp / "out"))

    def test_finalizer_wrong_sample(self) -> None:
        """Chunk from different sample → error mentions sample name."""
        treated_reads = self.endpoint_site_reads(0, 500, 8, "chr1")
        treated_bam = write_bam(
            self.tmp / "treated.bam",
            [("chr1", 2000)],
            treated_reads,
        )

        settings = make_settings("digenome", base=CHUNK_SETTINGS)

        # Create one chunk with sample name "Other"
        chunk_prefix = self.tmp / "chunk_000"
        call_chunk(
            settings,
            bam_path=str(treated_bam),
            sample="Other",
            index=0,
            count=1,
            out_prefix=str(chunk_prefix),
            control_path=None,
            control_sample="",
            vcf_path=None,
            blacklist_path=None,
        )

        # Try to finalize with different sample
        with self.assertRaisesRegex(ValueError, "belongs to sample"):
            finalize(
                settings,
                [str(chunk_prefix) + ".json"],
                "Sample",
                "",
                str(self.tmp / "out"),
            )

    def test_finalizer_different_settings(self) -> None:
        """Chunk created with different settings → error about settings mismatch."""
        treated_reads = self.endpoint_site_reads(0, 500, 8, "chr1")
        treated_bam = write_bam(
            self.tmp / "treated.bam",
            [("chr1", 2000)],
            treated_reads,
        )

        settings1 = make_settings("digenome", base=CHUNK_SETTINGS)

        # Create chunk with settings1
        chunk_prefix = self.tmp / "chunk_000"
        call_chunk(
            settings1,
            bam_path=str(treated_bam),
            sample="Sample",
            index=0,
            count=1,
            out_prefix=str(chunk_prefix),
            control_path=None,
            control_sample="",
            vcf_path=None,
            blacklist_path=None,
        )

        # Try to finalize with different settings
        settings2 = make_settings("digenome", base=CHUNK_SETTINGS, digenome_pair_score_cutoff=5.0)
        with self.assertRaisesRegex(ValueError, "different settings"):
            finalize(
                settings2,
                [str(chunk_prefix) + ".json"],
                "Sample",
                "",
                str(self.tmp / "out"),
            )

    def test_finalizer_record_outside_ownership(self) -> None:
        """Record moved outside chunk's ownership → error about ranges."""
        # Create reads on single contig at positions that will be split across chunks
        treated_reads = self.endpoint_site_reads(0, 400, 8, "chr1_pos1")
        treated_reads += self.endpoint_site_reads(0, 1600, 8, "chr1_pos2")
        treated_bam = write_bam(
            self.tmp / "treated.bam",
            [("chr1", 2000)],
            treated_reads,
        )

        settings = make_settings("digenome", base=CHUNK_SETTINGS)

        # Create 2 chunks; chunk 0 should own [0, 1000), chunk 1 should own [1000, 2000)
        chunk_prefixes = []
        for index in range(2):
            chunk_prefix = self.tmp / f"chunk_{index:03d}"
            call_chunk(
                settings,
                bam_path=str(treated_bam),
                sample="Sample",
                index=index,
                count=2,
                out_prefix=str(chunk_prefix),
                control_path=None,
                control_sample="",
                vcf_path=None,
                blacklist_path=None,
            )
            chunk_prefixes.append(chunk_prefix)

        # Corrupt chunk 0: change a position within its range to outside its range
        chunk_0_jsonl = f"{chunk_prefixes[0]}.jsonl.gz"
        with gzip.open(chunk_0_jsonl, "rt") as handle:
            lines = handle.readlines()
        if lines:
            # Move the first record outside chunk 0's range [0, 1000) to chunk 1's range [1000, 2000)
            record = json.loads(lines[0])
            original_pos = record["position_0based"]
            if original_pos < 1000:
                record["position_0based"] = 1500  # Move to chunk 1's range
            lines[0] = json.dumps(record) + "\n"

            with gzip.open(chunk_0_jsonl, "wt", compresslevel=1) as handle:
                handle.writelines(lines)

        # Try to finalize with corrupted chunk 0
        summaries = [str(cp) + ".json" for cp in chunk_prefixes]
        with self.assertRaisesRegex(ValueError, "outside the ranges it owns"):
            finalize(
                settings,
                summaries,
                "Sample",
                "",
                str(self.tmp / "out"),
            )

    def test_finalizer_missing_record(self) -> None:
        """One record line deleted from chunk → error about record count mismatch."""
        treated_reads = self.endpoint_site_reads(0, 500, 8, "chr1")
        treated_bam = write_bam(
            self.tmp / "treated.bam",
            [("chr1", 2000)],
            treated_reads,
        )

        settings = make_settings("digenome", base=CHUNK_SETTINGS)

        # Create chunk
        chunk_prefix = self.tmp / "chunk_000"
        call_chunk(
            settings,
            bam_path=str(treated_bam),
            sample="Sample",
            index=0,
            count=1,
            out_prefix=str(chunk_prefix),
            control_path=None,
            control_sample="",
            vcf_path=None,
            blacklist_path=None,
        )

        # Read the summary to get the expected record count
        summary = json.loads(Path(f"{chunk_prefix}.json").read_text())
        expected_records = summary["records"]

        # Remove one line from the jsonl.gz
        with gzip.open(f"{chunk_prefix}.jsonl.gz", "rt") as handle:
            lines = handle.readlines()
        if len(lines) > 1:
            lines = lines[:-1]  # Remove last line

            with gzip.open(f"{chunk_prefix}.jsonl.gz", "wt", compresslevel=1) as handle:
                handle.writelines(lines)

            # Try to finalize
            with self.assertRaisesRegex(
                ValueError, f"has {len(lines)} records, its summary says {expected_records}"
            ):
                finalize(
                    settings,
                    [str(chunk_prefix) + ".json"],
                    "Sample",
                    "",
                    str(self.tmp / "out"),
                )

    def test_finalizer_duplicate_record(self) -> None:
        """Duplicate record line → error about duplicate call."""
        treated_reads = self.endpoint_site_reads(0, 500, 8, "chr1")
        treated_bam = write_bam(
            self.tmp / "treated.bam",
            [("chr1", 2000)],
            treated_reads,
        )

        settings = make_settings("digenome", base=CHUNK_SETTINGS)

        # Create chunk
        chunk_prefix = self.tmp / "chunk_000"
        call_chunk(
            settings,
            bam_path=str(treated_bam),
            sample="Sample",
            index=0,
            count=1,
            out_prefix=str(chunk_prefix),
            control_path=None,
            control_sample="",
            vcf_path=None,
            blacklist_path=None,
        )

        # Read and duplicate one line
        with gzip.open(f"{chunk_prefix}.jsonl.gz", "rt") as handle:
            lines = handle.readlines()
        if lines:
            lines.insert(1, lines[0])  # Duplicate first line

            with gzip.open(f"{chunk_prefix}.jsonl.gz", "wt", compresslevel=1) as handle:
                handle.writelines(lines)

            # Try to finalize
            with self.assertRaisesRegex(ValueError, "Duplicate call"):
                finalize(
                    settings,
                    [str(chunk_prefix) + ".json"],
                    "Sample",
                    "",
                    str(self.tmp / "out"),
                )


if __name__ == "__main__":
    unittest.main()
