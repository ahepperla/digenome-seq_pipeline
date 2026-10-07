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
    reverse_background,
    forward_background,
    write_bam,
    write_vcf,
    run_caller,
    make_settings,
    CHUNK_SETTINGS,
)

sys.path.insert(0, str(Path(__file__).resolve().parent / "golden"))
from scenarios import multi_contig  # noqa: E402
from cleavage.call import call_chunk
from cleavage.finalize import finalize
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
        """The chunk count changes runtime, never results. Uses the
        multi-contig golden scenario: controls, a known indel, a blacklisted
        site, and chr10 before chr3 in the header."""
        scenario = multi_contig()
        inputs = scenario.write_inputs(self.tmp / "inputs")
        for analysis, expected_rows in (("digenome", 3), ("ndigenome", 7)):
            settings = make_settings(analysis, base=scenario.settings)
            serial_prefix, serial_qc = run_caller(settings, inputs, self.tmp / f"{analysis}_1", chunks=1)
            self.assertEqual(serial_qc["rows"], expected_rows)
            for chunks in (2, 3, 4, 6):
                with self.subTest(analysis=analysis, chunks=chunks):
                    prefix, qc = run_caller(settings, inputs, self.tmp / f"{analysis}_{chunks}", chunks=chunks)
                    for suffix in ("all.tsv", "high_confidence.tsv", "manual_review.tsv", "artifact.tsv", "bed"):
                        self.assertEqual(
                            Path(f"{prefix}.{analysis}.{suffix}").read_bytes(),
                            Path(f"{serial_prefix}.{analysis}.{suffix}").read_bytes(),
                            suffix,
                        )
                    self.assertEqual(
                        Path(f"{prefix}.{analysis}_mqc.tsv").read_bytes(),
                        Path(f"{serial_prefix}.{analysis}_mqc.tsv").read_bytes(),
                    )
                    for key in ("rows", "tiers", "candidates_before_filters"):
                        self.assertEqual(qc[key], serial_qc[key], key)
            self.check_sample_wide_q_values(Path(f"{serial_prefix}.{analysis}.all.tsv"), expected_rows)

    def check_sample_wide_q_values(self, path: Path, expected_rows: int) -> None:
        """Every row has a matched control, and q is BH over all of them."""
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["control_status"] for row in rows], ["MATCHED_CONTROL"] * expected_rows)
        p_values = [float(row["control_fisher_p"]) for row in rows]
        q_values = [float(row["control_fisher_q"]) for row in rows]
        for observed, wanted in zip(q_values, benjamini_hochberg(p_values), strict=True):
            # p is written with 8 significant digits, so allow that rounding.
            self.assertAlmostEqual(observed, wanted, delta=1e-7 * wanted)

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

    def test_partner_found_across_a_boundary_at_overhang_distance(self) -> None:
        """Chunks look |overhang| + pair_window past their ranges for partners.
        With window 0 the partner sits exactly |overhang| bases away, across
        the 2-chunk boundary at 1000, for a positive and a negative overhang."""
        for overhang, forward, reverse in ((4, 1002, 998), (-4, 998, 1002)):
            with self.subTest(overhang=overhang):
                reads = [make_read(f"fwd_{i}", forward, mapq=0) for i in range(8)]
                reads += [reverse_read_ending_at(f"rev_{i}", reverse, mapq=0) for i in range(8)]
                reads += forward_background("fbg", forward, 12, mapq=0)
                reads += reverse_background("rbg", reverse, 12, mapq=0)
                directory = self.tmp / f"overhang_{overhang}"
                directory.mkdir()
                bam = write_bam(directory / "treated.bam", [("chr1", 2000)], reads)
                settings = make_settings(
                    "digenome", base=CHUNK_SETTINGS, digenome_overhang=overhang, digenome_pair_window=0
                )
                inputs = {"bam": bam, "control_bam": None, "vcf": None, "blacklist": None}
                serial, serial_qc = run_caller(settings, inputs, directory / "serial", chunks=1)
                chunked, chunked_qc = run_caller(settings, inputs, directory / "chunked", chunks=2)
                self.assertEqual((serial_qc["rows"], chunked_qc["rows"]), (1, 1))
                self.assertEqual(
                    Path(f"{chunked}.digenome.all.tsv").read_bytes(),
                    Path(f"{serial}.digenome.all.tsv").read_bytes(),
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
