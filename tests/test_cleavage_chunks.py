from __future__ import annotations

import argparse
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
import call_cleavage as cleavage  # noqa: E402
import finalize_cleavage_chunks as finalizer  # noqa: E402
import plan_cleavage_chunks as planner  # noqa: E402


@unittest.skipIf(pysam is None, "pysam is not installed")
class CleavageChunkTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tempdir.name)
        self.header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [
                {"SN": "chr1", "LN": 2000},
                {"SN": "chr2", "LN": 2000},
                {"SN": "chr3", "LN": 2000},
            ],
        }

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def make_read(
        self,
        name: str,
        reference_id: int,
        start: int,
        *,
        reverse: bool = False,
        mapq: int = 60,
        extra_flag: int = 0,
    ):
        read = pysam.AlignedSegment()
        read.query_name = name
        read.query_sequence = "A" * 50
        read.flag = 1 | (16 if reverse else 0) | extra_flag
        read.reference_id = reference_id
        read.reference_start = start
        read.mapping_quality = mapq
        read.cigartuples = [(0, 50)]
        read.next_reference_id = reference_id
        read.next_reference_start = start + 100
        read.template_length = 150
        read.query_qualities = pysam.qualitystring_to_array("I" * 50)
        read.set_tag("NM", 0)
        return read

    def write_bam(self, name: str, reads: list) -> Path:
        bam = self.tmp / f"{name}.bam"
        reads.sort(
            key=lambda read: (
                read.reference_id,
                read.reference_start,
                read.flag,
                read.query_name,
            )
        )
        with pysam.AlignmentFile(bam, "wb", header=self.header) as output:
            for read in reads:
                output.write(read)
        pysam.index(str(bam))
        return bam

    def endpoint_site_reads(
        self,
        reference_id: int,
        position: int,
        signal_count: int,
        prefix: str,
        *,
        mapq: int = 0,
    ) -> list:
        reads = [
            self.make_read(
                f"{prefix}_forward_signal_{index}",
                reference_id,
                position,
                mapq=mapq,
            )
            for index in range(signal_count)
        ]
        reads += [
            self.make_read(
                f"{prefix}_reverse_signal_{index}",
                reference_id,
                position - 49,
                reverse=True,
                mapq=mapq,
            )
            for index in range(signal_count)
        ]
        reads += [
            self.make_read(
                f"{prefix}_forward_background_{index}",
                reference_id,
                position - 20 - index,
                mapq=mapq,
            )
            for index in range(12)
        ]
        reads += [
            self.make_read(
                f"{prefix}_reverse_background_{index}",
                reference_id,
                position - 48 + index,
                reverse=True,
                mapq=mapq,
            )
            for index in range(12)
        ]
        reads.append(
            self.make_read(
                f"{prefix}_secondary",
                reference_id,
                position,
                mapq=mapq,
                extra_flag=256,
            )
        )
        return reads

    def control_site_reads(
        self,
        reference_id: int,
        position: int,
        prefix: str,
    ) -> list:
        reads = [
            self.make_read(
                f"{prefix}_forward_{index}",
                reference_id,
                position - 20 - index,
                mapq=0,
            )
            for index in range(20)
        ]
        reads += [
            self.make_read(
                f"{prefix}_reverse_{index}",
                reference_id,
                position - 48 + index,
                reverse=True,
                mapq=0,
            )
            for index in range(20)
        ]
        return reads

    def write_variant_vcf(self) -> Path:
        plain = self.tmp / "variants.vcf"
        plain.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=2000>\n"
            "##contig=<ID=chr2,length=2000>\n"
            "##contig=<ID=chr3,length=2000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr2\t501\t.\tAA\tA\t60\tPASS\t.\n"
        )
        compressed = self.tmp / "variants.vcf.gz"
        pysam.tabix_compress(str(plain), str(compressed), force=True)
        pysam.tabix_index(str(compressed), preset="vcf", force=True)
        return compressed

    def caller_args(
        self,
        bam: Path,
        control: Path,
        vcf: Path,
        prefix: Path,
        analysis: str,
    ) -> argparse.Namespace:
        return argparse.Namespace(
            analysis=analysis,
            bam=str(bam),
            sample="Sample",
            output_prefix=str(prefix),
            control_bam=str(control),
            control_sample="Control",
            variant_vcf=str(vcf),
            keep_multimappers=True,
            intervals_file=None,
            chunk_id=None,
            raw_output=None,
            summary_output=None,
            ndigenome_min_count=5,
            ndigenome_min_fraction=0.20,
            ndigenome_min_mapq=0,
            ndigenome_opposite_window=2,
            ndigenome_ambiguous_min_count=2,
            ndigenome_ambiguous_min_fraction=0.05,
            digenome_overhang=0,
            digenome_pair_window=2,
            digenome_min_mapq=0,
            digenome_forward_cutoff=4,
            digenome_reverse_cutoff=4,
            digenome_depth_cutoff=4,
            digenome_fraction_cutoff=0.20,
            digenome_score_cutoff=0.20,
            artifact_window=10,
            max_softclip_fraction=0.20,
            max_indel_fraction=0.20,
            min_support_mean_mapq=0.0,
            control_min_depth=1,
            control_max_fraction=0.05,
            control_min_fold=3.0,
            control_max_q=0.05,
        )

    def test_planner_balances_contigs_and_caps_chunk_count(self) -> None:
        reads = []
        for reference_id, count in enumerate((8, 5, 3)):
            reads.extend(
                self.make_read(
                    f"contig_{reference_id}_{index}",
                    reference_id,
                    100 + index,
                )
                for index in range(count)
            )
        bam = self.write_bam("planner", reads)

        chunks = planner.plan_chunks(str(bam), 2)
        self.assertEqual(
            [
                (
                    chunk["mapped_records"],
                    [
                        (
                            interval["contig"],
                            interval["start"],
                            interval["end"],
                        )
                        for interval in chunk["intervals"]
                    ],
                )
                for chunk in chunks
            ],
            [
                (8, [("chr1", 0, 2000)]),
                (8, [("chr2", 0, 2000), ("chr3", 0, 2000)]),
            ],
        )

        capped = planner.plan_chunks(str(bam), 20)
        self.assertEqual(len(capped), 16)
        self.assertTrue(all(chunk["intervals"] for chunk in capped))
        self.assertTrue(
            any(
                interval["start"] > 0 or interval["end"] < 2000
                for chunk in capped
                for interval in chunk["intervals"]
            )
        )

        output_dir = self.tmp / "plan"
        plan_path = self.tmp / "chunks.tsv"
        padded_chunks = planner.plan_chunks(str(bam), 2, padding=10)
        planner.write_plan(padded_chunks, output_dir, plan_path)
        self.assertEqual(
            sorted(path.name for path in output_dir.iterdir()),
            ["chunk_000.intervals.tsv", "chunk_001.intervals.tsv"],
        )
        intervals = cleavage.read_intervals_file(
            str(output_dir / "chunk_000.intervals.tsv")
        )
        self.assertEqual(intervals[0].scan_start, 0)
        self.assertEqual(intervals[0].scan_end, 2000)

    def test_overlapping_owned_intervals_are_rejected(self) -> None:
        intervals = self.tmp / "overlap.intervals.tsv"
        intervals.write_text(
            "contig\towner_start\towner_end\tscan_start\tscan_end\n"
            "chr1\t0\t1100\t0\t1110\n"
            "chr1\t1000\t2000\t990\t2000\n"
        )
        with self.assertRaisesRegex(ValueError, "overlapping"):
            cleavage.read_intervals_file(str(intervals))

    def test_finalizer_rejects_missing_chunk_intervals(self) -> None:
        summary = self.tmp / "missing.chunk.json"
        summary.write_text(
            json.dumps(
                {
                    "sample": "Sample",
                    "analysis": "ndigenome",
                    "chunk_id": "chunk_000",
                    "contigs": ["chr1"],
                    "intervals": [
                        {"contig": "chr1", "start": 0, "end": 2000}
                    ],
                    "expected_contigs": ["chr1", "chr2"],
                    "expected_contig_lengths": {
                        "chr1": 2000,
                        "chr2": 2000,
                    },
                    "parameters": {},
                }
            )
        )
        with self.assertRaisesRegex(ValueError, "missing: chr2"):
            finalizer.load_summaries(
                [str(summary)],
                "Sample",
                "ndigenome",
            )

    def test_finalizer_rejects_overlapping_chunk_intervals(self) -> None:
        summaries = []
        for chunk_id, start, end in (
            ("chunk_000", 0, 1200),
            ("chunk_001", 1000, 2000),
        ):
            path = self.tmp / f"{chunk_id}.json"
            path.write_text(
                json.dumps(
                    {
                        "sample": "Sample",
                        "analysis": "ndigenome",
                        "chunk_id": chunk_id,
                        "contigs": ["chr1"],
                        "intervals": [
                            {
                                "contig": "chr1",
                                "start": start,
                                "end": end,
                            }
                        ],
                        "expected_contigs": ["chr1"],
                        "expected_contig_lengths": {"chr1": 2000},
                        "parameters": {},
                    }
                )
            )
            summaries.append(str(path))
        with self.assertRaisesRegex(ValueError, "overlaps"):
            finalizer.load_summaries(
                summaries,
                "Sample",
                "ndigenome",
            )

    def test_finalizer_rejects_a_row_outside_chunk_ownership(self) -> None:
        fragment = self.tmp / "outside.raw.jsonl.gz"
        with gzip.open(fragment, "wt") as handle:
            handle.write(
                json.dumps(
                    {
                        "_chunk_id": "chunk_000",
                        "contig": "chr1",
                        "position_0based": 1200,
                        "strand": "+",
                    }
                )
                + "\n"
            )
        summaries = [
            {
                "chunk_id": "chunk_000",
                "intervals": [
                    {"contig": "chr1", "start": 0, "end": 1000}
                ],
            }
        ]
        connection = finalizer.create_database(":memory:")
        try:
            with self.assertRaisesRegex(ValueError, "outside"):
                finalizer.load_fragments(
                    connection,
                    [str(fragment)],
                    summaries,
                )
        finally:
            connection.close()

    def test_boundary_padding_preserves_calls_without_duplicates(self) -> None:
        forward_position = 1001
        reverse_position = 999
        treated_reads = [
            self.make_read(
                f"boundary_forward_signal_{index}",
                0,
                forward_position,
                mapq=0,
            )
            for index in range(8)
        ]
        treated_reads += [
            self.make_read(
                f"boundary_reverse_signal_{index}",
                0,
                reverse_position - 49,
                reverse=True,
                mapq=0,
            )
            for index in range(8)
        ]
        treated_reads += [
            self.make_read(
                f"boundary_forward_background_{index}",
                0,
                forward_position - 20 - index,
                mapq=0,
            )
            for index in range(12)
        ]
        treated_reads += [
            self.make_read(
                f"boundary_reverse_background_{index}",
                0,
                reverse_position - 48 + index,
                reverse=True,
                mapq=0,
            )
            for index in range(12)
        ]
        treated = self.write_bam("boundary_treated", treated_reads)
        vcf = self.write_variant_vcf()
        chunks = planner.plan_chunks(str(treated), 2, padding=10)
        self.assertEqual(
            [
                (
                    interval["start"],
                    interval["end"],
                    interval["scan_start"],
                    interval["scan_end"],
                )
                for chunk in chunks
                for interval in chunk["intervals"]
            ],
            [(0, 1000, 0, 1010), (1000, 2000, 990, 2000)],
        )
        plan_dir = self.tmp / "boundary_plan"
        planner.write_plan(
            chunks,
            plan_dir,
            self.tmp / "boundary_plan.tsv",
        )

        for analysis in ("digenome", "ndigenome"):
            serial_prefix = self.tmp / f"boundary_serial_{analysis}"
            chunked_prefix = self.tmp / f"boundary_chunked_{analysis}"
            serial_args = self.caller_args(
                treated,
                Path(""),
                vcf,
                serial_prefix,
                analysis,
            )
            serial_args.control_bam = ""
            serial_args.control_sample = ""
            serial_args.digenome_overhang = 2
            serial_qc = cleavage.run_serial_calling(serial_args)

            raw_fragments = []
            summaries = []
            for interval_file in sorted(
                plan_dir.glob("*.intervals.tsv")
            ):
                chunk_id = interval_file.name.removesuffix(
                    ".intervals.tsv"
                )
                raw_path = self.tmp / (
                    f"boundary.{chunk_id}.{analysis}.raw.jsonl.gz"
                )
                summary_path = self.tmp / (
                    f"boundary.{chunk_id}.{analysis}.chunk.json"
                )
                args = self.caller_args(
                    treated,
                    Path(""),
                    vcf,
                    chunked_prefix,
                    analysis,
                )
                args.control_bam = ""
                args.control_sample = ""
                args.digenome_overhang = 2
                args.intervals_file = str(interval_file)
                args.chunk_id = chunk_id
                args.raw_output = str(raw_path)
                args.summary_output = str(summary_path)
                cleavage.run_chunk_calling(args)
                raw_fragments.append(str(raw_path))
                summaries.append(str(summary_path))

            chunked_qc = finalizer.finalize_chunks(
                argparse.Namespace(
                    sample="Sample",
                    analysis=analysis,
                    output_prefix=str(chunked_prefix),
                    raw_fragment=raw_fragments,
                    chunk_summary=summaries,
                )
            )
            self.assertEqual(
                Path(f"{serial_prefix}.{analysis}.all.tsv").read_bytes(),
                Path(f"{chunked_prefix}.{analysis}.all.tsv").read_bytes(),
            )
            self.assertEqual(
                serial_qc["candidate_endpoints_or_pairs_before_filters"],
                chunked_qc[
                    "candidate_endpoints_or_pairs_before_filters"
                ],
            )
            with Path(
                f"{chunked_prefix}.{analysis}.all.tsv"
            ).open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            keys = [
                (
                    row["contig"],
                    row["position_0based"],
                    row["strand"],
                )
                for row in rows
            ]
            self.assertEqual(len(keys), len(set(keys)))
            self.assertGreater(len(keys), 0)

    def test_serial_and_chunked_outputs_are_equivalent(self) -> None:
        treated_reads = self.endpoint_site_reads(0, 500, 8, "chr1")
        treated_reads += self.endpoint_site_reads(1, 700, 6, "chr2")
        treated_reads += [
            self.make_read(
                f"chr3_background_{index}",
                2,
                100 + index * 3,
                mapq=0,
            )
            for index in range(12)
        ]
        control_reads = self.control_site_reads(0, 500, "chr1_control")
        control_reads += self.control_site_reads(1, 700, "chr2_control")
        control_reads += [
            self.make_read(
                f"chr3_control_{index}",
                2,
                200 + index * 3,
                mapq=0,
            )
            for index in range(12)
        ]
        treated = self.write_bam("treated", treated_reads)
        control = self.write_bam("control", control_reads)
        vcf = self.write_variant_vcf()

        chunks = planner.plan_chunks(str(treated), 4, padding=10)
        plan_dir = self.tmp / "chunk_plan"
        planner.write_plan(chunks, plan_dir, self.tmp / "chunk_plan.tsv")
        interval_files = sorted(plan_dir.glob("*.intervals.tsv"))

        for analysis in ("digenome", "ndigenome"):
            serial_prefix = self.tmp / f"serial_{analysis}"
            chunked_prefix = self.tmp / f"chunked_{analysis}"
            serial_qc = cleavage.run_serial_calling(
                self.caller_args(
                    treated,
                    control,
                    vcf,
                    serial_prefix,
                    analysis,
                )
            )

            raw_fragments = []
            summaries = []
            for interval_file in interval_files:
                chunk_id = interval_file.name.removesuffix(".intervals.tsv")
                raw_path = (
                    self.tmp / f"{chunk_id}.{analysis}.raw.jsonl.gz"
                )
                summary_path = (
                    self.tmp / f"{chunk_id}.{analysis}.chunk.json"
                )
                args = self.caller_args(
                    treated,
                    control,
                    vcf,
                    chunked_prefix,
                    analysis,
                )
                args.intervals_file = str(interval_file)
                args.chunk_id = chunk_id
                args.raw_output = str(raw_path)
                args.summary_output = str(summary_path)
                cleavage.run_chunk_calling(args)
                raw_fragments.append(str(raw_path))
                summaries.append(str(summary_path))

            chunked_qc = finalizer.finalize_chunks(
                argparse.Namespace(
                    sample="Sample",
                    analysis=analysis,
                    output_prefix=str(chunked_prefix),
                    raw_fragment=raw_fragments,
                    chunk_summary=summaries,
                )
            )

            suffixes = [
                f".{analysis}.all.tsv",
                f".{analysis}.high_confidence.tsv",
                f".{analysis}.bed",
                f".{analysis}_mqc.tsv",
            ]
            for suffix in suffixes:
                self.assertEqual(
                    Path(f"{serial_prefix}{suffix}").read_bytes(),
                    Path(f"{chunked_prefix}{suffix}").read_bytes(),
                    suffix,
                )

            for key in (
                "candidate_endpoints_or_pairs_before_filters",
                "reported_candidates",
                "high_confidence_calls",
                "classifications",
                "signal_classifications",
            ):
                self.assertEqual(serial_qc[key], chunked_qc[key], key)

            with Path(
                f"{chunked_prefix}.{analysis}.all.tsv"
            ).open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            p_values = [float(row["control_fisher_p"]) for row in rows]
            q_values = [float(row["control_fisher_q"]) for row in rows]
            expected = cleavage.benjamini_hochberg(p_values)
            for observed, wanted in zip(q_values, expected):
                self.assertAlmostEqual(observed, wanted)

            qc_path = Path(f"{chunked_prefix}.{analysis}.qc.json")
            qc_document = json.loads(qc_path.read_text())
            self.assertEqual(qc_document["parameters"]["cleavage_chunks"], 4)


if __name__ == "__main__":
    unittest.main()
