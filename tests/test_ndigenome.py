from __future__ import annotations

import argparse
import csv
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


def query_length(cigar: list[tuple[int, int]]) -> int:
    return sum(length for operation, length in cigar if operation in (0, 1, 4, 7, 8))


class FilterReasonTests(unittest.TestCase):
    def test_control_filter_reports_each_failed_criterion(self) -> None:
        args = argparse.Namespace(
            analysis="ndigenome",
            max_softclip_fraction=0.20,
            max_indel_fraction=0.20,
            min_support_mean_mapq=10.0,
            control_max_fraction=0.05,
            control_min_fold=5.0,
            control_max_q=0.05,
        )
        cases = {
            "HIGH_CONTROL_FRACTION": {
                "control_fraction": 0.10,
                "control_fold_enrichment": 10.0,
                "control_fisher_q": 0.01,
            },
            "LOW_CONTROL_FOLD": {
                "control_fraction": 0.01,
                "control_fold_enrichment": 2.0,
                "control_fisher_q": 0.01,
            },
            "CONTROL_Q_FAIL": {
                "control_fraction": 0.01,
                "control_fold_enrichment": 10.0,
                "control_fisher_q": 0.10,
            },
        }
        for expected_reason, control_values in cases.items():
            with self.subTest(expected_reason=expected_reason):
                row = {
                    "softclip_fraction": 0.0,
                    "indel_fraction": 0.0,
                    "known_indel_overlap": "",
                    "support_mean_mapq": 60.0,
                    "control_status": "MATCHED_CONTROL",
                    "signal_classification": "SSB",
                    "_caller_filter_reasons": [],
                    **control_values,
                }

                cleavage.apply_filters_to_row(row, args)

                self.assertEqual(row["classification"], "ARTIFACT_RISK")
                self.assertEqual(row["filter_reasons"], expected_reason)
                self.assertNotIn(
                    "CONTROL_NOT_ENRICHED",
                    row["filter_reasons"],
                )


@unittest.skipIf(pysam is None, "pysam is not installed")
class NDigenomeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tempdir.name)
        self.header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 2000}],
        }

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def make_read(
        self,
        name: str,
        start: int,
        cigar: list[tuple[int, int]] | None = None,
        reverse: bool = False,
        extra_flag: int = 0,
        mapq: int = 60,
        nm: int = 0,
    ):
        read = pysam.AlignedSegment()
        read.query_name = name
        read.query_sequence = "A" * query_length(cigar or [(0, 50)])
        read.flag = 1 | (16 if reverse else 0) | extra_flag
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = mapq
        read.cigartuples = cigar or [(0, 50)]
        read.next_reference_id = 0
        read.next_reference_start = start + 100
        read.template_length = 150
        read.query_qualities = pysam.qualitystring_to_array(
            "I" * len(read.query_sequence)
        )
        read.set_tag("NM", nm)
        return read

    def write_bam(self, name: str, reads: list) -> Path:
        bam = self.tmp / f"{name}.bam"
        reads.sort(key=lambda read: (read.reference_id, read.reference_start, read.flag))
        with pysam.AlignmentFile(bam, "wb", header=self.header) as output:
            for read in reads:
                output.write(read)
        pysam.index(str(bam))
        return bam

    def background_covering(self, position: int, count: int) -> list:
        return [
            self.make_read(f"background_{position}_{index}", position - 20 - index)
            for index in range(count)
        ]

    def reverse_background_covering(self, position: int, count: int) -> list:
        return [
            self.make_read(
                f"reverse_background_{position}_{index}",
                position - 49 + index,
                reverse=True,
            )
            for index in range(count)
        ]

    def args(
        self,
        bam: Path,
        prefix: str,
        control: Path | None = None,
        vcf: Path | None = None,
        analysis: str = "ndigenome",
    ) -> argparse.Namespace:
        return argparse.Namespace(
            analysis=analysis,
            bam=str(bam),
            sample=prefix,
            output_prefix=str(self.tmp / prefix),
            control_bam=str(control) if control else None,
            control_sample="Control" if control else "",
            variant_vcf=str(vcf) if vcf else None,
            genome_blacklist=None,
            keep_multimappers=False,
            ndigenome_min_count=5,
            ndigenome_min_fraction=0.20,
            ndigenome_min_mapq=1,
            ndigenome_opposite_window=5,
            ndigenome_ambiguous_min_count=2,
            ndigenome_ambiguous_min_fraction=0.05,
            digenome_overhang=0,
            digenome_pair_window=2,
            digenome_min_mapq=1,
            digenome_forward_cutoff=4,
            digenome_reverse_cutoff=4,
            digenome_depth_cutoff=4,
            digenome_fraction_cutoff=0.20,
            digenome_score_cutoff=1.0,
            artifact_window=10,
            max_softclip_fraction=0.20,
            max_indel_fraction=0.20,
            min_support_mean_mapq=10.0,
            control_min_depth=1,
            control_max_fraction=0.05,
            control_min_fold=3.0,
            control_max_q=0.05,
        )

    def read_rows(
        self, prefix: str, analysis: str = "ndigenome"
    ) -> list[dict[str, str]]:
        path = self.tmp / f"{prefix}.{analysis}.all.tsv"
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def read_tier_rows(
        self,
        prefix: str,
        tier: str,
        analysis: str = "ndigenome",
    ) -> list[dict[str, str]]:
        path = self.tmp / f"{prefix}.{analysis}.{tier}.tsv"
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_endpoint_coordinates_include_complex_cigars(self) -> None:
        forward = self.make_read("forward", 100, [(4, 5), (0, 45)])
        reverse = self.make_read(
            "reverse", 100, [(5, 2), (0, 20), (2, 3), (0, 20), (4, 5)], reverse=True
        )
        self.assertEqual(cleavage.endpoint_position(forward), 100)
        self.assertEqual(cleavage.five_prime_softclip_length(forward), 5)
        self.assertEqual(cleavage.endpoint_position(reverse), 142)
        self.assertEqual(cleavage.five_prime_softclip_length(reverse), 5)

    def test_true_single_strand_signal_passes(self) -> None:
        reads = [
            self.make_read(f"signal_{index}", 100) for index in range(8)
        ] + self.background_covering(100, 12)
        bam = self.write_bam("ssb", reads)
        cleavage.run_serial_calling(self.args(bam, "ssb"))
        rows = self.read_rows("ssb")
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["classification"], "SSB")
        self.assertEqual(rows[0]["filter_status"], "PASS")
        self.assertEqual(
            len(self.read_tier_rows("ssb", "high_confidence")),
            1,
        )
        self.assertEqual(self.read_tier_rows("ssb", "manual_review"), [])
        self.assertEqual(self.read_tier_rows("ssb", "artifact"), [])

    def test_opposite_strand_signal_is_possible_dsb(self) -> None:
        reads = [self.make_read(f"forward_{index}", 100) for index in range(8)]
        reads += [
            self.make_read(f"reverse_{index}", 52, reverse=True) for index in range(8)
        ]
        reads += self.background_covering(100, 8)
        bam = self.write_bam("dsb", reads)
        cleavage.run_serial_calling(self.args(bam, "dsb"))
        rows = self.read_rows("dsb")
        self.assertTrue(rows)
        self.assertTrue(
            all(row["signal_classification"] == "POSSIBLE_DSB" for row in rows)
        )
        self.assertTrue(all(row["filter_status"] == "FILTERED" for row in rows))

    def test_blacklisted_opposite_endpoint_does_not_classify_focal_site(self) -> None:
        reads = [self.make_read(f"forward_{index}", 100) for index in range(8)]
        reads += [
            self.make_read(f"reverse_{index}", 53, reverse=True)
            for index in range(8)
        ]
        reads += self.background_covering(100, 8)
        bam = self.write_bam("blacklisted_opposite", reads)
        blacklist = self.tmp / "blacklisted_opposite.bed"
        blacklist.write_text("chr1\t102\t103\n")
        args = self.args(bam, "blacklisted_opposite")
        args.genome_blacklist = str(blacklist)

        cleavage.run_serial_calling(args)

        rows = self.read_rows("blacklisted_opposite")
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["position_0based"], "100")
        self.assertEqual(rows[0]["signal_classification"], "SSB")
        self.assertEqual(rows[0]["opposite_position_0based"], "")
        self.assertEqual(rows[0]["opposite_count"], "0")

    def test_opposite_signal_prefers_threshold_passing_fraction(self) -> None:
        reads = [
            self.make_read(f"forward_{index}", 100)
            for index in range(8)
        ]
        reads += [
            self.make_read(f"reverse_100_{index}", 51, reverse=True)
            for index in range(10)
        ]
        reads += [
            self.make_read(f"reverse_101_{index}", 52, reverse=True)
            for index in range(11)
        ]
        reads += [
            self.make_read(f"reverse_background_{index}", 101, reverse=True)
            for index in range(100)
        ]
        bam_path = self.write_bam("opposite_ranking", reads)

        with pysam.AlignmentFile(bam_path, "rb") as bam:
            position, metrics = cleavage.find_best_opposite_signal(
                bam=bam,
                contig="chr1",
                position=100,
                strand="+",
                window=5,
                artifact_window=10,
                min_mapq=1,
                primary_min_count=10,
                primary_min_fraction=0.20,
                ambiguous_min_count=2,
                ambiguous_min_fraction=0.05,
            )

        self.assertEqual(position, 100)
        self.assertEqual(metrics["endpoint_count"], 10)
        self.assertGreater(metrics["endpoint_fraction"], 0.20)

    def test_weak_opposite_signal_is_ambiguous(self) -> None:
        reads = [self.make_read(f"forward_{index}", 100) for index in range(8)]
        reads += [
            self.make_read(f"weak_reverse_{index}", 52, reverse=True)
            for index in range(2)
        ]
        reads += self.background_covering(100, 12)
        bam = self.write_bam("ambiguous", reads)
        cleavage.run_serial_calling(self.args(bam, "ambiguous"))
        rows = self.read_rows("ambiguous")
        self.assertEqual(rows[0]["signal_classification"], "AMBIGUOUS")
        self.assertEqual(
            len(self.read_tier_rows("ambiguous", "manual_review")),
            1,
        )
        self.assertEqual(self.read_tier_rows("ambiguous", "artifact"), [])

    def test_duplicate_secondary_and_low_mapq_reads_are_excluded(self) -> None:
        reads = [
            self.make_read(f"duplicate_{index}", 100, extra_flag=1024)
            for index in range(8)
        ]
        reads += [
            self.make_read(f"secondary_{index}", 100, extra_flag=256)
            for index in range(8)
        ]
        reads += [
            self.make_read(f"supplementary_{index}", 100, extra_flag=2048)
            for index in range(8)
        ]
        reads += [
            self.make_read(f"low_mapq_{index}", 100, mapq=0) for index in range(8)
        ]
        reads += self.background_covering(100, 10)
        bam = self.write_bam("filtered", reads)
        qc = cleavage.run_serial_calling(self.args(bam, "filtered"))
        self.assertEqual(qc["reported_candidates"], 0)

    def test_mapq_zero_primaries_are_counted_without_secondary_inflation(
        self,
    ) -> None:
        reads = [
            self.make_read(f"primary_{index}", 100, mapq=0)
            for index in range(8)
        ]
        reads += [
            self.make_read(
                f"secondary_{index}",
                100,
                extra_flag=256,
                mapq=0,
            )
            for index in range(8)
        ]
        reads += self.background_covering(100, 12)
        bam = self.write_bam("multimappers", reads)
        args = self.args(bam, "multimappers")
        args.ndigenome_min_mapq = 0
        args.min_support_mean_mapq = 0
        args.keep_multimappers = True

        qc = cleavage.run_serial_calling(args)

        rows = self.read_rows("multimappers")
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["endpoint_count"], "8")
        self.assertEqual(rows[0]["secondary_endpoint_count"], "8")
        self.assertEqual(rows[0]["filter_status"], "PASS")
        self.assertTrue(qc["alignment_counting"]["keep_multimappers"])
        self.assertTrue(qc["warnings"])

    def test_random_endpoints_do_not_create_candidates(self) -> None:
        reads = [
            self.make_read(f"random_{index}", 100 + index * 3)
            for index in range(30)
        ]
        bam = self.write_bam("random", reads)
        qc = cleavage.run_serial_calling(self.args(bam, "random"))
        self.assertEqual(qc["reported_candidates"], 0)

    def test_soft_clipping_creates_artifact_risk(self) -> None:
        reads = [
            self.make_read(f"softclip_{index}", 200, [(4, 5), (0, 45)])
            for index in range(8)
        ]
        reads += self.background_covering(200, 10)
        bam = self.write_bam("softclip", reads)
        cleavage.run_serial_calling(self.args(bam, "softclip"))
        row = self.read_rows("softclip")[0]
        self.assertEqual(row["classification"], "ARTIFACT_RISK")
        self.assertIn("HIGH_5P_SOFTCLIP", row["filter_reasons"])
        self.assertEqual(
            self.read_tier_rows("softclip", "manual_review"),
            [],
        )
        self.assertEqual(len(self.read_tier_rows("softclip", "artifact")), 1)

    def test_nearby_indel_and_vcf_are_annotated(self) -> None:
        reads = [
            self.make_read(f"indel_{index}", 300, [(0, 5), (1, 1), (0, 44)])
            for index in range(8)
        ]
        reads += self.background_covering(300, 8)
        bam = self.write_bam("indel", reads)

        plain_vcf = self.tmp / "variants.vcf"
        plain_vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=2000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t304\t.\tAA\tA\t60\tPASS\t.\n"
        )
        compressed_vcf = self.tmp / "variants.vcf.gz"
        pysam.tabix_compress(str(plain_vcf), str(compressed_vcf), force=True)
        pysam.tabix_index(str(compressed_vcf), preset="vcf", force=True)

        cleavage.run_serial_calling(
            self.args(bam, "indel", vcf=compressed_vcf)
        )
        row = self.read_rows("indel")[0]
        self.assertEqual(row["classification"], "ARTIFACT_RISK")
        self.assertIn("NEARBY_INDEL", row["filter_reasons"])
        self.assertIn("KNOWN_INDEL", row["filter_reasons"])

    def test_direct_caller_rejects_an_unindexed_vcf(self) -> None:
        reads = [
            self.make_read(f"signal_{index}", 350) for index in range(8)
        ] + self.background_covering(350, 8)
        bam = self.write_bam("unindexed_vcf", reads)
        vcf = self.tmp / "unindexed.vcf.gz"
        vcf.write_text("not-empty\n")
        with self.assertRaisesRegex(ValueError, "must have a .tbi or .csi"):
            cleavage.run_serial_calling(
                self.args(bam, "unindexed_vcf", vcf=vcf)
            )

    def test_matched_control_statistics_and_q_values(self) -> None:
        treated_reads = [
            self.make_read(f"treated_signal_{index}", 400) for index in range(10)
        ] + self.background_covering(400, 20)
        control_reads = self.background_covering(400, 30)
        treated = self.write_bam("treated", treated_reads)
        control = self.write_bam("control", control_reads)
        cleavage.run_serial_calling(
            self.args(treated, "controlled", control=control)
        )
        row = self.read_rows("controlled")[0]
        self.assertEqual(row["control_status"], "MATCHED_CONTROL")
        self.assertLess(float(row["control_fisher_q"]), 0.05)
        self.assertEqual(row["filter_status"], "PASS")

    def test_zero_control_depth_is_reported_as_insufficient(self) -> None:
        treated_reads = [
            self.make_read(f"treated_signal_{index}", 450)
            for index in range(10)
        ] + self.background_covering(450, 20)
        control_reads = [self.make_read("control_elsewhere", 100)]
        treated = self.write_bam("zero_depth_treated", treated_reads)
        control = self.write_bam("zero_depth_control", control_reads)

        cleavage.run_serial_calling(
            self.args(
                treated,
                "zero_control_depth",
                control=control,
            )
        )

        row = self.read_rows("zero_control_depth")[0]
        self.assertEqual(
            row["control_status"],
            "INSUFFICIENT_CONTROL_COVERAGE",
        )
        self.assertEqual(row["control_depth"], "0")
        self.assertEqual(row["control_fold_enrichment"], "")
        self.assertEqual(row["control_fisher_q"], "")
        self.assertIn(
            "INSUFFICIENT_CONTROL_COVERAGE",
            row["filter_reasons"],
        )

    def test_variant_vcf_contigs_must_match_the_bam(self) -> None:
        reads = [
            self.make_read(f"signal_{index}", 475)
            for index in range(8)
        ] + self.background_covering(475, 8)
        bam = self.write_bam("vcf_contig_mismatch", reads)

        plain_vcf = self.tmp / "mismatched.vcf"
        plain_vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=1,length=2000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t100\t.\tAA\tA\t60\tPASS\t.\n"
        )
        compressed_vcf = self.tmp / "mismatched.vcf.gz"
        pysam.tabix_compress(
            str(plain_vcf),
            str(compressed_vcf),
            force=True,
        )
        pysam.tabix_index(
            str(compressed_vcf),
            preset="vcf",
            force=True,
        )

        with self.assertRaisesRegex(
            ValueError,
            "missing analyzed BAM contig",
        ):
            cleavage.run_serial_calling(
                self.args(
                    bam,
                    "vcf_contig_mismatch",
                    vcf=compressed_vcf,
                )
            )

    def test_benjamini_hochberg_is_monotonic(self) -> None:
        adjusted = cleavage.benjamini_hochberg([0.01, 0.04, 0.03, 0.20])
        self.assertEqual(len(adjusted), 4)
        self.assertLessEqual(adjusted[0], adjusted[2])
        self.assertTrue(all(0 <= value <= 1 for value in adjusted))

    def test_direct_cli_uses_current_option_names(self) -> None:
        args = cleavage.build_parser().parse_args(
            [
                "--analysis",
                "ndigenome",
                "--bam",
                "input.bam",
                "--sample",
                "Sample",
                "--output-prefix",
                "Sample",
                "--ndigenome-min-count",
                "7",
                "--ndigenome-min-fraction",
                "0.3",
                "--ndigenome-min-mapq",
                "2",
                "--ndigenome-opposite-window",
                "4",
                "--ndigenome-ambiguous-min-count",
                "2",
                "--ndigenome-ambiguous-min-fraction",
                "0.1",
            ]
        )
        self.assertEqual(args.ndigenome_min_count, 7)
        self.assertEqual(args.ndigenome_min_fraction, 0.3)
        self.assertEqual(args.ndigenome_min_mapq, 2)
        self.assertEqual(args.ndigenome_opposite_window, 4)

    def test_digenome_paired_endpoints_pass(self) -> None:
        reads = [
            self.make_read(f"forward_{index}", 500) for index in range(8)
        ]
        reads += [
            self.make_read(f"reverse_{index}", 451, reverse=True)
            for index in range(8)
        ]
        bam = self.write_bam("digenome_pair", reads)
        qc = cleavage.run_serial_calling(
            self.args(bam, "digenome_pair", analysis="digenome")
        )
        rows = self.read_rows("digenome_pair", "digenome")
        self.assertEqual(qc["high_confidence_dsb"], 1)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["classification"], "DSB")
        self.assertEqual(rows[0]["filter_status"], "PASS")
        self.assertEqual(rows[0]["forward_position_0based"], "500")
        self.assertEqual(rows[0]["reverse_position_0based"], "500")

    def test_digenome_pair_is_excluded_when_one_endpoint_is_blacklisted(
        self,
    ) -> None:
        reads = [
            self.make_read(f"forward_{index}", 500) for index in range(8)
        ]
        reads += [
            self.make_read(f"reverse_{index}", 449, reverse=True)
            for index in range(8)
        ]
        bam = self.write_bam("digenome_one_blacklisted_endpoint", reads)
        cleavage.run_serial_calling(
            self.args(
                bam,
                "digenome_one_blacklisted_endpoint_baseline",
                analysis="digenome",
            )
        )
        baseline_rows = self.read_rows(
            "digenome_one_blacklisted_endpoint_baseline",
            "digenome",
        )
        self.assertEqual(len(baseline_rows), 1)
        self.assertEqual(
            baseline_rows[0]["reverse_position_0based"],
            "498",
        )

        blacklist = self.tmp / "digenome_one_blacklisted_endpoint.bed"
        blacklist.write_text("chr1\t498\t499\n")
        args = self.args(
            bam,
            "digenome_one_blacklisted_endpoint_masked",
            analysis="digenome",
        )
        args.genome_blacklist = str(blacklist)
        qc = cleavage.run_serial_calling(args)

        self.assertEqual(qc["reported_candidates"], 0)
        self.assertEqual(
            self.read_rows(
                "digenome_one_blacklisted_endpoint_masked",
                "digenome",
            ),
            [],
        )

    def test_digenome_rejects_one_strand_only_signal(self) -> None:
        reads = [
            self.make_read(f"forward_{index}", 550) for index in range(8)
        ]
        bam = self.write_bam("digenome_one_strand", reads)
        qc = cleavage.run_serial_calling(
            self.args(bam, "digenome_one_strand", analysis="digenome")
        )
        self.assertEqual(qc["reported_candidates"], 0)
        self.assertEqual(
            self.read_rows("digenome_one_strand", "digenome"), []
        )

    def test_digenome_overhang_pairs_expected_coordinates(self) -> None:
        reads = [
            self.make_read(f"forward_{index}", 600) for index in range(8)
        ]
        reads += [
            self.make_read(f"reverse_{index}", 547, reverse=True)
            for index in range(8)
        ]
        bam = self.write_bam("digenome_overhang", reads)
        args = self.args(bam, "digenome_overhang", analysis="digenome")
        args.digenome_overhang = 4
        args.digenome_pair_window = 0
        cleavage.run_serial_calling(args)
        row = self.read_rows("digenome_overhang", "digenome")[0]
        self.assertEqual(row["forward_position_0based"], "600")
        self.assertEqual(row["reverse_position_0based"], "596")

    def test_digenome_score_formula(self) -> None:
        score = cleavage.digenome_score(8, 0.5, 12, 0.25)
        self.assertAlmostEqual(score, 0.625)

    def test_digenome_pairing_prefers_a_threshold_passing_pair(self) -> None:
        reads = [
            self.make_read(f"forward_{index}", 1000) for index in range(8)
        ]
        reads += self.background_covering(1000, 12)
        reads += [
            self.make_read(f"reverse_exact_{index}", 951, reverse=True)
            for index in range(8)
        ]
        reads += [
            self.make_read(f"reverse_nearby_{index}", 949, reverse=True)
            for index in range(8)
        ]
        bam = self.write_bam("digenome_pair_priority", reads)
        args = self.args(
            bam, "digenome_pair_priority", analysis="digenome"
        )
        args.digenome_depth_cutoff = 10
        args.digenome_score_cutoff = 0.5
        cleavage.run_serial_calling(args)
        row = self.read_rows("digenome_pair_priority", "digenome")[0]
        self.assertEqual(row["reverse_position_0based"], "998")
        self.assertEqual(row["filter_status"], "PASS")

    def test_digenome_pairing_maximizes_nonconflicting_pairs(self) -> None:
        metrics = cleavage.empty_metrics()
        candidates = [
            cleavage.DigenomePairCandidate(
                "chr1",
                100,
                100,
                metrics,
                metrics,
                10.0,
                [],
            ),
            cleavage.DigenomePairCandidate(
                "chr1",
                100,
                101,
                metrics,
                metrics,
                9.0,
                [],
            ),
            cleavage.DigenomePairCandidate(
                "chr1",
                101,
                100,
                metrics,
                metrics,
                9.0,
                [],
            ),
        ]

        selected = cleavage.select_digenome_pairs(candidates)

        self.assertEqual(
            [
                (candidate.forward_position, candidate.reverse_position)
                for candidate in selected
            ],
            [(100, 101), (101, 100)],
        )

    def test_digenome_control_evidence_can_filter_a_pair(self) -> None:
        treated_reads = [
            self.make_read(f"treated_forward_{index}", 700)
            for index in range(10)
        ]
        treated_reads += [
            self.make_read(f"treated_reverse_{index}", 651, reverse=True)
            for index in range(10)
        ]
        control_reads = [
            self.make_read(f"control_forward_{index}", 700)
            for index in range(10)
        ]
        control_reads += [
            self.make_read(f"control_reverse_{index}", 651, reverse=True)
            for index in range(10)
        ]
        treated = self.write_bam("digenome_treated", treated_reads)
        control = self.write_bam("digenome_control", control_reads)
        cleavage.run_serial_calling(
            self.args(
                treated,
                "digenome_controlled",
                control=control,
                analysis="digenome",
            )
        )
        row = self.read_rows("digenome_controlled", "digenome")[0]
        self.assertEqual(row["control_status"], "MATCHED_CONTROL")
        self.assertEqual(row["classification"], "ARTIFACT_RISK")
        self.assertIn("HIGH_CONTROL_FRACTION", row["filter_reasons"])
        self.assertIn("LOW_CONTROL_FOLD", row["filter_reasons"])
        self.assertIn("CONTROL_Q_FAIL", row["filter_reasons"])

    def test_digenome_indel_and_vcf_artifacts_are_shared(self) -> None:
        forward_cigar = [(0, 5), (1, 1), (0, 44)]
        reads = [
            self.make_read(
                f"forward_indel_{index}", 800, cigar=forward_cigar
            )
            for index in range(8)
        ]
        reads += [
            self.make_read(f"reverse_{index}", 751, reverse=True)
            for index in range(8)
        ]
        bam = self.write_bam("digenome_indel", reads)

        plain_vcf = self.tmp / "digenome_variants.vcf"
        plain_vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=2000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t805\t.\tAA\tA\t60\tPASS\t.\n"
        )
        compressed_vcf = self.tmp / "digenome_variants.vcf.gz"
        pysam.tabix_compress(
            str(plain_vcf), str(compressed_vcf), force=True
        )
        pysam.tabix_index(str(compressed_vcf), preset="vcf", force=True)

        cleavage.run_serial_calling(
            self.args(
                bam,
                "digenome_indel",
                vcf=compressed_vcf,
                analysis="digenome",
            )
        )
        row = self.read_rows("digenome_indel", "digenome")[0]
        self.assertEqual(row["classification"], "ARTIFACT_RISK")
        self.assertIn("NEARBY_INDEL", row["filter_reasons"])
        self.assertIn("KNOWN_INDEL", row["filter_reasons"])

if __name__ == "__main__":
    unittest.main()
