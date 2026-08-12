from __future__ import annotations

import csv
from collections import Counter
import json
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class WorkflowStaticTests(unittest.TestCase):
    @staticmethod
    def read_fastq(path: Path) -> list[tuple[str, str, str]]:
        lines = path.read_text().splitlines()
        if len(lines) % 4:
            raise AssertionError(f"Malformed FASTQ fixture: {path}")
        return [
            (lines[index], lines[index + 1], lines[index + 3])
            for index in range(0, len(lines), 4)
        ]

    def test_analysis_modes_and_runtime_contract(self) -> None:
        main = (ROOT / "main.nf").read_text()
        config = (ROOT / "nextflow.config").read_text()
        base = (ROOT / "conf" / "base.config").read_text()
        self.assertIn("process PLAN_CLEAVAGE_CHUNKS", main)
        self.assertIn("process CLEAVAGE_CALL_CHUNK", main)
        self.assertIn("process FINALIZE_CLEAVAGE_CALL", main)
        self.assertIn("process PREFLIGHT", main)
        preflight_block = main.split(
            "process PREFLIGHT", 1
        )[1].split("process VALIDATE_SAMPLESHEET", 1)[0]
        self.assertIn("cache false", preflight_block)
        self.assertIn("preflight_ready_ch = PREFLIGHT.out.ready.first()", main)
        self.assertEqual(main.count("preflight_ready_ch"), 4)
        self.assertIn("withName: PREFLIGHT", base)
        self.assertIn(
            'path "${meta.sample}.${selected_analysis}.manual_review.tsv"',
            main,
        )
        self.assertIn(
            'path "${meta.sample}.${selected_analysis}.artifact.tsv"',
            main,
        )
        self.assertNotIn("process CLEAVAGE_CALL {", main)
        self.assertNotIn("process NDIGENOME_CALL", main)
        self.assertNotIn("process DIGENOME_DSB", main)
        self.assertIn("--analysis ${shellQuote(selected_analysis)}", main)
        self.assertIn("matched_requests_ch", main)
        self.assertIn("control_sample: row.control as String", main)
        self.assertIn("tuple val(meta)", main)
        self.assertNotIn("row.role", main)
        self.assertNotIn("row.analysis_group", main)
        self.assertNotIn("--analysis-group", main)
        self.assertNotIn("startsWith('NO_CONTROL')", main)
        self.assertNotIn("startsWith('NO_VARIANT')", main)
        self.assertNotIn("digenome_runtime", main + config)
        self.assertIn("container = params.containers.cleavage", base)
        self.assertNotIn("params.containers.digenome", base)
        self.assertNotIn("params.containers.ndigenome", base)
        self.assertIn("cleavage_pysam_v0.23.3.sif", config)
        self.assertNotIn("digenome.sif", config)

    def test_cleavage_chunks_default_and_cpu_scaling(self) -> None:
        main = (ROOT / "main.nf").read_text()
        config = (ROOT / "nextflow.config").read_text()
        base = (ROOT / "conf" / "base.config").read_text()
        self.assertIn("cleavage_chunks = 8", config)
        self.assertIn(
            "int cleavage_chunks = params.cleavage_chunks as int",
            main,
        )
        self.assertIn(
            "int cleavage_chunk_padding = [",
            main,
        )
        self.assertIn("--padding ${chunk_padding}", main)
        self.assertIn("--intervals-file", main)
        self.assertIn("genome_blacklist = null", config)
        self.assertIn("--genome-blacklist", main)
        self.assertIn("genome_blacklist_ch", main)
        self.assertIn("bin/genome_blacklist.py", main)
        self.assertIn('path("chunk_*.intervals.tsv")', main)
        self.assertIn("withName: CLEAVAGE_CALL_CHUNK", base)
        self.assertIn(
            "maxForks = params.cleavage_chunks as int",
            base,
        )
        chunk_block = base.split(
            "withName: CLEAVAGE_CALL_CHUNK", 1
        )[1].split("}", 1)[0]
        self.assertIn("cpus = 1", chunk_block)

    def test_parameter_reference_covers_every_configured_parameter(
        self,
    ) -> None:
        config = (ROOT / "nextflow.config").read_text()
        reference = (ROOT / "docs" / "parameters.md").read_text()
        readme = (ROOT / "README.md").read_text()
        schema = json.loads(
            (ROOT / "nextflow_schema.json").read_text()
        )
        params_block = config.split("params {", 1)[1].split("\n}", 1)[0]
        configured = set(
            re.findall(
                r"^    ([a-z][a-z0-9_]*)\s*=",
                params_block,
                flags=re.MULTILINE,
            )
        )
        documented_names = re.findall(
            r"^\| `(?:--|params\.)([a-z][a-z0-9_]*)` \|",
            reference,
            flags=re.MULTILINE,
        )
        documented = set(documented_names)
        schema_names = set(schema["properties"])

        self.assertEqual(documented, configured)
        self.assertEqual(schema_names, configured)
        self.assertEqual(len(documented_names), len(documented))
        self.assertEqual(schema["required"], ["input", "genome"])
        self.assertFalse(schema["additionalProperties"])
        for line in reference.splitlines():
            if re.match(
                r"^\| `(?:--|params\.)[a-z][a-z0-9_]*` \|",
                line,
            ):
                cells = line.split("|")
                self.assertTrue(cells[3].strip())
        self.assertIn(
            "[docs/parameters.md](docs/parameters.md)",
            readme,
        )

    def test_old_digenome_score_names_are_rejected(self) -> None:
        combined = "\n".join(
            (ROOT / path).read_text()
            for path in (
                "main.nf",
                "nextflow.config",
                "nextflow_schema.json",
                "docs/parameters.md",
            )
        )
        self.assertIn("digenome_pair_score_cutoff", combined)
        self.assertIn("--digenome-pair-score-cutoff", combined)
        self.assertNotIn("digenome_score_cutoff", combined)
        self.assertNotIn("--digenome-score-cutoff", combined)

    def test_chunk_plan_is_joined_before_one_to_many_expansion(self) -> None:
        main = (ROOT / "main.nf").read_text()
        request_block = main.split(
            "chunk_requests_ch = cleavage_requests_ch", 1
        )[1].split("CLEAVAGE_CALL_CHUNK", 1)[0]

        self.assertIn(
            ".join(PLAN_CLEAVAGE_CHUNKS.out.chunks, by: 0)",
            request_block,
        )
        self.assertIn(".flatMap", request_block)
        self.assertIn("files.sort { it.name }.collect", request_block)
        self.assertNotIn("planned_chunks_ch", main)

    def test_corrected_alignment_container_version(self) -> None:
        config = (ROOT / "nextflow.config").read_text()
        self.assertIn("bwa-mem2_v2.3_samtools_v1.22.sif", config)
        self.assertNotIn("bwa-mem2_v2.2.1", config)

    def test_alignment_retains_index_ready_dependency(self) -> None:
        main = (ROOT / "main.nf").read_text()
        pe_block = main.split(
            "process ALIGN_MARKDUP_PE", 1
        )[1].split("process ALIGN_MARKDUP_SE", 1)[0]
        se_block = main.split(
            "process ALIGN_MARKDUP_SE", 1
        )[1].split("process PLAN_CLEAVAGE_CHUNKS", 1)[0]

        self.assertIn("path(index_ready)", pe_block)
        self.assertIn("path(index_ready)", se_block)
        self.assertIn("test -s ${shellQuote(index_ready)}", pe_block)
        self.assertIn("test -s ${shellQuote(index_ready)}", se_block)
        self.assertIn(
            "tuple(values.index_prefix as String, ready)",
            main,
        )
        self.assertIn(
            "meta, read1, read2, index_prefix, index_ready ->",
            main,
        )
        self.assertIn(
            "meta, read1, index_prefix, index_ready ->",
            main,
        )

    def test_alignment_thread_budget_does_not_oversubscribe_pipe(self) -> None:
        main = (ROOT / "main.nf").read_text()
        self.assertIn(
            "total_cpus - piped_sort_threads - 1",
            main,
        )
        self.assertIn(
            "int samtools_threads = Math.max(0, total_cpus - 1)",
            main,
        )
        self.assertIn(
            "ALIGN_MARKDUP_PE requires at least 2 CPUs",
            main,
        )
        self.assertNotIn("| samtools view", main)
        self.assertNotIn("-t ${task.cpus}", main)

    def test_keep_multimappers_applies_coordinated_mapq_overrides(self) -> None:
        main = (ROOT / "main.nf").read_text()
        self.assertIn(
            "int effective_digenome_min_mapq = keep_multimappers ?",
            main,
        )
        self.assertIn(
            "int effective_ndigenome_min_mapq = keep_multimappers ?",
            main,
        )
        self.assertIn(
            "double effective_min_support_mean_mapq = keep_multimappers ?",
            main,
        )
        self.assertIn(
            "--digenome-min-mapq ${effective_digenome_min_mapq}",
            main,
        )
        self.assertIn(
            "--ndigenome-min-mapq ${effective_ndigenome_min_mapq}",
            main,
        )
        self.assertIn(
            "--min-support-mean-mapq ${effective_min_support_mean_mapq}",
            main,
        )

    def test_longleaf_profile_is_explicit(self) -> None:
        config = (ROOT / "nextflow.config").read_text()
        longleaf = (ROOT / "conf" / "longleaf.config").read_text()
        main = (ROOT / "main.nf").read_text()
        self.assertIn("longleaf {", config)
        for bind_path in ("/proj", "/work", "/users", "/overflow", "/nas"):
            self.assertIn(f"--bind {bind_path}", longleaf)
        self.assertNotIn(
            "/proj/jmsimon/Zylka/digenome-seq_pipeline",
            longleaf,
        )
        self.assertIn(
            'ref_cache = "${projectDir}/reference_cache"',
            longleaf,
        )
        self.assertIn(
            'cacheDir = "${projectDir}/containers"',
            longleaf,
        )
        self.assertIn("bash ${shellQuote(index_helper)}", main)

    def test_cleavage_definition_uses_apptainer_compatible_digest(self) -> None:
        definition = (
            ROOT / "containers" / "cleavage_pysam_v0.23.3.def"
        ).read_text()
        from_line = next(
            line for line in definition.splitlines() if line.startswith("From:")
        )
        self.assertIn("python@sha256:", from_line)
        self.assertNotIn("-bookworm@sha256:", from_line)
        self.assertIn("apt-get install -y --no-install-recommends procps", definition)

    def test_container_provenance_manifest_is_complete(self) -> None:
        manifest_path = ROOT / "containers" / "sources.tsv"
        with manifest_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))

        expected_fields = [
            "artifact",
            "source_tag",
            "oci_index_digest",
            "linux_amd64_manifest_digest",
            "construction",
            "sif_build_arch",
            "apptainer_version",
            "sif_build_time_utc",
            "registry_verified_utc",
        ]
        self.assertEqual(list(rows[0]), expected_fields)
        self.assertEqual(len(rows), 4)

        checksum_artifacts = {
            line.split()[1]
            for line in (
                ROOT / "containers" / "checksums.sha256"
            ).read_text().splitlines()
            if line.strip()
        }
        manifest_artifacts = {row["artifact"] for row in rows}
        self.assertEqual(manifest_artifacts, checksum_artifacts)

        digest_pattern = re.compile(r"sha256:[0-9a-f]{64}")
        timestamp_pattern = re.compile(
            r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z"
        )
        for row in rows:
            self.assertTrue(row["source_tag"].startswith("docker.io/"))
            self.assertRegex(row["oci_index_digest"], digest_pattern)
            self.assertRegex(
                row["linux_amd64_manifest_digest"],
                digest_pattern,
            )
            self.assertEqual(row["sif_build_arch"], "amd64")
            self.assertEqual(row["apptainer_version"], "1.5.0-1.el8")
            self.assertRegex(row["sif_build_time_utc"], timestamp_pattern)
            self.assertRegex(row["registry_verified_utc"], timestamp_pattern)

        cleavage_row = next(
            row
            for row in rows
            if row["artifact"] == "cleavage_pysam_v0.23.3.sif"
        )
        definition = (
            ROOT / "containers" / "cleavage_pysam_v0.23.3.def"
        ).read_text()
        self.assertIn(cleavage_row["oci_index_digest"], definition)
        self.assertEqual(
            cleavage_row["construction"],
            "definition_file:cleavage_pysam_v0.23.3.def",
        )

    def test_legacy_digenome_runtime_is_removed(self) -> None:
        sources = (ROOT / "containers" / "sources.tsv").read_text()
        main = (ROOT / "main.nf").read_text()
        self.assertNotIn("digenome_standalone_v1.0", sources)
        self.assertNotIn("digenome.output.txt", main)
        self.assertFalse((ROOT / "bin" / "call_ndigenome.py").exists())
        self.assertFalse((ROOT / "bin" / "compare_digenome_calls.py").exists())
        self.assertFalse((ROOT / "containers" / "build_ndigenome.sh").exists())

    def test_smoke_fixture_contains_unique_paired_alignments(self) -> None:
        fixture = ROOT / "tests" / "fixtures"
        samplesheet_header = (
            fixture / "tiny_samplesheet.csv"
        ).read_text().splitlines()[0]
        self.assertEqual(samplesheet_header, "sample,fastq_1,fastq_2")
        reference = "".join(
            line
            for line in (fixture / "tiny.fa").read_text().splitlines()
            if not line.startswith(">")
        )
        read1 = self.read_fastq(fixture / "tiny_R1.fastq")
        read2 = self.read_fastq(fixture / "tiny_R2.fastq")
        self.assertEqual(len(read1), 33)
        self.assertEqual(len(read1), len(read2))

        complement = str.maketrans("ACGT", "TGCA")
        forward_starts = Counter()
        reverse_endpoints = Counter()
        fragment_coordinates = set()
        for first, second in zip(read1, read2):
            first_name, first_sequence, first_quality = first
            second_name, second_sequence, second_quality = second
            self.assertEqual(first_name.removesuffix("/1"), second_name.removesuffix("/2"))
            self.assertEqual(len(first_sequence), len(first_quality))
            self.assertEqual(len(second_sequence), len(second_quality))
            self.assertEqual(reference.count(first_sequence), 1)
            reverse_complement = second_sequence.translate(complement)[::-1]
            self.assertEqual(reference.count(reverse_complement), 1)
            forward_start = reference.index(first_sequence)
            reverse_start = reference.index(reverse_complement)
            reverse_endpoint = reverse_start + len(reverse_complement) - 1
            forward_starts[forward_start] += 1
            reverse_endpoints[reverse_endpoint] += 1
            fragment_coordinates.add((forward_start, reverse_endpoint))

        self.assertEqual(len(fragment_coordinates), len(read1))
        self.assertEqual(forward_starts[250], 11)
        self.assertEqual(reverse_endpoints[250], 11)
        self.assertEqual(forward_starts[300], 11)

    def test_smoke_config_fits_a_two_cpu_executor(self) -> None:
        smoke_config = (
            ROOT / "tests" / "fixtures" / "smoke.config"
        ).read_text()
        self.assertIn("withName: PREPARE_BWAMEM2_INDEX", smoke_config)
        self.assertIn("cpus = 2", smoke_config)
        self.assertIn("withName: ALIGN_MARKDUP_PE", smoke_config)
        self.assertIn("withName: ALIGN_MARKDUP_SE", smoke_config)
        self.assertNotIn("|", smoke_config)


if __name__ == "__main__":
    unittest.main()
