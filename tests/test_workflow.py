"""Checks on the Nextflow workflow, its configuration, and its inputs.

The lint and stub-run tests use $NEXTFLOW if set (for example
NEXTFLOW=nextflow-25.04.7 to match Longleaf), else `nextflow` on PATH, and are
skipped when neither exists. Stub runs set NXF_SYNTAX_PARSER=v1 so that newer
Nextflow parses command-line parameters the way Longleaf's 25.04 does.
"""

from __future__ import annotations

import ast
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from collections import Counter
from dataclasses import fields
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests" / "fixtures"))
sys.path.insert(0, str(ROOT / "bin"))
import build_smoke_fixture  # noqa: E402
from cleavage.settings import CallerSettings  # noqa: E402

NEXTFLOW = shutil.which(os.environ.get("NEXTFLOW", "nextflow"))


def config_parameters() -> dict[str, str]:
    """Parameter name -> default as written in nextflow.config's params block."""
    block = (ROOT / "nextflow.config").read_text().split("\nparams {", 1)[1].split("\n}", 1)[0]
    return dict(re.findall(r"^    ([a-z][a-z0-9_]*) = (.+)$", block, flags=re.MULTILINE))


def schema_parameters() -> dict[str, dict]:
    schema = json.loads((ROOT / "nextflow_schema.json").read_text())
    return {
        name: definition
        for group in schema["$defs"].values()
        for name, definition in group["properties"].items()
    }


def documented_parameters() -> list[str]:
    text = (ROOT / "docs" / "parameters.md").read_text()
    return re.findall(r"^\| `--([a-z][a-z0-9_]*)` \|", text, flags=re.MULTILINE)


def config_literal(text: str):
    """A config default as a Python value, or None for expressions."""
    text = {"true": "True", "false": "False", "null": "None"}.get(text, text)
    if "$" in text or text.startswith("["):
        return None
    return ast.literal_eval(text)


class ParameterTests(unittest.TestCase):
    def test_config_schema_and_docs_list_the_same_parameters(self) -> None:
        configured = set(config_parameters()) - {"genomes"}
        documented = documented_parameters()
        self.assertEqual(set(schema_parameters()), configured)
        self.assertEqual(set(documented), configured)
        self.assertEqual(len(documented), len(set(documented)))

    def test_schema_defaults_match_config(self) -> None:
        schema = schema_parameters()
        for name, text in config_parameters().items():
            value = config_literal(text)
            if name in ("genomes", "ref_cache") or value is None:
                continue
            with self.subTest(parameter=name):
                self.assertEqual(schema[name].get("default"), value)

    def test_every_schema_parameter_is_described(self) -> None:
        for name, definition in schema_parameters().items():
            with self.subTest(parameter=name):
                self.assertTrue(definition["description"].strip())


class CallerSettingsTests(unittest.TestCase):
    def test_main_nf_builds_exactly_the_caller_settings(self) -> None:
        """callerSettings() in main.nf and CallerSettings must list the same names."""
        block = (ROOT / "main.nf").read_text().split("def callerSettings()", 1)[1].split("\n}\n", 1)[0]
        keys = re.findall(r"^        ([a-z_]+):", block, flags=re.MULTILINE)
        self.assertEqual(sorted(keys), sorted(field.name for field in fields(CallerSettings)))


class ContainerTests(unittest.TestCase):
    def test_provenance_records_are_complete(self) -> None:
        with (ROOT / "containers" / "sources.tsv").open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(
            list(rows[0]),
            [
                "artifact", "source_tag", "oci_index_digest", "linux_amd64_manifest_digest",
                "construction", "sif_build_arch", "apptainer_version", "sif_build_time_utc",
                "registry_verified_utc",
            ],
        )
        checksummed = {
            line.split()[1]
            for line in (ROOT / "containers" / "checksums.sha256").read_text().splitlines()
            if line.strip()
        }
        self.assertEqual({row["artifact"] for row in rows}, checksummed)
        for row in rows:
            with self.subTest(artifact=row["artifact"]):
                self.assertRegex(row["oci_index_digest"], r"^sha256:[0-9a-f]{64}$")
                self.assertRegex(row["linux_amd64_manifest_digest"], r"^sha256:[0-9a-f]{64}$")

    def test_cleavage_image_is_built_from_the_recorded_digest(self) -> None:
        definition = (ROOT / "containers" / "cleavage_pysam_v0.23.3.def").read_text()
        with (ROOT / "containers" / "sources.tsv").open(newline="") as handle:
            row = next(row for row in csv.DictReader(handle, delimiter="\t") if row["artifact"] == "cleavage_pysam_v0.23.3.sif")
        self.assertIn(f"From: python@{row['oci_index_digest']}", definition)
        self.assertIn("pysam==0.23.3", definition)

    def test_every_configured_image_has_provenance(self) -> None:
        configured = set(re.findall(r"containers/([\w.-]+\.sif)", (ROOT / "conf" / "base.config").read_text()))
        checksummed = {
            line.split()[1]
            for line in (ROOT / "containers" / "checksums.sha256").read_text().splitlines()
            if line.strip()
        }
        self.assertEqual(configured, checksummed)


class SmokeFixtureTests(unittest.TestCase):
    def test_fixture_reads_align_uniquely_with_the_expected_pileups(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            samplesheet = build_smoke_fixture.build(directory)
            self.assertEqual(
                samplesheet.read_text().splitlines()[0],
                "sample,fastq_1,fastq_2",
            )
            read1 = self.sequences(directory / "tiny_R1.fastq.gz")
            read2 = self.sequences(directory / "tiny_R2.fastq.gz")
        reference = build_smoke_fixture.reference()
        forward_starts: Counter[int] = Counter()
        reverse_endpoints: Counter[int] = Counter()
        fragments = set()
        for first, second in zip(read1, read2):
            mate = build_smoke_fixture.reverse_complement(second)
            self.assertEqual(reference.count(first), 1)
            self.assertEqual(reference.count(mate), 1)
            start = reference.index(first)
            endpoint = reference.index(mate) + len(mate) - 1
            forward_starts[start] += 1
            reverse_endpoints[endpoint] += 1
            fragments.add((start, endpoint))
        self.assertEqual(len(read1), 33)
        self.assertEqual(len(fragments), 33)
        self.assertEqual(forward_starts[250], 11)
        self.assertEqual(reverse_endpoints[250], 11)
        self.assertEqual(forward_starts[300], 11)

    @staticmethod
    def sequences(path: Path) -> list[str]:
        with gzip.open(path, "rt") as handle:
            return handle.read().splitlines()[1::4]


@unittest.skipIf(NEXTFLOW is None, "nextflow is not installed")
class NextflowTests(unittest.TestCase):
    def test_lint_finds_no_errors(self) -> None:
        result = subprocess.run(
            [NEXTFLOW, "lint", "main.nf", "nextflow.config", "conf"],
            cwd=ROOT, capture_output=True, text=True, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertNotIn("had errors", result.stdout)

    SHORT_READ_PROCESSES = {
        "SAMPLESHEET", "SHORT_READS:PREPARE_INDEX", "SHORT_READS:FASTP", "SHORT_READS:ALIGN",
        "CALL_CHUNK", "FINALIZE", "MULTIQC",
    }

    def test_digenome_stub_run_publishes_every_output(self) -> None:
        samples = ["Treated", "Control", "Uncontrolled", "SingleEnd"]
        called = ["Treated", "Uncontrolled", "SingleEnd"]
        out, info = self.stub_run("digenome", self.fastq_sheet, samples, ["--keep_multimappers"])
        self.assertEqual(self.published(out), self.expected_files("digenome", samples, called))
        self.assertEqual(self.processes(out), self.SHORT_READ_PROCESSES)
        self.assertEqual(self.chunk_tasks(out), Counter({sample: 3 for sample in called}))
        # --keep_multimappers lowers every MAPQ threshold together.
        settings = info["caller_settings"]
        self.assertEqual(
            (settings["digenome_min_mapq"], settings["ndigenome_min_mapq"], settings["cleavage_min_support_mean_mapq"]),
            (0, 0, 0),
        )

    def test_ndigenome_stub_run_publishes_every_output(self) -> None:
        samples = ["Treated", "Control", "Uncontrolled"]
        called = ["Treated", "Uncontrolled"]
        out, info = self.stub_run("ndigenome", self.fastq_sheet, samples, [])
        self.assertEqual(self.published(out), self.expected_files("ndigenome", samples, called))
        self.assertEqual(self.processes(out), self.SHORT_READ_PROCESSES)
        self.assertEqual(self.chunk_tasks(out), Counter({sample: 3 for sample in called}))
        self.assertEqual((info["genome"], Path(info["fasta"]).name), ("tiny", "tiny.fa"))
        settings = info["caller_settings"]
        self.assertEqual((settings["ndigenome_min_mapq"], settings["cleavage_min_support_mean_mapq"]), (1, 10))

    def test_long_read_stub_run_calls_the_given_bams(self) -> None:
        samples = ["Treated", "Control", "Uncontrolled"]
        called = ["Treated", "Uncontrolled"]
        # No --genome: long reads need no reference.
        out, info = self.stub_run("ndigenome", self.bam_sheet, samples, ["--long_reads", "-params-file", "no_genome.json"])
        aligned_outputs = ("bam/", "fastp/", "qc/")
        self.assertEqual(
            self.published(out),
            {path for path in self.expected_files("ndigenome", samples, called) if not path.startswith(aligned_outputs)},
        )
        self.assertEqual(self.processes(out), {"SAMPLESHEET", "CALL_CHUNK", "FINALIZE", "MULTIQC"})
        self.assertEqual(self.chunk_tasks(out), Counter({sample: 3 for sample in called}))
        self.assertEqual((info["genome"], info["fasta"], info["ref_cache"]), (None, None, None))
        # --long_reads loosens the clip and indel limits together.
        settings = info["caller_settings"]
        self.assertEqual(
            (settings["long_reads"], settings["cleavage_max_softclip_fraction"], settings["cleavage_max_indel_fraction"]),
            (True, 1.0, 1.0),
        )

    def stub_run(self, analysis: str, write_sheet, samples: list[str], options: list[str]) -> tuple[Path, dict]:
        """Run the stub workflow on empty inputs from `write_sheet`; Treated
        uses Control. Returns the output directory and analysis_parameters.json."""
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        directory = Path(temporary.name)
        (directory / "samplesheet.csv").write_text(write_sheet(directory, samples))
        (directory / "no_genome.json").write_text('{"genome": null}\n')
        out = directory / "results"
        environment = {
            **os.environ,
            "NXF_SYNTAX_PARSER": "v1",
            "PATH": f"{Path(sys.executable).parent}:{os.environ['PATH']}",
        }
        result = subprocess.run(
            [
                NEXTFLOW, "run", str(ROOT), "-profile", "test", "-stub-run",
                "--input", str(directory / "samplesheet.csv"),
                "--analysis", analysis,
                "--outdir", str(out),
                "-work-dir", str(directory / "work"),
                *options,
            ],
            cwd=directory, env=environment, capture_output=True, text=True, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout[-3000:] + result.stderr[-3000:])
        return out, json.loads((out / "pipeline_info" / "analysis_parameters.json").read_text())

    @staticmethod
    def fastq_sheet(directory: Path, samples: list[str]) -> str:
        """Empty FASTQs. Treated has two lanes; SingleEnd is single-end."""
        rows = ["sample,fastq_1,fastq_2,control"]
        for sample in samples:
            lanes = 2 if sample == "Treated" else 1
            control = "Control" if sample == "Treated" else ""
            for lane in range(lanes):
                read1 = directory / f"{sample}_L{lane}_R1.fastq.gz"
                read2 = directory / f"{sample}_L{lane}_R2.fastq.gz"
                read1.touch()
                if sample == "SingleEnd":
                    rows.append(f"{sample},{read1},,{control}")
                else:
                    read2.touch()
                    rows.append(f"{sample},{read1},{read2},{control}")
        return "\n".join(rows) + "\n"

    @staticmethod
    def bam_sheet(directory: Path, samples: list[str]) -> str:
        """Empty indexed BAMs, all named reads.bam, so the treated and control
        BAMs of one task share a file name."""
        rows = ["sample,bam,control"]
        for sample in samples:
            bam = directory / sample / "reads.bam"
            bam.parent.mkdir()
            bam.touch()
            Path(f"{bam}.bai").touch()
            rows.append(f"{sample},{bam},{'Control' if sample == 'Treated' else ''}")
        return "\n".join(rows) + "\n"

    @staticmethod
    def processes(out: Path) -> set[str]:
        with (out / "pipeline_info" / "trace.txt").open(newline="") as handle:
            return {row["name"].split(" (")[0] for row in csv.DictReader(handle, delimiter="\t")}

    @staticmethod
    def published(out: Path) -> set[str]:
        return {str(path.relative_to(out)) for path in out.rglob("*") if path.is_file()}

    @staticmethod
    def expected_files(analysis: str, samples: list[str], called: list[str]) -> set[str]:
        files = {
            "multiqc/multiqc_report.html",
            "pipeline_info/analysis_parameters.json",
            "pipeline_info/container_checksums.sha256",
            "pipeline_info/container_sources.tsv",
            "pipeline_info/dag.html",
            "pipeline_info/nextflow_report.html",
            "pipeline_info/timeline.html",
            "pipeline_info/trace.txt",
        }
        for sample in samples:
            files |= {
                f"bam/{sample}.sorted.markdup.bam",
                f"bam/{sample}.sorted.markdup.bam.bai",
                f"fastp/{sample}.fastp.html",
                f"fastp/{sample}.fastp.json",
                f"qc/{sample}.flagstat.txt",
                f"qc/{sample}.stats.txt",
                f"qc/{sample}.markdup.metrics.txt",
            }
        for sample in called:
            stem = f"{analysis}/{sample}.{analysis}"
            files |= {
                f"{stem}.all.tsv",
                f"{stem}.high_confidence.tsv",
                f"{stem}.manual_review.tsv",
                f"{stem}.artifact.tsv",
                f"{stem}.bed",
                f"{stem}.qc.json",
                f"{stem}_mqc.tsv",
                f"pipeline_info/cleavage_chunks/{sample}.cleavage_chunks.tsv",
            }
        return files

    @staticmethod
    def chunk_tasks(out: Path) -> Counter[str]:
        with (out / "pipeline_info" / "trace.txt").open(newline="") as handle:
            names = [row["name"] for row in csv.DictReader(handle, delimiter="\t")]
        return Counter(
            re.match(r"CALL_CHUNK \((\w+):\d+\)", name).group(1)
            for name in names
            if name.startswith("CALL_CHUNK")
        )


if __name__ == "__main__":
    unittest.main()
