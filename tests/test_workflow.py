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

import pysam

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


def unique_keys(pairs: list[tuple[str, object]]) -> dict:
    """A JSON object hook that rejects duplicate keys, as nf-schema does."""
    keys = [key for key, _value in pairs]
    duplicates = sorted({key for key in keys if keys.count(key) > 1})
    if duplicates:
        raise ValueError(f"duplicate JSON keys: {', '.join(duplicates)}")
    return dict(pairs)


def schema_parameters() -> dict[str, dict]:
    schema = json.loads((ROOT / "nextflow_schema.json").read_text(), object_pairs_hook=unique_keys)
    return {
        name: definition
        for group in schema["$defs"].values()
        for name, definition in group["properties"].items()
    }


def documented_parameters() -> list[tuple[str, str]]:
    """(name, default cell) for each row of the README's parameter tables."""
    text = (ROOT / "README.md").read_text()
    return re.findall(r"^\| `--([a-z][a-z0-9_]*)` \| ([^|]*) \|", text, flags=re.MULTILINE)


def config_literal(text: str):
    """A config default as a Python value, or None for expressions."""
    text = {"true": "True", "false": "False", "null": "None"}.get(text, text)
    if "$" in text or text.startswith("["):
        return None
    return ast.literal_eval(text)


class ParameterTests(unittest.TestCase):
    def test_config_schema_and_docs_list_the_same_parameters(self) -> None:
        configured = set(config_parameters()) - {"genomes"}
        documented = [name for name, _default in documented_parameters()]
        self.assertEqual(set(schema_parameters()), configured)
        self.assertEqual(set(documented), configured)
        self.assertEqual(len(documented), len(set(documented)))

    def test_documented_defaults_match_config(self) -> None:
        """Numbers and booleans in the README's Default column; words such as
        "required" or "none" are not checked."""
        config = config_parameters()
        for name, cell in documented_parameters():
            try:
                documented = json.loads(cell.strip().strip("`"))
            except ValueError:
                continue
            with self.subTest(parameter=name):
                self.assertEqual(documented, config_literal(config[name]))

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

    INPUT_MTIME = 1_000_000_000  # set on input BAMs, to show the run never writes them

    FASTQ_PROCESSES = {
        "SAMPLESHEET", "FROM_FASTQS:PREPARE_INDEX", "FROM_FASTQS:FASTP", "FROM_FASTQS:ALIGN",
        "CALL_CHUNK", "FINALIZE", "MULTIQC",
    }
    BAM_PROCESSES = {"SAMPLESHEET", "FROM_BAMS:SORT_BAM", "CALL_CHUNK", "FINALIZE", "MULTIQC"}
    ALIGNED_OUTPUTS = ("bam/", "fastp/", "qc/")

    def test_digenome_stub_run_publishes_every_output(self) -> None:
        samples = ["Treated", "Control", "Uncontrolled", "SingleEnd"]
        called = ["Treated", "Uncontrolled", "SingleEnd"]
        # SMOKE is an alias of the test profile's genome, in another case.
        out, info = self.stub_run("digenome", self.fastq_sheet, samples, ["--keep_multimappers", "--genome", "SMOKE"])
        self.assertEqual(self.published(out), self.expected_files("digenome", samples, called))
        self.assertEqual(self.processes(out), self.FASTQ_PROCESSES)
        self.assertEqual(self.chunk_tasks(out), Counter({sample: 3 for sample in called}))
        # The index and the run record use the genome's name, not the alias.
        self.assertIn("FROM_FASTQS:PREPARE_INDEX (tiny)", (out / "pipeline_info" / "trace.txt").read_text())
        self.assertEqual(info["genome"], "tiny")
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
        self.assertEqual(self.processes(out), self.FASTQ_PROCESSES)
        self.assertEqual(self.chunk_tasks(out), Counter({sample: 3 for sample in called}))
        self.assertEqual((info["genome"], Path(info["fasta"]).name), ("tiny", "tiny.fa"))
        settings = info["caller_settings"]
        self.assertEqual((settings["ndigenome_min_mapq"], settings["cleavage_min_support_mean_mapq"]), (1, 10))

    def test_long_read_stub_run_calls_the_given_bams(self) -> None:
        samples = ["Treated", "Control", "Uncontrolled"]
        called = ["Treated", "Uncontrolled"]
        # No --genome: long reads need no reference.
        out, info = self.stub_run("ndigenome", self.bam_sheet, samples, ["--long_reads", "-params-file", "no_genome.json"])
        self.assertEqual(
            self.published(out),
            {path for path in self.expected_files("ndigenome", samples, called) if not path.startswith(self.ALIGNED_OUTPUTS)},
        )
        self.assertEqual(self.processes(out), self.BAM_PROCESSES)
        # Only the unsorted BAM is sorted, and the input file is left untouched.
        trace = (out / "pipeline_info" / "trace.txt").read_text()
        self.assertEqual(re.findall(r"FROM_BAMS:SORT_BAM \((\w+)\)", trace), ["Uncontrolled"])
        unsorted = out.parent / "Uncontrolled" / "Uncontrolled.sorted.bam"
        self.assertEqual(unsorted.stat().st_mtime, self.INPUT_MTIME)
        self.assertEqual(self.chunk_tasks(out), Counter({sample: 3 for sample in called}))
        self.assertEqual((info["genome"], info["fasta"], info["ref_cache"]), (None, None, None))
        # --long_reads loosens the clip and indel limits together.
        settings = info["caller_settings"]
        self.assertEqual(
            (settings["long_reads"], settings["cleavage_max_softclip_fraction"], settings["cleavage_max_indel_fraction"]),
            (True, 1.0, 1.0),
        )

    def test_short_read_bam_stub_run_skips_alignment(self) -> None:
        samples = ["Treated", "Control", "Uncontrolled"]
        called = ["Treated", "Uncontrolled"]
        out, info = self.stub_run("digenome", self.bam_sheet, samples, ["-params-file", "no_genome.json"])
        self.assertEqual(
            self.published(out),
            {path for path in self.expected_files("digenome", samples, called) if not path.startswith(self.ALIGNED_OUTPUTS)},
        )
        self.assertEqual(self.processes(out), self.BAM_PROCESSES)
        self.assertEqual(self.chunk_tasks(out), Counter({sample: 3 for sample in called}))
        self.assertEqual((info["genome"], info["fasta"]), (None, None))
        # Short-read rules: the long-read limits stay at their defaults.
        settings = info["caller_settings"]
        self.assertEqual(
            (settings["long_reads"], settings["cleavage_max_softclip_fraction"], settings["cleavage_max_indel_fraction"]),
            (False, 0.2, 0.2),
        )

    def test_long_reads_need_a_bam_samplesheet(self) -> None:
        result = self.launch("digenome", self.fastq_sheet, ["Treated", "Control"], ["--long_reads"])[1]
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--long_reads needs a samplesheet of aligned BAMs", result.stdout + result.stderr)

    def stub_run(self, analysis: str, write_sheet, samples: list[str], options: list[str]) -> tuple[Path, dict]:
        """Run the stub workflow on empty inputs from `write_sheet`; Treated
        uses Control. Returns the output directory and analysis_parameters.json."""
        out, result = self.launch(analysis, write_sheet, samples, options)
        self.assertEqual(result.returncode, 0, result.stdout[-3000:] + result.stderr[-3000:])
        return out, json.loads((out / "pipeline_info" / "analysis_parameters.json").read_text())

    def launch(self, analysis: str, write_sheet, samples: list[str], options: list[str]):
        """Write the inputs and start a stub run; returns the output directory
        and the finished process."""
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
        return out, result

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
        """Header-only BAMs. Treated and Control are sorted, indexed, and both
        named reads.bam, so the treated and control BAMs of one task share a
        file name. Uncontrolled is unsorted and already named like SORT_BAM's
        output, which must not overwrite it."""
        rows = ["sample,bam,control"]
        for sample in samples:
            unsorted = sample == "Uncontrolled"
            bam = directory / sample / (f"{sample}.sorted.bam" if unsorted else "reads.bam")
            bam.parent.mkdir()
            header = {"HD": {"VN": "1.6", "SO": "unsorted" if unsorted else "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
            with pysam.AlignmentFile(str(bam), "wb", header=header):
                pass
            if not unsorted:
                Path(f"{bam}.bai").touch()
            os.utime(bam, (NextflowTests.INPUT_MTIME, NextflowTests.INPUT_MTIME))
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
