from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
from validate_pipeline_params import validate_preflight  # noqa: E402


class PreflightTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tempdir.name)
        self.schema = json.loads(
            (ROOT / "nextflow_schema.json").read_text()
        )
        self.parameters = {
            name: definition.get("default")
            for name, definition in self.schema["properties"].items()
        }
        samplesheet = self.tmp / "samplesheet.csv"
        samplesheet.write_text("sample,fastq_1\n")
        fasta = self.tmp / "reference.fa"
        fasta.write_text(">chr1\nA\n")
        containers = {}
        for name in ("python", "fastp", "align", "cleavage", "multiqc"):
            path = self.tmp / f"{name}.sif"
            path.write_text("container\n")
            containers[name] = str(path)
        self.parameters.update(
            {
                "input": str(samplesheet),
                "genome": "tiny",
                "genomes": {
                    "tiny": {
                        "fasta": str(fasta),
                        "aliases": [],
                    }
                },
                "containers": containers,
                "container_bind_paths": [str(self.tmp)],
                "ref_cache": str(self.tmp / "reference_cache"),
                "outdir": str(self.tmp / "results"),
            }
        )
        self.document = {
            "parameters": self.parameters,
            "paths": {
                "input": str(samplesheet),
                "fasta": str(fasta),
                "genome_blacklist": "",
                "ref_cache": str(self.tmp / "reference_cache"),
                "outdir": str(self.tmp / "results"),
            },
        }

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_valid_configuration_passes(self) -> None:
        self.assertEqual(
            validate_preflight(self.document, self.schema),
            [],
        )

    def test_unknown_parameter_is_rejected(self) -> None:
        self.parameters["typo_parameter"] = 1
        with self.assertRaisesRegex(ValueError, "Unknown pipeline"):
            validate_preflight(self.document, self.schema)

    def test_fraction_outside_range_is_rejected(self) -> None:
        self.parameters["cleavage_control_max_q"] = 1.5
        with self.assertRaisesRegex(ValueError, "must be at most 1"):
            validate_preflight(self.document, self.schema)

    def test_duplicate_genome_alias_is_rejected(self) -> None:
        self.parameters["genomes"]["other"] = {
            "fasta": self.document["paths"]["fasta"],
            "aliases": ["tiny"],
        }
        with self.assertRaisesRegex(ValueError, "shared by tiny and other"):
            validate_preflight(self.document, self.schema)

    def test_empty_genome_map_is_rejected(self) -> None:
        self.parameters["genomes"] = {}
        with self.assertRaisesRegex(
            ValueError,
            "genomes must contain at least 1 entry",
        ):
            validate_preflight(self.document, self.schema)

    def test_missing_local_container_is_rejected(self) -> None:
        self.parameters["containers"]["cleavage"] = str(
            self.tmp / "missing.sif"
        )
        with self.assertRaisesRegex(ValueError, "Container 'cleavage'"):
            validate_preflight(self.document, self.schema)

    def test_missing_resolved_path_is_rejected(self) -> None:
        del self.document["paths"]["fasta"]
        with self.assertRaisesRegex(
            ValueError,
            "missing resolved path.*fasta",
        ):
            validate_preflight(self.document, self.schema)


if __name__ == "__main__":
    unittest.main()
