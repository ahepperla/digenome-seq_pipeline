#!/usr/bin/env python3
"""Validate pipeline parameters and required host paths before execution."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path
from typing import Any


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate pipeline parameters and runtime paths."
    )
    parser.add_argument("--parameters", required=True)
    parser.add_argument("--schema", required=True)
    parser.add_argument("--ready", required=True)
    return parser


def _matches_type(value: Any, expected: str) -> bool:
    if expected == "null":
        return value is None
    if expected == "boolean":
        return isinstance(value, bool)
    if expected == "integer":
        return isinstance(value, int) and not isinstance(value, bool)
    if expected == "number":
        return (
            isinstance(value, (int, float))
            and not isinstance(value, bool)
        )
    if expected == "string":
        return isinstance(value, str)
    if expected == "array":
        return isinstance(value, list)
    if expected == "object":
        return isinstance(value, dict)
    raise ValueError(f"Unsupported schema type: {expected}")


def _validate_schema_value(
    name: str,
    value: Any,
    definition: dict[str, Any],
) -> None:
    expected_types = definition.get("type")
    if isinstance(expected_types, str):
        expected_types = [expected_types]
    if expected_types and not any(
        _matches_type(value, expected)
        for expected in expected_types
    ):
        raise ValueError(
            f"--{name} has type {type(value).__name__}; expected "
            + " or ".join(expected_types)
        )
    if value is None:
        return
    if "enum" in definition and value not in definition["enum"]:
        choices = ", ".join(str(choice) for choice in definition["enum"])
        raise ValueError(f"--{name} must be one of: {choices}")
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if "minimum" in definition and value < definition["minimum"]:
            raise ValueError(
                f"--{name} must be at least {definition['minimum']}"
            )
        if "maximum" in definition and value > definition["maximum"]:
            raise ValueError(
                f"--{name} must be at most {definition['maximum']}"
            )
    if isinstance(value, str):
        if len(value) < int(definition.get("minLength", 0)):
            raise ValueError(f"--{name} must not be empty")
    if isinstance(value, list):
        item_definition = definition.get("items")
        if item_definition:
            for index, item in enumerate(value):
                _validate_schema_value(
                    f"{name}[{index}]",
                    item,
                    item_definition,
                )
    if isinstance(value, dict):
        if (
            "minProperties" in definition
            and len(value) < definition["minProperties"]
        ):
            raise ValueError(
                f"--{name} must contain at least "
                f"{definition['minProperties']} entr"
                + ("y" if definition["minProperties"] == 1 else "ies")
            )
        required = set(definition.get("required", []))
        missing = sorted(required - set(value))
        if missing:
            raise ValueError(
                f"--{name} is missing required key(s): "
                + ", ".join(missing)
            )
        properties = definition.get("properties", {})
        if definition.get("additionalProperties") is False:
            unknown = sorted(set(value) - set(properties))
            if unknown:
                raise ValueError(
                    f"--{name} has unknown key(s): "
                    + ", ".join(unknown)
                )
        additional = definition.get("additionalProperties")
        for key, item in value.items():
            item_definition = properties.get(key)
            if item_definition is None and isinstance(additional, dict):
                item_definition = additional
            if item_definition is not None:
                _validate_schema_value(
                    f"{name}.{key}",
                    item,
                    item_definition,
                )


def _require_readable_file(label: str, value: str) -> Path:
    path = Path(value)
    if not path.is_file():
        raise ValueError(f"{label} does not exist or is not a file: {value}")
    if not os.access(path, os.R_OK):
        raise ValueError(f"{label} is not readable: {value}")
    return path


def _nearest_existing_parent(path: Path) -> Path:
    candidate = path
    while not candidate.exists():
        parent = candidate.parent
        if parent == candidate:
            break
        candidate = parent
    return candidate


def _require_creatable_path(label: str, value: str) -> list[str]:
    path = Path(value)
    warnings: list[str] = []
    if path.exists():
        if not path.is_dir():
            raise ValueError(f"{label} is not a directory: {value}")
        if not os.access(path, os.R_OK):
            raise ValueError(f"{label} is not readable: {value}")
        if not os.access(path, os.W_OK):
            warnings.append(
                f"{label} is read-only; existing content can be reused but "
                "new files cannot be created"
            )
        return warnings
    parent = _nearest_existing_parent(path.parent)
    if not parent.is_dir() or not os.access(parent, os.W_OK):
        raise ValueError(
            f"{label} cannot be created because its nearest existing parent "
            f"is not writable: {parent}"
        )
    return warnings


def _validate_genome_aliases(genomes: dict[str, Any]) -> None:
    owners: dict[str, str] = {}
    for genome, config in genomes.items():
        names = [genome, *config.get("aliases", [])]
        for name in names:
            normalized = str(name).lower()
            previous = owners.get(normalized)
            if previous is not None and previous != genome:
                raise ValueError(
                    f"Genome name or alias '{name}' is shared by "
                    f"{previous} and {genome}"
                )
            owners[normalized] = genome


def _looks_like_local_container(value: str) -> bool:
    return (
        "://" not in value
        and (
            value.startswith(("/", "./", "../"))
            or value.endswith(".sif")
        )
    )


def validate_preflight(
    document: dict[str, Any],
    schema: dict[str, Any],
) -> list[str]:
    parameters = document.get("parameters")
    paths = document.get("paths")
    if not isinstance(parameters, dict) or not isinstance(paths, dict):
        raise ValueError(
            "Preflight document must contain parameters and paths objects"
        )
    required_paths = {
        "input",
        "fasta",
        "genome_blacklist",
        "outdir",
        "ref_cache",
    }
    missing_paths = sorted(required_paths - set(paths))
    if missing_paths:
        raise ValueError(
            "Preflight document is missing resolved path(s): "
            + ", ".join(missing_paths)
        )
    for name in required_paths:
        if not isinstance(paths[name], str):
            raise ValueError(f"Resolved path '{name}' must be a string")

    definitions = schema.get("properties", {})
    if schema.get("additionalProperties") is False:
        unknown = sorted(set(parameters) - set(definitions))
        if unknown:
            raise ValueError(
                "Unknown pipeline parameter(s): "
                + ", ".join(f"--{name}" for name in unknown)
            )
    missing = sorted(set(schema.get("required", [])) - set(parameters))
    if missing:
        raise ValueError(
            "Missing required pipeline parameter(s): "
            + ", ".join(f"--{name}" for name in missing)
        )
    for name, value in parameters.items():
        definition = definitions.get(name)
        if definition is not None:
            _validate_schema_value(name, value, definition)

    for required_name in ("input", "genome"):
        value = parameters.get(required_name)
        if value is None or not str(value).strip():
            raise ValueError(f"--{required_name} is required")

    _validate_genome_aliases(parameters["genomes"])
    _require_readable_file("Input samplesheet", paths["input"])
    _require_readable_file("Reference FASTA", paths["fasta"])
    blacklist = paths.get("genome_blacklist")
    if blacklist:
        _require_readable_file("Genome blacklist", blacklist)

    for name, value in parameters["containers"].items():
        text = str(value)
        if _looks_like_local_container(text):
            _require_readable_file(f"Container '{name}'", text)
    for bind_path in parameters["container_bind_paths"]:
        path = Path(bind_path)
        if not path.is_dir():
            raise ValueError(
                f"Container bind path does not exist: {bind_path}"
            )
        if not os.access(path, os.R_OK):
            raise ValueError(
                f"Container bind path is not readable: {bind_path}"
            )

    warnings = []
    warnings.extend(
        _require_creatable_path("Output directory", paths["outdir"])
    )
    warnings.extend(
        _require_creatable_path("Reference cache", paths["ref_cache"])
    )
    return warnings


def main() -> None:
    args = build_parser().parse_args()
    try:
        parameter_path = Path(args.parameters)
        schema_path = Path(args.schema)
        document = json.loads(parameter_path.read_text())
        schema = json.loads(schema_path.read_text())
        warnings = validate_preflight(document, schema)
        ready = {
            "schema_version": 1,
            "status": "PASS",
            "parameters_sha256": hashlib.sha256(
                parameter_path.read_bytes()
            ).hexdigest(),
            "resolved_paths": document["paths"],
            "schema_sha256": hashlib.sha256(
                schema_path.read_bytes()
            ).hexdigest(),
            "warnings": warnings,
        }
        Path(args.ready).write_text(
            json.dumps(ready, indent=2, sort_keys=True) + "\n"
        )
        for warning in warnings:
            print(f"WARNING: {warning}", file=sys.stderr)
        print("Preflight validation passed.")
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
