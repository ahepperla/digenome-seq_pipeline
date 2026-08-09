#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
output="${script_dir}/cleavage_pysam_v0.23.3.sif"
definition="${script_dir}/cleavage_pysam_v0.23.3.def"

command -v apptainer >/dev/null || {
    echo "ERROR: apptainer is required to build $output" >&2
    exit 1
}

apptainer build --force "$output" "$definition"
reported_pysam=$(
    apptainer exec "$output" \
        python3 -c 'import importlib.metadata; print(importlib.metadata.version("pysam"))'
)
if [[ "$reported_pysam" != "0.23.3" ]]; then
    echo "ERROR: built image reports pysam $reported_pysam, expected 0.23.3" >&2
    exit 1
fi
apptainer exec "$output" sh -c 'command -v ps >/dev/null'

if command -v sha256sum >/dev/null; then
    checksum=$(sha256sum "$output" | awk '{print $1}')
elif command -v shasum >/dev/null; then
    checksum=$(shasum -a 256 "$output" | awk '{print $1}')
else
    echo "ERROR: sha256sum or shasum is required to checksum $output" >&2
    exit 1
fi
manifest="${script_dir}/checksums.sha256"
temporary=$(mktemp "${manifest}.XXXXXX")
awk \
    '$2 != "cleavage_pysam_v0.23.3.sif"' \
    "$manifest" > "$temporary"
printf '%s  %s\n' "$checksum" "$(basename "$output")" >> "$temporary"
mv "$temporary" "$manifest"
printf '%s  %s\n' "$checksum" "$(basename "$output")"
