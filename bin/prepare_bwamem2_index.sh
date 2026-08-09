#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: prepare_bwamem2_index.sh \
  --genome NAME --fasta FASTA --cache-dir DIR --ready FILE \
  [--lock-timeout-seconds 172800] [--stale-lock-seconds 172800]
EOF
    exit 2
}

genome=''
fasta=''
cache_dir=''
ready_file=''
lock_timeout_seconds=172800
stale_lock_seconds=172800

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genome) genome=$2; shift 2 ;;
        --fasta) fasta=$2; shift 2 ;;
        --cache-dir) cache_dir=$2; shift 2 ;;
        --ready) ready_file=$2; shift 2 ;;
        --lock-timeout-seconds) lock_timeout_seconds=$2; shift 2 ;;
        --stale-lock-seconds) stale_lock_seconds=$2; shift 2 ;;
        *) usage ;;
    esac
done

[[ -n "$genome" && -n "$fasta" && -n "$cache_dir" && -n "$ready_file" ]] || usage
[[ -s "$fasta" ]] || { echo "ERROR: FASTA does not exist or is empty: $fasta" >&2; exit 1; }
command -v bwa-mem2 >/dev/null || { echo "ERROR: bwa-mem2 is not available" >&2; exit 1; }

sha256_file() {
    if command -v sha256sum >/dev/null; then
        sha256sum "$1" | awk '{print $1}'
    elif command -v shasum >/dev/null; then
        shasum -a 256 "$1" | awk '{print $1}'
    else
        echo "ERROR: sha256sum or shasum is required" >&2
        exit 1
    fi
}

file_size() {
    if stat -c %s "$1" >/dev/null 2>&1; then
        stat -c %s "$1"
    else
        stat -f %z "$1"
    fi
}

file_mtime() {
    if stat -c %Y "$1" >/dev/null 2>&1; then
        stat -c %Y "$1"
    else
        stat -f %m "$1"
    fi
}

index_complete() {
    local prefix=$1
    local suffix
    for suffix in .0123 .amb .ann .bwt.2bit.64 .pac; do
        [[ -s "${prefix}${suffix}" ]] || return 1
    done
}

manifest_value() {
    local manifest=$1
    local key=$2
    awk -F '\t' -v wanted="$key" '$1 == wanted { print $2; exit }' "$manifest"
}

parse_bwa_version() {
    local output=$1
    local line

    while IFS= read -r line; do
        line=${line//$'\r'/}
        line=${line//$'\t'/ }
        if [[ "$line" =~ ^[0-9]+([.][0-9]+)+([-+._[:alnum:]]*)?$ ]]; then
            printf '%s\n' "$line"
            return 0
        fi
    done <<< "$output"

    return 1
}

fasta=$(cd "$(dirname "$fasta")" && pwd -P)/$(basename "$fasta")
cache_dir=$(mkdir -p "$cache_dir" && cd "$cache_dir" && pwd -P)
fasta_sha256=$(sha256_file "$fasta")
fasta_size=$(file_size "$fasta")
fasta_mtime=$(file_mtime "$fasta")
bwa_version_output=$(bwa-mem2 version 2>&1)
if ! bwa_version=$(parse_bwa_version "$bwa_version_output"); then
    echo "ERROR: could not find a semantic version in 'bwa-mem2 version' output" >&2
    exit 1
fi

genome_root="${cache_dir}/${genome}"
fingerprint_dir="${genome_root}/${fasta_sha256}"
final_dir="${fingerprint_dir}/bwamem2"
index_prefix="${final_dir}/${genome}"
manifest="${final_dir}/index.complete.tsv"
lock_dir="${genome_root}/.${fasta_sha256}.build.lock"
owner_file="${lock_dir}/owner.tsv"

mkdir -p "$genome_root"
[[ -w "$genome_root" ]] || {
    echo "ERROR: reference cache is not writable: $genome_root" >&2
    exit 1
}

manifest_matches() {
    [[ -s "$manifest" ]] || return 1
    [[ "$(manifest_value "$manifest" genome)" == "$genome" ]] || return 1
    [[ "$(manifest_value "$manifest" fasta_sha256)" == "$fasta_sha256" ]] || return 1
    [[ "$(manifest_value "$manifest" bwa_mem2_version)" == "$bwa_version" ]] || return 1
    index_complete "$index_prefix"
}

migrate_legacy_manifest() {
    [[ -s "$manifest" ]] || return 1
    [[ "$(manifest_value "$manifest" genome)" == "$genome" ]] || return 1
    [[ "$(manifest_value "$manifest" fasta_sha256)" == "$fasta_sha256" ]] || return 1
    index_complete "$index_prefix" || return 1

    local recorded_version
    recorded_version=$(manifest_value "$manifest" bwa_mem2_version)
    [[ "$recorded_version" == "Looking to launch executable"* ]] || return 1
    [[ "$recorded_version" == *"/bwa-mem2-${bwa_version}_x64-linux/"* ]] || return 1

    local original_created_utc
    local temporary_manifest
    original_created_utc=$(manifest_value "$manifest" created_utc)
    temporary_manifest=$(mktemp "${manifest}.XXXXXX") || {
        echo "ERROR: could not create a temporary index manifest" >&2
        exit 1
    }
    if ! cat > "$temporary_manifest" <<EOF
schema_version	2
genome	${genome}
fasta_path	${fasta}
fasta_size	${fasta_size}
fasta_mtime_epoch	${fasta_mtime}
fasta_sha256	${fasta_sha256}
bwa_mem2_version	${bwa_version}
created_utc	${original_created_utc}
version_metadata_migrated_utc	$(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF
    then
        rm -f "$temporary_manifest"
        echo "ERROR: could not write migrated index metadata" >&2
        exit 1
    fi
    mv "$temporary_manifest" "$manifest" || {
        rm -f "$temporary_manifest"
        echo "ERROR: could not publish migrated index metadata: $manifest" >&2
        exit 1
    }
    echo "Updated legacy bwa-mem2 version metadata: $manifest"
}

write_ready() {
    cat > "$ready_file" <<EOF
genome	${genome}
fasta	${fasta}
fasta_sha256	${fasta_sha256}
index_prefix	${index_prefix}
manifest	${manifest}
bwa_mem2_version	${bwa_version}
EOF
}

if [[ ! -d "$lock_dir" ]]; then
    migrate_legacy_manifest || true
fi
if manifest_matches; then
    echo "Using validated bwa-mem2 index: $index_prefix"
    write_ready
    exit 0
fi

lock_owned=false
tmp_dir=''
cleanup() {
    if [[ -n "$tmp_dir" && -d "$tmp_dir" ]]; then
        rm -rf "$tmp_dir"
    fi
    if [[ "$lock_owned" == true && -d "$lock_dir" ]]; then
        rm -rf "$lock_dir"
    fi
}
trap cleanup EXIT

waited=0
while ! mkdir "$lock_dir" 2>/dev/null; do
    now=$(date +%s)
    lock_epoch=0
    lock_host=unknown
    lock_pid=unknown
    if [[ -s "$owner_file" ]]; then
        lock_epoch=$(manifest_value "$owner_file" created_epoch)
        lock_host=$(manifest_value "$owner_file" hostname)
        lock_pid=$(manifest_value "$owner_file" pid)
    fi
    lock_age=$((now - lock_epoch))
    if (( lock_epoch > 0 && lock_age >= stale_lock_seconds )); then
        if [[ "$lock_host" == "$(hostname)" && "$lock_pid" =~ ^[0-9]+$ ]] \
            && ! kill -0 "$lock_pid" 2>/dev/null; then
            echo "Removing stale local lock from dead PID $lock_pid: $lock_dir" >&2
            rm -rf "$lock_dir"
            continue
        fi
        echo "ERROR: stale or unverifiable index lock detected: $lock_dir" >&2
        echo "Owner host=$lock_host pid=$lock_pid age_seconds=$lock_age" >&2
        echo "Verify that no build is active before removing this lock." >&2
        exit 1
    fi
    if (( waited >= lock_timeout_seconds )); then
        echo "ERROR: timed out waiting for index lock: $lock_dir" >&2
        exit 1
    fi
    sleep 60
    waited=$((waited + 60))
    if manifest_matches; then
        write_ready
        exit 0
    fi
done

lock_owned=true
cat > "$owner_file" <<EOF
hostname	$(hostname)
pid	$$
created_epoch	$(date +%s)
genome	${genome}
fasta_sha256	${fasta_sha256}
EOF

if manifest_matches; then
    write_ready
    exit 0
fi
migrate_legacy_manifest || true
if manifest_matches; then
    write_ready
    exit 0
fi

tmp_dir=$(mktemp -d "${genome_root}/.build.${fasta_sha256}.XXXXXX")
tmp_prefix="${tmp_dir}/${genome}"
echo "Building bwa-mem2 index for $genome in $tmp_dir"
bwa-mem2 index -p "$tmp_prefix" "$fasta"
index_complete "$tmp_prefix" || {
    echo "ERROR: bwa-mem2 did not create a complete index for $genome" >&2
    exit 1
}

cat > "${tmp_dir}/index.complete.tsv" <<EOF
schema_version	2
genome	${genome}
fasta_path	${fasta}
fasta_size	${fasta_size}
fasta_mtime_epoch	${fasta_mtime}
fasta_sha256	${fasta_sha256}
bwa_mem2_version	${bwa_version}
created_utc	$(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF

mkdir -p "$fingerprint_dir"
if [[ -e "$final_dir" ]]; then
    quarantine="${fingerprint_dir}/bwamem2.incomplete.$(date -u +%Y%m%dT%H%M%SZ).$$"
    echo "Moving incomplete index aside: $quarantine" >&2
    mv "$final_dir" "$quarantine"
fi
mv "$tmp_dir" "$final_dir"
tmp_dir=''

manifest_matches || {
    echo "ERROR: published bwa-mem2 index failed final validation: $final_dir" >&2
    exit 1
}

write_ready
echo "Published validated bwa-mem2 index: $index_prefix"
