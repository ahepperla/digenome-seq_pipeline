#!/usr/bin/env bash
# Find or build the bwa-mem2 index for a FASTA in the shared reference cache,
# then print the index prefix on stdout. Progress messages go to stderr.
#
# Layout: <cache>/<genome>/<fasta_sha256>/<bwa_mem2_version>/bwamem2/<genome>.*
# An index is reused only when its manifest matches the genome name, FASTA
# SHA-256, and bwa-mem2 version, and every index file is present.
#
# Several users and runs share the cache, so:
# - one builder at a time holds a mkdir lock with an owner record; a stale
#   lock is removed automatically only when its owner is a dead PID on this host;
# - the index is built in a temporary directory, validated, then moved into
#   place; an incomplete index is moved aside, never deleted;
# - modes are set explicitly, not left to umask: namespace directories through
#   the version level are 1777, finished index directories 755, index files 644,
#   and lock directories and owner records are world-readable. An extra setgid
#   bit is accepted.
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: prepare_bwamem2_index.sh --genome NAME --fasta FASTA --cache-dir DIR
       [--lock-timeout-seconds 172800] [--stale-lock-seconds 172800]
EOF
    exit 2
}

genome=''
fasta=''
cache_dir=''
lock_timeout_seconds=172800
stale_lock_seconds=172800

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genome) genome=$2; shift 2 ;;
        --fasta) fasta=$2; shift 2 ;;
        --cache-dir) cache_dir=$2; shift 2 ;;
        --lock-timeout-seconds) lock_timeout_seconds=$2; shift 2 ;;
        --stale-lock-seconds) stale_lock_seconds=$2; shift 2 ;;
        *) usage ;;
    esac
done

[[ -n "$genome" && -n "$fasta" && -n "$cache_dir" ]] || usage
[[ -s "$fasta" ]] || { echo "ERROR: FASTA does not exist or is empty: $fasta" >&2; exit 1; }
command -v bwa-mem2 >/dev/null || { echo "ERROR: bwa-mem2 is not available" >&2; exit 1; }

log() {
    echo "$@" >&2
}

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

sha256_file() {
    if command -v sha256sum >/dev/null; then
        sha256sum "$1" | awk '{print $1}'
    else
        shasum -a 256 "$1" | awk '{print $1}'
    fi
}

# GNU stat (Linux containers) or BSD stat (the macOS test machine).
path_mode() {
    stat -c %a "$1" 2>/dev/null || stat -f %Lp "$1"
}

mode_satisfies_requirement() {
    local current=$((8#$1)) required=$((8#$2)) setgid=$((8#2000))
    (( (current & ~setgid) == required ))
}

ensure_directory_mode() {
    local path=$1 mode=$2 label=$3 current
    mkdir -p "$path" || fail "could not create $label: $path"
    current=$(path_mode "$path")
    if ! mode_satisfies_requirement "$current" "$mode" && ! chmod "$mode" "$path"; then
        fail "$label must have mode $mode but is $current: $path"
    fi
}

make_index_world_readable() {
    find "$1" -type d -exec chmod 755 {} + || fail "could not make index directories traversable: $1"
    find "$1" -type f -exec chmod 644 {} + || fail "could not make index files readable: $1"
}

index_complete() {
    local prefix=$1 suffix
    for suffix in .0123 .amb .ann .bwt.2bit.64 .pac; do
        [[ -s "${prefix}${suffix}" ]] || return 1
    done
}

manifest_value() {
    awk -F '\t' -v wanted="$2" '$1 == wanted { print $2; exit }' "$1"
}

# `bwa-mem2 version` may print CPU-dispatch messages before the version line.
bwa_version() {
    local line
    while IFS= read -r line; do
        line=${line//$'\r'/}
        if [[ "$line" =~ ^[0-9]+([.][0-9]+)+([-+._[:alnum:]]*)?$ ]]; then
            printf '%s\n' "$line"
            return 0
        fi
    done < <(bwa-mem2 version 2>&1)
    return 1
}

fasta=$(cd "$(dirname "$fasta")" && pwd -P)/$(basename "$fasta")
ensure_directory_mode "$cache_dir" 1777 "reference cache root"
cache_dir=$(cd "$cache_dir" && pwd -P)
fasta_sha256=$(sha256_file "$fasta")
version=$(bwa_version) || fail "could not find a semantic version in 'bwa-mem2 version' output"

genome_dir="${cache_dir}/${genome}"
fasta_dir="${genome_dir}/${fasta_sha256}"
version_dir="${fasta_dir}/${version}"
final_dir="${version_dir}/bwamem2"
index_prefix="${final_dir}/${genome}"
manifest="${final_dir}/index.complete.tsv"
lock_dir="${genome_dir}/.${fasta_sha256}.${version}.build.lock"
owner_file="${lock_dir}/owner.tsv"

ensure_directory_mode "$genome_dir" 1777 "genome cache namespace"
[[ -w "$genome_dir" ]] || fail "reference cache is not writable: $genome_dir"

index_ready() {
    [[ -s "$manifest" ]] || return 1
    [[ "$(manifest_value "$manifest" genome)" == "$genome" ]] || return 1
    [[ "$(manifest_value "$manifest" fasta_sha256)" == "$fasta_sha256" ]] || return 1
    [[ "$(manifest_value "$manifest" bwa_mem2_version)" == "$version" ]] || return 1
    index_complete "$index_prefix"
}

finish() {
    log "Using bwa-mem2 index: $index_prefix"
    printf '%s\n' "$index_prefix"
    exit 0
}

index_ready && finish

lock_owned=false
tmp_dir=''
cleanup() {
    [[ -z "$tmp_dir" || ! -d "$tmp_dir" ]] || rm -rf "$tmp_dir"
    [[ "$lock_owned" != true || ! -d "$lock_dir" ]] || rm -rf "$lock_dir"
}
trap cleanup EXIT

# Wait for another builder, or remove its lock if it is provably dead.
waited=0
while ! mkdir "$lock_dir" 2>/dev/null; do
    lock_epoch=0 lock_host=unknown lock_pid=unknown
    if [[ -s "$owner_file" ]]; then
        lock_epoch=$(manifest_value "$owner_file" created_epoch)
        lock_host=$(manifest_value "$owner_file" hostname)
        lock_pid=$(manifest_value "$owner_file" pid)
    fi
    lock_age=$(( $(date +%s) - lock_epoch ))
    if (( lock_epoch > 0 && lock_age >= stale_lock_seconds )); then
        if [[ "$lock_host" == "$(hostname)" && "$lock_pid" =~ ^[0-9]+$ ]] \
            && ! kill -0 "$lock_pid" 2>/dev/null; then
            log "Removing stale local lock from dead PID $lock_pid: $lock_dir"
            rm -rf "$lock_dir"
            continue
        fi
        log "Owner host=$lock_host pid=$lock_pid age_seconds=$lock_age"
        log "Verify that no build is active before removing this lock."
        fail "stale or unverifiable index lock detected: $lock_dir"
    fi
    (( waited < lock_timeout_seconds )) || fail "timed out waiting for index lock: $lock_dir"
    sleep 60
    waited=$((waited + 60))
    index_ready && finish
done

lock_owned=true
chmod 755 "$lock_dir" || fail "could not make the index lock world-readable: $lock_dir"
cat > "$owner_file" <<EOF
hostname	$(hostname)
pid	$$
created_epoch	$(date +%s)
genome	${genome}
fasta_sha256	${fasta_sha256}
EOF
chmod 644 "$owner_file" || fail "could not make the lock owner record readable: $owner_file"

index_ready && finish
ensure_directory_mode "$fasta_dir" 1777 "FASTA cache namespace"
ensure_directory_mode "$version_dir" 1777 "bwa-mem2 version namespace"

tmp_dir=$(mktemp -d "${genome_dir}/.build.${fasta_sha256}.${version}.XXXXXX")
log "Building bwa-mem2 index for $genome in $tmp_dir"
bwa-mem2 index -p "${tmp_dir}/${genome}" "$fasta" >&2
index_complete "${tmp_dir}/${genome}" || fail "bwa-mem2 did not create a complete index for $genome"
cat > "${tmp_dir}/index.complete.tsv" <<EOF
schema_version	3
genome	${genome}
fasta_path	${fasta}
fasta_sha256	${fasta_sha256}
bwa_mem2_version	${version}
created_utc	$(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF
make_index_world_readable "$tmp_dir"

if [[ -e "$final_dir" ]]; then
    quarantine="${version_dir}/bwamem2.incomplete.$(date -u +%Y%m%dT%H%M%SZ).$$"
    log "Moving incomplete index aside: $quarantine"
    mv "$final_dir" "$quarantine"
fi
mv "$tmp_dir" "$final_dir"
tmp_dir=''
index_ready || fail "published bwa-mem2 index failed final validation: $final_dir"
finish
