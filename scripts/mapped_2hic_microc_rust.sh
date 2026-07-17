#!/bin/bash

set -eo pipefail

dir=$(cd "$(dirname "$0")" && pwd)
repo_dir=$(cd "$dir/.." && pwd)

while [ "$#" -gt 0 ]; do
    case "$1" in
        -c) conf_file=$2; shift 2 ;;
        -h|--help) echo "usage: $0 -c CONFIG"; exit 0 ;;
        *) echo "$0: unknown option $1" >&2; exit 2 ;;
    esac
done

if [ -z "${conf_file:-}" ]; then
    echo "$0: -c CONFIG is required" >&2
    exit 2
fi

CONF=$conf_file . "$dir/hic.inc.sh"

if [ -n "${GENOME_FRAGMENT:-}" ]; then
    die "mapped_2hic_microc_rust.sh is only for Micro-C/DNase mode (GENOME_FRAGMENT must be empty)"
fi

rust_bin=${HICPRO_MICROC_RUST_BIN:-$dir/hicpro-microc-rs}
breakpoints=${HICPRO_BREAKPOINTS:-$repo_dir/annotations/GRCh38_noalt_decoy_as.breakpoints.high_confidence.tsv}

if [ ! -x "$rust_bin" ]; then
    die "Rust launcher not found: $rust_bin"
fi
if [ ! -f "$breakpoints" ]; then
    die "Breakpoint file not found: $breakpoints"
fi

opts=(--threads "${N_CPU:-1}" --breakpoints "$breakpoints" -v)
if [ "${GET_ALL_INTERACTION_CLASSES:-0}" -eq 1 ]; then opts+=(-a); fi
if [ -n "${MIN_CIS_DIST:-}" ] && [ "$MIN_CIS_DIST" -ge 0 ]; then opts+=(-d "$MIN_CIS_DIST"); fi
if [ -n "${ALLELE_SPECIFIC_SNP:-}" ]; then opts+=(-g XA); fi

for r in $(get_paired_bam); do
    sample_dir=$(get_sample_dir "$r")
    datadir=${MAPC_OUTPUT}/data/${sample_dir}
    ldir=${LOGS_DIR}/${sample_dir}
    mkdir -p "$datadir" "$ldir"
    logfile=${ldir}/mapped_2hic_dnase.log

    shard_opts=()
    if [ "${RUST_SHARD_VALIDPAIRS:-0}" -eq 1 ]; then
        shard_opts+=(--shard-dir "$datadir/shards")
    fi

    echo "Logs: $logfile"
    "$rust_bin" "${opts[@]}" "${shard_opts[@]}" -r "$r" -o "$datadir" >"$logfile" 2>&1

    out_valid=$(basename "$r" | sed -e 's/.bam$/.validPairs/')
    LC_ALL=C sort --parallel="${N_CPU:-1}" -T "$TMP_DIR" \
        -k2,2V -k3,3n -k5,5V -k6,6n \
        -o "$datadir/$out_valid" "$datadir/$out_valid" >>"$logfile" 2>&1
done
