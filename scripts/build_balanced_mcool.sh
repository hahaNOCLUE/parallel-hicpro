#!/bin/bash

set -eo pipefail

dir=$(cd "$(dirname "$0")" && pwd)

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

die() { echo "Exit: $*" >&2; exit 1; }
while IFS="=" read -r key value; do
    key=$(echo "$key" | sed -e "s/^[[:space:]]*//" -e "s/[[:space:]]*$//")
    value=$(echo "$value" | sed -e "s/#.*//" -e "s/^[[:space:]]*//" -e "s/[[:space:]]*$//")
    if [[ "$key" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]]; then printf -v "$key" "%s" "$value"; export "$key"; fi
done < "$conf_file"

cooler_bin=${COOLER_BIN:-cooler}
bin_size=${MCOOL_BIN_SIZE:-${BIN_SIZE%% *}}
resolutions=${MCOOL_RESOLUTIONS:-${BIN_SIZE// /,}}
assembly=${MCOOL_ASSEMBLY:-GRCh38}
out_dir=${MCOOL_OUTPUT_DIR:-${MAPC_OUTPUT}/mcool}
chunksize=${MCOOL_CHUNKSIZE:-5000000}

if ! "$cooler_bin" --version >/dev/null 2>&1; then
    die "cooler is unavailable: $cooler_bin"
fi

chromsizes=${MCOOL_CHROMSIZES:-${GENOME_SIZE:-}}
tmp_chromsizes=
if [ -z "$chromsizes" ]; then
    reference_dir=${BOWTIE2_IDX_PATH}
    dict=${reference_dir}/${REFERENCE_GENOME}.dict
    fai=$(find "$reference_dir" -maxdepth 1 -name '*.fai' -print -quit)
    tmp_chromsizes=${TMP_DIR}/.${REFERENCE_GENOME}.$$.chrom.sizes
    if [ -n "$fai" ]; then
        cut -f1,2 "$fai" > "$tmp_chromsizes"
    elif [ -f "$dict" ]; then
        awk -F'\t' '$1=="@SQ" {sub(/^SN:/,"",$2); sub(/^LN:/,"",$3); print $2"\t"$3}' "$dict" > "$tmp_chromsizes"
    else
        die "set MCOOL_CHROMSIZES or provide a FASTA .fai / ${REFERENCE_GENOME}.dict"
    fi
    chromsizes=$tmp_chromsizes
fi
trap 'if [ -n "$tmp_chromsizes" ]; then rm -f "$tmp_chromsizes"; fi' EXIT

mkdir -p "$out_dir" "$TMP_DIR"
found=0
for pairs in "${MAPC_OUTPUT}"/data/*/*.allValidPairs; do
    [ -f "$pairs" ] || continue
    found=1
    sample=$(basename "$pairs" .allValidPairs)
    base_cool=${out_dir}/${sample}.${bin_size}.cool
    mcool=${out_dir}/${sample}.mcool
    logfile=${out_dir}/${sample}.mcool.log

    echo "Building balanced mcool: $mcool"
    rm -f "$base_cool" "$mcool"
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{keep[$1]=1; next} ($2 in keep) && ($5 in keep)' "$chromsizes" "$pairs" | \
        "$cooler_bin" cload pairs \
            -c1 2 -p1 3 -c2 5 -p2 6 \
            --input-copy-status unique --assembly "$assembly" \
            --chunksize "$chunksize" --temp-dir "$TMP_DIR" \
            "${chromsizes}:${bin_size}" - "$base_cool" >"$logfile" 2>&1

    "$cooler_bin" zoomify --nproc "${N_CPU:-1}" \
        --resolutions "$resolutions" --balance \
        --balance-args "--nproc ${N_CPU:-1} --ignore-diags ${MCOOL_IGNORE_DIAGS:-2} --convergence-policy store_final" \
        --out "$mcool" "$base_cool" >>"$logfile" 2>&1
    rm -f "$base_cool"
    "$cooler_bin" ls "$mcool" >>"$logfile" 2>&1
done

if [ "$found" -eq 0 ]; then
    die "no sample-level .allValidPairs files found under ${MAPC_OUTPUT}/data"
fi
