#!/bin/bash
set -euo pipefail
command -v bowtie >/dev/null 2>&1 || { echo >&2 "ERROR: bowtie not found. Please install bowtie first."; exit 1; }
command -v bowtie-build >/dev/null 2>&1 || { echo >&2 "ERROR: bowtie-build not found. Please install bowtie first."; exit 1; }

# Automatic path based on script location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${ROOT_DIR}/output/bowtie"
mkdir -p "${OUTPUT_DIR}"

# Usage
if [ $# -ne 3 ]; then
    echo "Usage: $0 <reference_fasta> <input_spacer_fasta> <output_prefix>"
    exit 1
fi

REF_FA="$1"
INPUT_FA="$2"
PREFIX="$3"
INDEX="${OUTPUT_DIR}/bowtie_index"

echo "=== Build bowtie index ==="
bowtie-build --quiet "${REF_FA}" "${INDEX}"

echo "=== Run bowtie alignment ==="
for V in 0 1 2; do
    echo "Processing ${V} mismatches..."
    bowtie -v ${V} -a --norc --quiet "${INDEX}" "${INPUT_FA}" \
        | awk '{print $1 "\t" $NF}' \
        | sort \
        | uniq -c \
        | awk '{print $1 "\t" $2}' \
        > "${OUTPUT_DIR}/${PREFIX}_mm${V}_counts.txt"
done

echo "=== Done ==="
echo "Outputs: ${OUTPUT_DIR}"

