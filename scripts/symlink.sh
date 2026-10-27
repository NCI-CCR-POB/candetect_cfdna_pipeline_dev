#!/usr/bin/env bash
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # don't hide errors within pipes

# Script variables
input_read1="${1}"
output_read1="${2}"
output_read2="${3}"
outdir="${4}"

# Create output directory
mkdir -p "$outdir"

# Get the R2 file path
input_read2="$(echo $input_read1 | sed "s/_R1/_R2/g")"

# Check if input files exist
if [ ! -f "$input_read1" ]; then
    echo "Error: Input file 1 not found: $input_read1"
    exit 1
fi

if [ ! -f "$input_read2" ]; then
    echo "Error: Input file 2 not found: $input_read2"
    exit 1
fi

# Create symlinks using absolute paths
ln -sf "$input_read1" "$output_read1"
ln -sf "$input_read2" "$output_read2"
