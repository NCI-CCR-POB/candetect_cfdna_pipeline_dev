#!/usr/bin/env bash
set -euo pipefail

# Get input arguments
input_files="$1"    
prefix="$2"         
outdir="$3"         

# Create output directory if it doesn't exist
mkdir -p "${outdir}"

# Debug information
echo "=== Debug Information ==="
echo "Working directory: $(pwd)"
echo "Input files: ${input_files}"
echo "Prefix: ${prefix}"
echo "Output dir: ${outdir}"
echo "======================="

# Run MultiQC with the input files
multiqc \
    --filename "${prefix}" \
    --force \
    --outdir "${outdir}" \
    --verbose \
    ${input_files}

# List created files
echo "=== Created Files ==="
ls -la "${outdir}"
echo "=== Data Directory Contents ==="
ls -la "${outdir}/${prefix}_data" || echo "Data directory not found"
echo "====================="

# Check if main HTML file exists
if [ ! -f "${outdir}/${prefix}.html" ]; then
    echo "ERROR: MultiQC HTML report not created"
    exit 1
fi

echo "MultiQC completed successfully"
