#!/usr/bin/env bash
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # don't hide errors within pipes

# Script variables
input_bam=$1
blacklist=$2
threads=$3
output_tsv=$4
output_pdf=$5

# Run plotCoverage for single BAM file
plotCoverage \
    --bamfiles "$input_bam" \
    --numberOfProcessors "$threads" \
    --blackListFileName "$blacklist" \
    --plotFile "$output_pdf" \
    --outRawCounts "$output_tsv" \
    --skipZeros \
    --centerReads
