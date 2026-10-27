#!/usr/bin/env bash
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # don't hide errors within pipes

# Script variables
input_bam="${1}"
output_bg="${2}"
binsize="${3}"
blacklist="${4}"
threads="${5}"
read_type="${6:-PE}"     # Default to PE if not specified
read_length="${7:-150}"  # Default to 150bp if not specified

# Set up extension parameters based on read type
if [ "$read_type" == "PE" ]; then
    extend_param="--centerReads"
else
    extend_param="--extendReads ${read_length}"
fi

# Run bamCoverage with appropriate parameters
bamCoverage \
    --bam "$input_bam" \
    --outFileName "$output_bg" \
    --binSize "$binsize" \
    --numberOfProcessors "$threads" \
    --blackListFileName "$blacklist" \
    $extend_param
