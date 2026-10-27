#!/usr/bin/env bash
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # don't hide errors within pipes

# Script variables
raw_bam="${1}"
dedup_bam="${2}"
threads="${3}"

# Create temporary directory
tmp_dir=$(mktemp -d -p $(dirname "$dedup_bam"))
trap 'rm -rf "$tmp_dir"' EXIT

# Step 1: Sort by name
echo "Sorting by name..."
samtools sort -@ "$threads" -n -o "$tmp_dir/namesort.bam" "$raw_bam"

# Step 2: Fix mates
echo "Running fixmate..."
samtools fixmate -m "$tmp_dir/namesort.bam" "$tmp_dir/fixmate.bam"

# Step 3: Sort by position
echo "Sorting by position..."
samtools sort -@ "$threads" -o "$tmp_dir/positionsort.bam" "$tmp_dir/fixmate.bam"

# Step 4: Mark and remove duplicates
echo "Marking and removing duplicates..."
samtools markdup -@ "$threads" -r "$tmp_dir/positionsort.bam" "$dedup_bam"

# Step 5: Index the final BAM
echo "Indexing final BAM..."
samtools index "$dedup_bam"
