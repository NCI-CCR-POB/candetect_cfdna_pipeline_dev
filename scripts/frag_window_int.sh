#!/usr/bin/env bash
input=$1
keep_bed=$2
output=$3
genome=$4

#module load bedtools

# Function to check file validity
check_bed() {
    local bed_file=$1
    local desc=$2
    
    if [[ ! -f $bed_file ]]; then
        echo "ERROR: $desc file does not exist: $bed_file"
        exit 1
    fi
    
    if [[ ! -s $bed_file ]]; then
        echo "ERROR: $desc file is empty: $bed_file"
        exit 1
    fi
    
    # Check if file is readable and valid BED format
    head -n 1 $bed_file || {
        echo "ERROR: Cannot read $desc file: $bed_file"
        exit 1
    }
}

# Check input files
check_bed "$input" "Input BED"
check_bed "$keep_bed" "Keep BED"
check_bed "$genome" "Genome file"

# Create temporary directory
tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

# Sort both input files
echo "Sorting input files..."

# Sort input BED file
echo "Sorting input BED file..."
sort -k1,1 -k2,2n "$input" > "${tmpdir}/pre_sorted_input.bed"
bedtools sort -i "${tmpdir}/pre_sorted_input.bed" -g "$genome" > "${tmpdir}/sorted_input.bed" || {
    echo "ERROR: Failed to sort input file"
    exit 1
}

# Sort reference BED file
echo "Sorting reference BED file..."
sort -k1,1 -k2,2n "$keep_bed" > "${tmpdir}/pre_sorted_ref.bed"
bedtools sort -i "${tmpdir}/pre_sorted_ref.bed" -g "$genome" > "${tmpdir}/sorted_ref.bed" || {
    echo "ERROR: Failed to sort reference file"
    exit 1
}

# Run intersection with sorted files
echo "Running intersection..."
bedtools intersect -sorted -g "$genome" -c \
    -a "${tmpdir}/sorted_ref.bed" \
    -b "${tmpdir}/sorted_input.bed" > "$output" || {
    echo "ERROR: Failed to run intersection"
    exit 1
}

# Verify output
if [[ ! -s "$output" ]]; then
    echo "ERROR: Output file is empty: $output"
    exit 1
fi

echo "Successfully created: $output"
echo "First few lines of output:"
head -n 3 "$output"

