#!/usr/bin/env bash
set -o errexit
set -o nounset
set -o pipefail

# Script variables
input_r1="${1}"
output_r1="${2}"
output_r2="${3}"
outdir="${4}"

# Create output directory
mkdir -p "$outdir"

# Get potential R2 and R3 paths
input_r2="$(echo $input_r1 | sed "s/_R1/_R2/g")"
input_r3="$(echo $input_r1 | sed "s/_R1/_R3/g")"

# Log file setup
log_file="${outdir}/read_pair_detection.log"
echo "=== Read Pair Detection Log ===" > "$log_file"
echo "Timestamp: $(date)" >> "$log_file"
echo "Input R1: $input_r1" >> "$log_file"
echo "Potential R2: $input_r2" >> "$log_file"
echo "Potential R3: $input_r3" >> "$log_file"

# Check file existence and sizes
if [ -f "$input_r2" ] && [ -f "$input_r3" ]; then
    echo "Both R2 and R3 files exist" >> "$log_file"
    size_r1=$(stat -c%s "$input_r1")
    size_r2=$(stat -c%s "$input_r2")
    size_r3=$(stat -c%s "$input_r3")
    
    echo "File sizes:" >> "$log_file"
    echo "R1: $size_r1 bytes" >> "$log_file"
    echo "R2: $size_r2 bytes" >> "$log_file"
    echo "R3: $size_r3 bytes" >> "$log_file"
    
    # Calculate size differences
    diff_r2=$((size_r1 - size_r2))
    diff_r3=$((size_r1 - size_r3))
    echo "Size differences from R1:" >> "$log_file"
    echo "R2 difference: $diff_r2 bytes" >> "$log_file"
    echo "R3 difference: $diff_r3 bytes" >> "$log_file"
    
    # If R3 is closer in size to R1, use R3 as the pair
    if [ ${diff_r3#-} -lt ${diff_r2#-} ]; then
        echo "Decision: Using R3 as pair for R1 (smaller size difference)" >> "$log_file"
        ln -sf "$input_r1" "$output_r1"
        ln -sf "$input_r3" "$output_r2"
        echo "Created symlinks:" >> "$log_file"
        ls -lh "$output_r1" "$output_r2" >> "$log_file"
    else
        echo "Decision: Using R2 as pair for R1 (smaller size difference)" >> "$log_file"
        ln -sf "$input_r1" "$output_r1"
        ln -sf "$input_r2" "$output_r2"
        echo "Created symlinks:" >> "$log_file"
        ls -lh "$output_r1" "$output_r2" >> "$log_file"
    fi
elif [ -f "$input_r2" ]; then
    echo "Only R2 exists" >> "$log_file"
    echo "Decision: Using R2 as pair for R1" >> "$log_file"
    ln -sf "$input_r1" "$output_r1"
    ln -sf "$input_r2" "$output_r2"
    echo "Created symlinks:" >> "$log_file"
    ls -lh "$output_r1" "$output_r2" >> "$log_file"
elif [ -f "$input_r3" ]; then
    echo "Only R3 exists" >> "$log_file"
    echo "Decision: Using R3 as pair for R1" >> "$log_file"
    ln -sf "$input_r1" "$output_r1"
    ln -sf "$input_r3" "$output_r2"
    echo "Created symlinks:" >> "$log_file"
    ls -lh "$output_r1" "$output_r2" >> "$log_file"
else
    echo "Error: No matching read file found for $input_r1" >> "$log_file"
    exit 1
fi

echo "=== End Read Pair Detection Log ===" >> "$log_file"