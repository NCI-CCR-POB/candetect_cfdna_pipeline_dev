#!/bin/bash

# Set the base directory where all FRAGLE results are stored
FRAGLE_DIR="/Volumes/ShernLiquidBx/SARC031_4-3-25/FRAGLE"

# Get current date for the master file
CURRENT_DATE=$(date +%Y-%m-%d)

# First, process each subdirectory to create individual compiled files
# Using find to get all directories and exclude the root directory itself
find "$FRAGLE_DIR" -mindepth 1 -maxdepth 1 -type d ! -name ".*" | while read dir; do
    # Get the sample ID from the directory name
    SAMPLE_ID=$(basename "$dir")
    
    # Check if the required CSV files exist
    if [ ! -f "$dir/Fragle.csv" ] || [ ! -f "$dir/HT.csv" ] || [ ! -f "$dir/LT.csv" ]; then
        echo "Warning: Missing required CSV files in $dir"
        continue
    fi
    
    # Extract ctDNA burden values using awk (skip header line)
    FRAGLE_BURDEN=$(awk -F',' 'NR==2 {print $2}' "$dir/Fragle.csv")
    HT_BURDEN=$(awk -F',' 'NR==2 {print $2}' "$dir/HT.csv")
    LT_BURDEN=$(awk -F',' 'NR==2 {print $2}' "$dir/LT.csv")
    
    # Create the compiled CSV for this sample
    echo "sample_id,fragle_ctdna_burden,ht_ctdna_burden,lt_ctdna_burden" > "$dir/${SAMPLE_ID}_compiled_FRAGLE.csv"
    echo "$SAMPLE_ID,$FRAGLE_BURDEN,$HT_BURDEN,$LT_BURDEN" >> "$dir/${SAMPLE_ID}_compiled_FRAGLE.csv"
done

# Create the master compiled CSV
MASTER_CSV="$FRAGLE_DIR/master_compiled_csv_${CURRENT_DATE}.csv"

# Write header to master CSV
echo "sample_id,fragle_ctdna_burden,ht_ctdna_burden,lt_ctdna_burden" > "$MASTER_CSV"

# Combine all individual compiled CSVs (skip headers)
find "$FRAGLE_DIR" -mindepth 2 -maxdepth 2 -name "*_compiled_FRAGLE.csv" | while read compiled; do
    if [ -f "$compiled" ]; then
        tail -n +2 "$compiled" >> "$MASTER_CSV"
    fi
done

echo "Processing complete. Master file created at: $MASTER_CSV" 