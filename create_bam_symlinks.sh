#!/bin/bash

# Set directories
FRAG_BAMS="/data/ShernLiquidBx/GBC/EMseq/analysis/cfdna_wgs/bams"  # Directory where symlinks will be created
LIBRARIES_FILE="/data/ShernLiquidBx/GBC/EMseq/inputs/libraries_withhealthy.tsv"

# Create output directory if it doesn't exist
mkdir -p "$FRAG_BAMS"

# Read the libraries file and create symlinks
while IFS=$'\t' read -r library cohort file library_name library_type; do
    # Skip header line
    if [[ "$library" == "library" ]]; then
        continue
    fi
    
    # Create symlink
    target="$FRAG_BAMS/${library}_filt.bam"
    echo "Creating symlink: $target -> $file"
    ln -sf "$file" "$target"
    
    # Verify symlink was created
    if [ -L "$target" ]; then
        echo "Successfully created symlink for library $library"
    else
        echo "Failed to create symlink for library $library"
    fi
done < "$LIBRARIES_FILE"

echo "Symlink creation completed" 