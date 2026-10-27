#!/usr/bin/env python3
"""
Script to visualize f-profiles in stacked barplots e.g. https://www.pnas.org/doi/10.1073/pnas.2220982120#sec-3 Fig7a
"""
import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt

def load_inputs():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Visualize f-profile contributions across samples as a stacked bar chart.")
    parser.add_argument("--input_file", type=str, required=True, help="Path to the f-profile TSV file")
    parser.add_argument("--metadata_file", type=str, required=False, help="Optional metadata TSV file with sample name, cohort, and cycle. Use 'none' to ignore.")
    parser.add_argument("--output_file", type=str, required=True, help="Path to save the output plot")
    parser.add_argument("--title", type=str, required=False, default="Stacked Bar Graph of f-Profile Contributions by Sample", help="Title of the plot")
    parser.add_argument("--colormap", type=str, required=False, default="viridis", help="Colormap to use for the plot (default: 'viridis')")

    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        raise FileNotFoundError(f"Input TSV file not found: {args.input_file}")
    
    if args.metadata_file and args.metadata_file.lower() != "none" and not os.path.exists(args.metadata_file):
        raise FileNotFoundError(f"Metadata TSV file not found: {args.metadata_file}")

    return args

def preprocess_metadata(df, metadata_file):
    """Preprocess metadata to match and order samples."""
    metadata = pd.read_csv(metadata_file, sep="\t")
    
    # Support either 'library' or 'sample' as the sample identifier column
    id_col = "library" if "library" in metadata.columns else "sample"
    metadata["library_motifs"] = metadata[id_col] + "_motifs"

    # Ensure metadata sample names match with input file columns
    sort_cols = [c for c in ["cohort", "cycle", id_col] if c in metadata.columns]
    sample_order = metadata.sort_values(by=sort_cols)["library_motifs"]
    
    # Filter input data to match ordered sample names
    df = df[sample_order]
    return df

def visualize_fprofiles(input_file, metadata_file=None, output_file=None, title="Stacked Bar Graph", colormap="viridis"):
    """Generate a stacked bar chart from f-profile data, optionally using metadata for sorting and grouping."""
    print(f"Loading input file: {input_file}")
    df = pd.read_csv(input_file, sep="\t", index_col=0)

    if metadata_file and metadata_file.lower() != "none":
        print(f"Processing metadata file: {metadata_file}")
        df = preprocess_metadata(df, metadata_file)

    print("Generating plot...")
    ax = df.T.plot(kind='bar', stacked=True, figsize=(12, 7), colormap=colormap, edgecolor='black', width=0.8)

    plt.ylabel('Percentage Contribution (%)')
    plt.xlabel('Samples')
    plt.title(title)
    plt.ylim(0, 100)
    plt.legend(title='f-Profile', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig(output_file, format=output_file.split('.')[-1], dpi=600)
    print(f"Plot saved to {output_file}")

if __name__ == "__main__":
    args = load_inputs()
    visualize_fprofiles(args.input_file, args.metadata_file, args.output_file, args.title, args.colormap)
