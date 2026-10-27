#!/usr/bin/env python3
import os
import glob
import pandas as pd
import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser(description='Extract tumor fraction features from ichorCNA output')
    parser.add_argument('--ichor_dir', required=True, help='Directory containing ichorCNA output')
    parser.add_argument('--output', help='Output file path (takes precedence if both specified)')
    parser.add_argument('--out_dir', help='Output directory for features (deprecated)')
    return parser.parse_args()

def extract_tumor_fraction(ichor_dir):
    """Extract tumor fraction from ichorCNA output files"""
    print(f"\nStarting tumor fraction extraction from directory: {ichor_dir}", file=sys.stderr)
    results = []
    
    # Find all .params.txt files in the ichorCNA output directory and its subdirectories
    param_files = glob.glob(os.path.join(ichor_dir, '*.params.txt')) + \
                 glob.glob(os.path.join(ichor_dir, '*', '*.params.txt'))
    
    print(f"\nFound {len(param_files)} parameter files:", file=sys.stderr)
    for pf in param_files:
        print(f"  - {pf}", file=sys.stderr)
    
    if not param_files:
        print(f"ERROR: No .params.txt files found in {ichor_dir}", file=sys.stderr)
        sys.exit(1)
    
    for param_file in param_files:
        # Extract sample ID from the path
        if os.path.dirname(param_file) == ichor_dir:
            # If file is directly in ichor_dir, use filename without extension
            sample_id = os.path.basename(param_file).replace('.params.txt', '')
            print(f"\nProcessing file directly in ichor_dir: {param_file}", file=sys.stderr)
            print(f"  Sample ID: {sample_id}", file=sys.stderr)
        else:
            # If file is in a subdirectory, use the subdirectory name
            sample_id = os.path.basename(os.path.dirname(param_file))
            print(f"\nProcessing file in subdirectory: {param_file}", file=sys.stderr)
            print(f"  Sample ID: {sample_id}", file=sys.stderr)
        
        try:
            # Read the params file
            print(f"  Reading file contents...", file=sys.stderr)
            with open(param_file, 'r') as f:
                lines = f.readlines()
            print(f"  File contains {len(lines)} lines", file=sys.stderr)
            
            # Find the line containing tumor fraction - look for both formats
            tf_line = None
            print("  Searching for tumor fraction line...", file=sys.stderr)
            for i, line in enumerate(lines):
                if 'Tumor Fraction:' in line:
                    print(f"    Found tumor fraction line (format 1) at line {i+1}: {line.strip()}", file=sys.stderr)
                    tf_line = line
                    break
                elif 'Tumor Fraction\t' in line and not line.startswith('Sample'):
                    print(f"    Found tumor fraction line (format 2) at line {i+1}: {line.strip()}", file=sys.stderr)
                    tf_line = line
                    break
            
            if tf_line:
                # Extract tumor fraction value
                if 'Tumor Fraction:' in tf_line:
                    tumor_fraction = float(tf_line.split(':')[1].strip())
                    print(f"    Extracted tumor fraction (format 1): {tumor_fraction}", file=sys.stderr)
                else:
                    tumor_fraction = float(tf_line.split('\t')[1].strip())
                    print(f"    Extracted tumor fraction (format 2): {tumor_fraction}", file=sys.stderr)
                
                results.append({
                    'Sample_ID': sample_id,
                    'Tumor_Fraction': tumor_fraction
                })
                print(f"    Successfully added sample {sample_id} with tumor fraction {tumor_fraction}", file=sys.stderr)
            else:
                print(f"  WARNING: Could not find tumor fraction in {sample_id}", file=sys.stderr)
                print("  File contents:", file=sys.stderr)
                for i, line in enumerate(lines):
                    print(f"    Line {i+1}: {line.strip()}", file=sys.stderr)
        except Exception as e:
            print(f"  ERROR processing {sample_id}: {str(e)}", file=sys.stderr)
            print("  File contents:", file=sys.stderr)
            try:
                for i, line in enumerate(lines):
                    print(f"    Line {i+1}: {line.strip()}", file=sys.stderr)
            except:
                print("    Unable to print file contents", file=sys.stderr)
    
    print(f"\nProcessed {len(param_files)} files, found {len(results)} valid tumor fractions", file=sys.stderr)
    
    if not results:
        print("ERROR: No valid tumor fractions found", file=sys.stderr)
        sys.exit(1)
        
    return pd.DataFrame(results)

def main():
    args = parse_args()
    
    print(f"\nStarting tumor fraction extraction script", file=sys.stderr)
    print(f"Input directory: {args.ichor_dir}", file=sys.stderr)
    
    # Determine output path
    if args.output:
        output_file = args.output
        out_dir = os.path.dirname(output_file)
    elif args.out_dir:
        out_dir = args.out_dir
        output_file = os.path.join(out_dir, 'tf.csv')
    else:
        print("ERROR: Either --output or --out_dir must be specified", file=sys.stderr)
        sys.exit(1)
    
    print(f"Output will be written to: {output_file}", file=sys.stderr)
    
    # Create output directory if it doesn't exist
    os.makedirs(out_dir, exist_ok=True)
    print(f"Created/verified output directory: {out_dir}", file=sys.stderr)
    
    # Extract tumor fraction features
    print("\nExtracting tumor fraction features...", file=sys.stderr)
    df = extract_tumor_fraction(args.ichor_dir)
    
    # Save to CSV
    df.to_csv(output_file, index=False)
    print(f"\nFeatures saved to {output_file}", file=sys.stderr)
    print("\nFinal results:", file=sys.stderr)
    print(df, file=sys.stderr)

if __name__ == '__main__':
    main()