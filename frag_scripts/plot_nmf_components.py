#!/usr/bin/env python3
"""
Script to plot overlapping histograms of NMF components.
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
import seaborn as sns

def load_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=str, required=True, help="Path to the NMF model pickle file")
    parser.add_argument("--freq_csv", type=str, required=True, help="Path to the frequency histogram CSV file")
    parser.add_argument("--output", type=str, required=True, help="Path to output plot file")
    parser.add_argument("--cohort1", type=str, required=True, help="First cohort name")
    parser.add_argument("--cohort2", type=str, required=True, help="Second cohort name")
    parser.add_argument("--n_components", type=int, required=True, help="Number of NMF components")
    return parser.parse_args()

def plot_nmf_components(model_path, freq_csv, output_path, cohort1, cohort2, n_components):
    # Load the NMF model
    with open(model_path, 'rb') as f:
        model = pickle.load(f)
    
    # Load the frequency data
    freq_df = pd.read_csv(freq_csv, index_col=0)
    
    # Get the components from the model
    components = model.named_steps['nmf'].components_
    
    # Create a single figure
    plt.figure(figsize=(12, 8))
    
    # Set style
    plt.style.use('seaborn')
    
    # Get the fragment sizes from the frequency data
    fragment_sizes = freq_df.columns.astype(int)
    
    # Create a color palette
    colors = plt.cm.Set1(np.linspace(0, 1, n_components))
    
    # Plot each component
    for i, component in enumerate(components):
        plt.plot(fragment_sizes, component, 
                label=f'Component {i+1}',
                color=colors[i],
                linewidth=2)
    
    # Customize the plot
    plt.title('NMF Components', fontsize=14, pad=20)
    plt.xlabel('Fragment Size (bp)', fontsize=12)
    plt.ylabel('Component Weight', fontsize=12)
    plt.grid(True, alpha=0.3)
    
    # Add legend
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    args = load_arguments()
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    # Generate the plot
    plot_nmf_components(
        args.model,
        args.freq_csv,
        args.output,
        args.cohort1,
        args.cohort2,
        args.n_components
    )

if __name__ == "__main__":
    main() 