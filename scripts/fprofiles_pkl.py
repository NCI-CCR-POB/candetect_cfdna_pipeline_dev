#!/usr/bin/env python3
"""
Script to perform Non-negative Matrix Factorization (NMF) on motif frequency data from a TSV file,
followed by deconvolution analysis using Non-Negative Least Squares (NNLS).
"""

import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from scipy.optimize import nnls
import logging
import argparse
import pickle
import os

# --- Load Arguments --- #
def load_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data_file", type=str, required=True, help="Path to the TSV file containing motif frequencies")
    parser.add_argument("--output_fprof_per_lib", type=str, required=True, help="Path to output file for F-profiles per library")
    parser.add_argument("--output_motif_per_fprof", type=str, required=True, help="Path to output file for normalized motif F-profiles")
    parser.add_argument("--save_model", type=str, default=None, help="Path to save the F-profile NMF model")
    parser.add_argument("--load_model", type=str, default=None, help="Path to load a pre-trained F-profile NMF model")
    return parser.parse_args()

# --- Helper Functions --- #
def load_data(filepath):
    """Load motif data from a TSV file."""
    try:
        # Read the data, keeping the first column as index
        df = pd.read_csv(filepath, delimiter='\t', index_col=0)
        
        # Check if the data is empty
        if df.empty:
            raise ValueError("Input data file is empty")
            
        # Convert all columns to numeric, coercing errors to NaN
        df = df.apply(pd.to_numeric, errors='coerce')
        
        # Check for any NaN values
        if df.isna().any().any():
            logging.warning("Found NaN values in the data. These will be replaced with 0.")
            df = df.fillna(0)
            
        return df
    except Exception as e:
        logging.error(f"Error loading data from {filepath}: {str(e)}")
        raise

def perform_nmf(data_matrix, n_components=6):
    """Perform NMF on the provided data matrix."""
    try:
        model = NMF(n_components=n_components, init='random', random_state=0, max_iter=100000)
        W = model.fit_transform(data_matrix)
        F = model.components_
        return W, F
    except Exception as e:
        logging.error(f"Error performing NMF: {str(e)}")
        raise

def save_nmf_model(W, F, filepath):
    """Save the NMF components (W and F) to a pickle file."""
    with open(filepath, 'wb') as f:
        pickle.dump({'W': W, 'F': F}, f)
    logging.info(f"Model components saved to {filepath}.")

def load_nmf_model(filepath):
    """Load the NMF components (W and F) from a pickle file."""
    with open(filepath, 'rb') as f:
        model = pickle.load(f)
    logging.info(f"Model components loaded from {filepath}.")
    print(type(model))
    print(model)
    #print(model.named_steps)
    return model['W'], model['F']

def perform_nnls(F, sample_motif_frequencies):
    """Perform NNLS to decompose sample motif frequencies into F-profile contributions."""
    P, _ = nnls(F.T, sample_motif_frequencies)
    return 100 * P / np.sum(P)

# --- Main Function --- #
def main():
    logging.basicConfig(level=logging.INFO)
    args = load_arguments()

    try:
        logging.info("Loading motif frequency data...")
        motifs_df = load_data(args.data_file)
        
        # Transpose the data for NMF
        M = motifs_df.T
        
        logging.info(f"Data loaded with shape {M.shape}")
        
        if args.load_model and os.path.exists(args.load_model):
            logging.info("Loading pre-trained NMF F-profile model...")
            W, F = load_nmf_model(args.load_model)
        else:
            logging.info("Performing NMF for F-profiles...")
            W, F = perform_nmf(M, n_components=6)
            
            if args.save_model:
                save_nmf_model(W, F, args.save_model)

        # Calculate NNLS contributions
        nnls_contributions_df = pd.DataFrame()
        for lib in motifs_df.columns:
            logging.info(f"Processing sample {lib}...")
            sample_motif_frequencies = motifs_df[lib].values
            P_normalized = perform_nnls(F, sample_motif_frequencies)
            nnls_contributions_df[lib] = P_normalized

        nnls_contributions_df.index = [f'fprof{i+1}' for i in range(F.shape[0])]
        nnls_contributions_df.to_csv(args.output_fprof_per_lib, sep='\t')
        logging.info(f"NNLS contributions saved to {args.output_fprof_per_lib}.")

        # Normalize and save F-profiles
        F_normalized = F / F.sum(axis=1, keepdims=True)
        F_df_normalized = pd.DataFrame(F_normalized, columns=motifs_df.index, index=[f'fprof{i+1}' for i in range(F.shape[0])])
        F_df_normalized.to_csv(args.output_motif_per_fprof, sep='\t', index=True)
        logging.info(f"Normalized motif profiles saved to {args.output_motif_per_fprof}.")

    except Exception as e:
        logging.error(f"Error in main function: {str(e)}")
        raise

# --- Main Guard --- #
if __name__ == "__main__":
    main()
