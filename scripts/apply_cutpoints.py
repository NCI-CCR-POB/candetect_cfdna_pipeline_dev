#!/usr/bin/env python3
"""
Apply a pre-established Youden J cutpoint to new sample probabilities.
Mode 3: existing model + existing cutpoint — no new AUC computed.

Reads:
  --probs_csv      : probs CSV from model application (repeat_0..repeat_N cols)
  --scores_csv     : original run's scores CSV containing med_score_cutoff
  --samplelist     : sample list TSV (library, cohort columns)
  --cohort1/2      : cohort names
  --out_csv        : output path

Output columns: repeat_0..N, med_score, med_score_cutoff, med_score_pred_class,
                 cohort, cohort_int
"""

import argparse
import pandas as pd
import sys

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--probs_csv",  required=True, help="New sample probs CSV")
parser.add_argument("--scores_csv", required=True, help="Original scores CSV with med_score_cutoff")
parser.add_argument("--samplelist", required=True, help="Sample list TSV")
parser.add_argument("--cohort1",    required=True, help="First cohort (label 0)")
parser.add_argument("--cohort2",    required=True, help="Second cohort (label 1)")
parser.add_argument("--out_csv",    required=True, help="Output scores CSV path")
args = parser.parse_args()

print("Running apply_cutpoints.py with:")
for k, v in sorted(vars(args).items()):
    print(f"  {k}: {v}")

# Load new sample probs and compute med_score
probs_df = pd.read_csv(args.probs_csv, index_col=0)
#no repleat when using existing model and cutpoint
#repeat_cols = [c for c in probs_df.columns if c.startswith("repeat_")]
#if not repeat_cols:
#    sys.exit(f"ERROR: No repeat_* columns found in {args.probs_csv}")
#probs_df["med_score"] = probs_df[repeat_cols].median(axis=1)

# Load cutpoint from the original run's scores file
scores_df = pd.read_csv(args.scores_csv, index_col=0)
if "med_score_cutoff" not in scores_df.columns:
    sys.exit(f"ERROR: med_score_cutoff column not found in {args.scores_csv}")
cutoff_values = scores_df["med_score_cutoff"].dropna()
if cutoff_values.empty:
    sys.exit(f"ERROR: med_score_cutoff is entirely null in {args.scores_csv}")
cutoff = float(cutoff_values.iloc[0])
print(f"Applying cutpoint from original run: {cutoff:.6f}")

# Apply cutpoint
probs_df["med_score_cutoff"] = cutoff
#probs_df["med_score_pred_class"] = (probs_df["med_score"] >= cutoff).astype(int)
probs_df["med_score_pred_class"] = (probs_df["score"] >= cutoff).astype(int)

# Annotate with cohort info
dx_df = pd.read_csv(args.samplelist, sep="\t", index_col=0)
dx_df.index = dx_df.index.astype(str).str.strip().str.upper()
probs_df.index = probs_df.index.astype(str).str.strip().str.upper()
common = probs_df.index.intersection(dx_df.index)
if common.empty:
    print("WARNING: No matching indices between probs and samplelist — cohort columns will be NaN")
probs_df["cohort"] = dx_df.loc[common, "cohort"] if not common.empty else None
#probs_df["cohort_int"] = probs_df["cohort"].map({args.cohort1: 0, args.cohort2: 1})
probs_df["cohort_int"] = probs_df["cohort"]

probs_df.to_csv(args.out_csv)
print(f"Scores written to: {args.out_csv}")

# Summary
n_cohort1 = (probs_df["med_score_pred_class"] == 0).sum()
n_cohort2 = (probs_df["med_score_pred_class"] == 1).sum()
print(f"Predicted {args.cohort1}: {n_cohort1}, {args.cohort2}: {n_cohort2}")
