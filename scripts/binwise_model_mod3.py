import argparse
import model
import numpy as np
import os
import pandas as pd
import pickle
import snakemake
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--armz_csv", help='Path to arm zscore data (input)', required=True)
parser.add_argument("--ratio_csv", help='Path to ratio data (input)', required=True)
parser.add_argument("--samplelist", help='Path to samplelist (input)', required=True)
parser.add_argument("--cohort1", help='First cohort in comparison (input)', required=True)
parser.add_argument("--cohort2", help='Second cohort in comparison (input)', required=True)
parser.add_argument("--out_dir", help='Path to output directory', required=True)
parser.add_argument("--existing_model", help='Path to an existing model (.pkl)', default=None)

parser.add_argument("--name", type=str, required=True, help="Name for the output files")

args = parser.parse_args()

print("Running with the following argument values:")
for k,v in sorted(vars(args).items()):
    print(f"\t{k}: {v}")

EXISTING_MODEL = args.existing_model
ARMZ_CSV = args.armz_csv
RATIO_CSV = args.ratio_csv
SAMPLELIST = args.samplelist
COHORT1 = args.cohort1
COHORT2 = args.cohort2
OUT_DIR = args.out_dir
name = args.name

# create output directories
models_dir = f"{OUT_DIR}/model"
coeffs_dir = f"{OUT_DIR}/coeffs"
probs_dir = f"{OUT_DIR}/probs"
os.makedirs(models_dir, exist_ok=True)
os.makedirs(coeffs_dir, exist_ok=True)
os.makedirs(probs_dir, exist_ok=True)

###############################################################################
# Prepare data
print("\n\n---------------- Preparing Data ----------------\n")
cohorts_df = pd.read_csv(SAMPLELIST, sep='\t', index_col=0)
#cohorts_df = cohorts_df.loc[cohorts_df["cohort"].isin([COHORT1,COHORT2])]
# enforce ordering
#cohorts_df["cohort_int"] = cohorts_df["cohort"].map({COHORT1: 0, COHORT2: 1})
cohorts_df["cohort_int"] = 0
#cohorts_df['index_col'] = cohorts_df.index
#cohorts_df = cohorts_df.sort_values(['cohort_int', 'index_col'])
cohorts_df = cohorts_df[["cohort_int"]].rename(columns={'cohort_int':'cohort'})

z_df = pd.read_csv(ARMZ_CSV, index_col=0)
ratio_df = pd.read_csv(RATIO_CSV, index_col=0)
df = pd.merge(z_df, ratio_df, how="inner", left_index=True, right_index=True)
df.sort_index(inplace=True)
df = df.loc[cohorts_df.index]

print(cohorts_df)
print(df)

###############################################################################
# Run model
# Use Existing Model or Create a New One
if EXISTING_MODEL:
    print("\n\n---------------- Using Existing Model ----------------\n")
    with open(EXISTING_MODEL, 'rb') as f:
        best_clf = pickle.load(f)
else:
    print("\n\n---------------- Training New Model ----------------\n")
    name = f"{COHORT1}_v_{COHORT2}"

    pca_pipe = make_pipeline(
        StandardScaler(),
        PCA(n_components=0.9, random_state=model.RANDOM_SEED if hasattr(model, 'RANDOM_SEED') else 42),
    )
    col_trans = ColumnTransformer(
        [("zscore", "passthrough", z_df.columns.to_list()),
         ("pcs", pca_pipe, ratio_df.columns.to_list())],
    )
    pipe = make_pipeline(
        col_trans,
        LogisticRegression(),
    )

    # Run cross-validation
    best_clf = model.cross_validation(
        df,
        cohorts_df,
        pipe,
        models_dir,
        name,
        verbose=1,
    )
    # Save best classifier as a pickle object
    with open(f'{models_dir}/{name}.pkl', 'wb') as f:
        pickle.dump(best_clf, f)

# Re-run best classifier CV to get scores per fold
#best_coeffs_df, prob_df = model.get_coeffs_and_probs(df, cohorts_df, best_clf)
#best_coeffs_df.to_csv(f"{coeffs_dir}/{name}.csv")
#prob_df.to_csv(f"{probs_dir}/{name}.csv")

# Get scores for ALL samples at once
scores = best_clf.predict_proba(df)[:, 1]

# Save to DataFrame
result_df = pd.DataFrame({
    "score": scores
}, index=df.index)

print(result_df)
result_df .to_csv(f"{probs_dir}/{name}.csv")
