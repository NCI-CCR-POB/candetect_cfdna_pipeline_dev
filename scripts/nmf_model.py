import argparse
import model
import numpy as np
import os
import pandas as pd
import pickle
from sklearn.decomposition import NMF
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--freq_csv", help="Path to input frequency histogram data (CSV file)")
parser.add_argument("--samplelist", help='Path to samplelist (input)', required=True)
parser.add_argument("--n_components", help='Number of nmf components (input)', type=int, required=True)
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
FREQ_CSV = args.freq_csv
SAMPLELIST = args.samplelist
N_COMPONENTS = args.n_components
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
#cohorts_df = pd.read_csv(SAMPLELIST, index_col=0)
cohorts_df = pd.read_csv(SAMPLELIST, sep="\t", index_col=0)
cohorts_df = cohorts_df.loc[cohorts_df["cohort"].isin([COHORT1,COHORT2])]

if cohorts_df.empty:
    raise ValueError(f"No samples found for cohorts: {COHORT1}, {COHORT2}")


# Enforce ordering
print(cohorts_df.columns)
cohorts_df["cohort_int"] = cohorts_df["cohort"].map({COHORT1: 0, COHORT2: 1})
cohorts_df['index_col'] = cohorts_df.index
cohorts_df = cohorts_df.sort_values(['cohort_int', 'index_col'])
cohorts_df = cohorts_df[["cohort_int"]].rename(columns={'cohort_int':'cohort'})

df = pd.read_csv(FREQ_CSV, index_col=0)
df = df.loc[cohorts_df.index]

print(cohorts_df)
print(df)

###############################################################################
# Run model
# Check if an existing model is provided
if EXISTING_MODEL and os.path.exists(EXISTING_MODEL):
    print(f"Loading existing model from {EXISTING_MODEL}")
    with open(EXISTING_MODEL, 'rb') as f:
        best_clf = pickle.load(f)
    model_output_file = f'{models_dir}/{name}.pkl'
    #Save the loaded model to the output directory for consistency, even if its the same as the existing model
    with open(model_output_file, 'wb') as f:
        pickle.dump(best_clf,f)
    print(f"Save the pre-trained model and saved: {model_output_file}")
else:
    # Handle case when existing_model is empty or doesn't exist
    print("\n\n---------------- No Existing Model Found ----------------\n")
    print(f"Training a new model with {N_COMPONENTS} components...")
    nmf = NMF(n_components=N_COMPONENTS,
              init='random',
              random_state=model.RANDOM_SEED,
              solver="mu",
              beta_loss="kullback-leibler",
              max_iter=5000)

    pipe = make_pipeline(
        nmf,
        LogisticRegression()
    )

    # Run cross-validation and obtain the best classifier
    best_clf = model.cross_validation(
        df,
        cohorts_df,
        pipe,
        models_dir,
        name,
        verbose=1,
    )

    # Save the best classifier as a pickle object
    model_output_file = f'{models_dir}/{name}.pkl'
    with open(model_output_file, 'wb') as f:
        pickle.dump(best_clf, f)
    print(f"New model trained and saved: {model_output_file}")

# Regardless of whether the model was loaded or created, get coefficients and probabilities
best_coeffs_df, prob_df = model.get_coeffs_and_probs(
    df,
    cohorts_df,
    best_clf,
)

# Save outputs
best_coeffs_df.to_csv(f"{coeffs_dir}/{name}.csv")
prob_df.to_csv(f"{probs_dir}/{name}.csv")

print("\n---------------- Model Training Complete ----------------")
print(f"Saved coefficients to: {coeffs_dir}/{name}.csv")
print(f"Saved probabilities to: {probs_dir}/{name}.csv")
