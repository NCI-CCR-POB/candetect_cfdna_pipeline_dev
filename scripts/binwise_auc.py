import argparse
import aucs
import itertools as it
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import pandas as pd
import seaborn as sns
import numpy as np

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--probs_dir", help='Path to end motifs probs dir (input)', required=True)
parser.add_argument("--samplelist", help='Path to samplelist (input)', required=True)
parser.add_argument("--out_dir", help='Path to output directory', required=True)
parser.add_argument("--cohorts", help='Disease cohorts', required=True)

args = parser.parse_args()
print("Running with the following argument values:")
for k,v in sorted(vars(args).items()):
    print(f"\t{k}: {v}")

PROBS_DIR = args.probs_dir
SAMPLELIST = args.samplelist
OUT_DIR = args.out_dir
COHORTS = args.cohorts.split(',')

print(globals().keys())

print(f"Loading dx_df from samplelist: {SAMPLELIST}")
try:
    dx_df = pd.read_csv(SAMPLELIST, sep='\t', header=0)
except Exception as e:
    raise RuntimeError(f"Error loading samplelist file: {e}")
required_columns = ['library', 'cohort'] # Retain only `library` and `cohort` columns
missing_columns = [col for col in required_columns if col not in dx_df.columns]
if missing_columns:
    raise ValueError(f"Missing required columns in samplelist file: {missing_columns}")
dx_df = dx_df[required_columns]
dx_df = dx_df.apply(lambda col: col.map(lambda x: x.strip() if isinstance(x, str) else x))
print("dx_df:")
print(dx_df.head())

COHORT_1 = []
COHORT_2 = []
for c1, c2 in it.combinations(COHORTS, 2):
    COHORT_1.append(c1)
    COHORT_2.append(c2)
    print(c1, c2)

# create output directories
scores_dir = f"{OUT_DIR}/scores"
aucs_dir = f"{OUT_DIR}/aucs/aucs"  # Changed to match expected path
aucs_plots_dir = f"{OUT_DIR}/plots/aucs"
aucs_v_comp_plots_dir = f"{OUT_DIR}/plots/aucs_v_components"
os.makedirs(scores_dir, exist_ok=True)
os.makedirs(aucs_dir, exist_ok=True)
os.makedirs(aucs_plots_dir, exist_ok=True)
os.makedirs(aucs_v_comp_plots_dir, exist_ok=True)

cohorts_df = pd.read_csv(SAMPLELIST, sep='\t', index_col=0)

################################################################################
# Compute final score and thresholds for each comparison
# Summarize scores across all comparisons

best_rows = []
rows = []
for c1, c2 in it.combinations(COHORTS, 2):
    print("c1:", c1, "c2:", c2)
    name = f"{c1}_v_{c2}"
    row = [c1, c2]
        
    # get scores df
    probs_df = pd.read_csv(
        f"{PROBS_DIR}/{name}.csv", index_col=0)
    
    # Calculate median score from repeat columns
    repeat_cols = [col for col in probs_df.columns if col.startswith('repeat_')]
    probs_df['med_score'] = probs_df[repeat_cols].median(axis=1)
    
    scores_df = aucs.get_scores_and_thresholds(probs_df,cohorts_df,c1,c2)
    print(scores_df)
    scores_df.to_csv(f"{scores_dir}/{name}.csv")

    # get aucs for each fold, avg_score, med_score
    y_true = scores_df["cohort_int"]
    print(y_true)
    scores_df = scores_df.iloc[:,0:11]
    row += aucs.get_aucs(scores_df, y_true)
    rows.append(row)

aucs_df = pd.DataFrame(data=rows,
                    columns=["cohort1", "cohort2"]+scores_df.columns.to_list())
print(aucs_df)
aucs_df.to_csv(f"{aucs_dir}/aucs.csv")

############################################################################
# Plot AUC summaries
df = aucs_df.pivot(index="cohort1", columns="cohort2",
                values="med_score")
df = df.reindex(COHORTS[:-1])
df = df[COHORTS[1:]]
aucs.plot_aucs_heatmap(df,aucs_plots_dir,"med_score")