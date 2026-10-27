import argparse
import itertools as it
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
from sklearn.metrics import roc_auc_score, roc_curve

print("Global variables at start:", globals().keys()) 

def get_optimal_threshold(y_true, y_score):
    fpr, tpr, thresh = roc_curve(y_true, y_score)
    optimal_idx = np.argmax(tpr - fpr)
    optimal_threshold = thresh[optimal_idx]
    return optimal_threshold

def get_scores_and_thresholds(probs_df, dx_df, c1, c2):
    # Align indices to ensure compatibility
    dx_df.index = dx_df.index.str.strip().str.upper()   
    probs_df.index = probs_df.index.str.strip().str.upper()
    print("probs_df index:", probs_df.index)
    print("dx_df index:", dx_df.index)

    common_indices = probs_df.index.intersection(dx_df.index)
    if common_indices.empty:
        raise KeyError("No matching indices between probs_df and dx_df.")
    dx_df = dx_df.loc[common_indices]
    # Filter dx_df to only include indices present in probs_df
    #if not probs_df.index.isin(dx_df.index).all():
    #    missing_indices = set(probs_df.index) - set(dx_df.index)
    #    raise KeyError(f"Indices missing in dx_df: {missing_indices}")
    #dx_df = dx_df.loc[probs_df.index]

    # Assign cohort
    probs_df["cohort"] = dx_df["cohort"]
    probs_df["cohort_int"] = probs_df["cohort"].map({c1: 0, c2: 1})

    # Get median optimal threshold
    optimal_threshold = get_optimal_threshold(probs_df["cohort_int"], probs_df["med_score"])
    probs_df["med_score_cutoff"] = optimal_threshold
    probs_df["med_score_pred_class"] = probs_df["med_score"].map(lambda x: 1 if x >= optimal_threshold else 0)
    probs_df["med_score_correct_pred"] = probs_df["med_score_pred_class"] == probs_df["cohort_int"]

    # Get avg optimal threshold
    #optimal_threshold = get_optimal_threshold(probs_df["cohort_int"], probs_df["avg_score"])
    #probs_df["avg_score_cutoff"] = optimal_threshold
    #probs_df["avg_score_pred_class"] = probs_df["avg_score"].map(lambda x: 1 if x >= optimal_threshold else 0)
    #probs_df["avg_score_correct_pred"] = probs_df["avg_score_pred_class"] == probs_df["cohort_int"]

    return probs_df

def get_aucs(scores_df, y_true):
    row = []
    cols = scores_df.columns.to_list()
    for col in cols:
        y_pred = scores_df[col]
        row.append(roc_auc_score(y_true, y_pred))
    return row

# Plot AUC heatmap
def plot_aucs_heatmap(df, out_dir, name):
    ax = sns.heatmap(df,
                     linewidth=0.5,
                     annot=True,
                     cbar=True,
                     square=True,
                     yticklabels=True,
                     vmin=0.5,
                     vmax=1,
                     cmap="vlag",
                     annot_kws={"fontsize": 20},
                     )
    ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=16)
    ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=16)
    ax.set_title("AUC", fontsize=24)
    ax.set_ylabel("")
    ax.set_xlabel("")
    plt.tight_layout()
    plt.savefig(f"{out_dir}/auc.{name}.pdf")
    plt.clf()
    plt.close()