import argparse
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
import itertools as it
import glob
import sys

def plot_coefficient_heatmap(best_coeffs_df, name, out_dir):
    ########################################################
    # heatmap of coefficients in each training fold
    try:
        # some magic for scaling the heatmap
        # https://stackoverflow.com/questions/35127920/overlapping-yticklabels-is-it-possible-to-control-cell-size-of-heatmap-in-seabo
        # get the tick label font size
        fontsize_pt = 6
        plt.yticks(fontsize=fontsize_pt)
        plt.xticks(fontsize=fontsize_pt)
        dpi = 72.27

        # compute the matrix height in points and inches
        matrix_height_pt = fontsize_pt * best_coeffs_df.shape[0]
        matrix_height_in = matrix_height_pt / dpi

        # compute the required figure height
        top_margin = 0.04  # in percentage of the figure height
        bottom_margin = 0.12  # in percentage of the figure height
        figure_height = matrix_height_in / (1 - top_margin - bottom_margin)
        if figure_height < 4:
            figure_height = 4

        coeff_min = abs(best_coeffs_df.min().min())
        coeff_max = abs(best_coeffs_df.max().max())
        max_bound = max([coeff_min,coeff_max])

        # build the figure instance with the desired height
        fig, ax = plt.subplots(
            figsize=(8, figure_height),
            gridspec_kw=dict(top=1-top_margin, bottom=bottom_margin))

        # plt.figure(figsize=(8,18))
        sns.set(font_scale = 0.7)
        ax = sns.heatmap(best_coeffs_df,
                         ax=ax,
                         cmap="vlag",
                         center=0,
                         vmin=-max_bound,
                         vmax=max_bound,
                         yticklabels=True,
                         square=True,
                         annot_kws={"fontsize":fontsize_pt},
                         )
        plt.xlabel("Cross-Validation Fold")
        plt.ylabel("Feature")
        plt.title(name)
        
        # Create all parent directories
        os.makedirs(os.path.dirname(f"{out_dir}/{name}.pdf"), exist_ok=True)
        plt.savefig(f"{out_dir}/{name}.pdf")
        plt.clf()
        plt.close()
        print(f"Successfully saved plot to {out_dir}/{name}.pdf")
    except Exception as e:
        print(f"Error in plot_coefficient_heatmap: {str(e)}", file=sys.stderr)
        raise

################################################################################
# Setup
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--data_dir", help='Path to data dir (input)', required=True)
parser.add_argument("--armz_csv", help='Path to arm zscore data (input)', required=True)

args = parser.parse_args()
print("Running with the following argument values:")
for k,v in sorted(vars(args).items()):
    print(f"\t{k}: {v}")

DATA_DIR = args.data_dir
ARMZ_CSV = args.armz_csv

FEATURES = ["binwise"]

############################################################################
# Plot coefficient heatmaps across each OVO comparison
try:
    print(f"Reading armz data from {ARMZ_CSV}")
    armz_df = pd.read_csv(ARMZ_CSV, index_col=0)
    arm_names = armz_df.columns.to_list()
    print(f"Found {len(arm_names)} arm names")

    for feature in FEATURES:
        # Update paths to match snakemake workflow structure
        coeffs_dir = f"{DATA_DIR}/model/coeffs"
        coeff_heatmap_dir = f"{DATA_DIR}/binwise/plots/coeff_heatmap"
        
        print(f"Looking for coefficient files in {coeffs_dir}")
        csv_files = sorted(glob.glob(f"{coeffs_dir}/*.csv"))
        print(f"Found {len(csv_files)} coefficient files: {csv_files}")

        for csv in csv_files:
            # parse path to get "<c1>_v_<c2>" OR "<c1>_v_<c2>.n=<k>"
            name = "binwise_armz.clean"  # Use the expected output name
            print(f"Processing {name}", flush=True)

            print(f"Reading coefficients from {csv}")
            coeffs_df = pd.read_csv(csv, index_col=0)
            print(f"Coefficient dataframe shape: {coeffs_df.shape}")

            # label features
            if feature == "binwise":
                pcs = ["pc" + str(i) for i in range(coeffs_df.shape[0] - len(arm_names))]
                coeffs_df.set_index(pd.Index(arm_names+pcs), inplace=True)
            elif feature == "motifs":
                coeffs_df.set_index(pd.Index(motif_names), inplace=True)
            print("Coefficient dataframe:")
            print(coeffs_df)

            plot_coefficient_heatmap(coeffs_df, name, coeff_heatmap_dir)

except Exception as e:
    print(f"Error in main: {str(e)}", file=sys.stderr)
    raise
