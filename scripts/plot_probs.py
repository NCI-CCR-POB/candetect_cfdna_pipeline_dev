import argparse
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
import itertools as it
import glob

def box_strip_plot(*, df, x, y, labels_order, hue, hue_order,
                   plot_name, out_dir, figsize=None):
    """Creates a boxplot with datapoints"""
    if figsize:
        plt.figure(figsize=figsize)
    ax = sns.boxplot(
        data=df,
        x=x,
        y=y,
        order=labels_order,
        hue=hue,
        hue_order=hue_order,
        boxprops={'alpha': 0.4},
        fliersize=0,
        dodge=True,
    )
    sns.stripplot(
        data=df,
        x=x,
        y=y,
        ax=ax,
        order=labels_order,
        hue=hue,
        hue_order=hue_order,
        jitter=True,
        linewidth=1,
        edgecolor='gray',
        dodge=True,
    )

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
    plt.title(plot_name)
    plt.tight_layout()
    plt.savefig(f"{out_dir}/{plot_name}.pdf")
    plt.clf()
    plt.close()

################################################################################
# Setup
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--probs_dir", help='Path to probs dir (input)', required=True)
    parser.add_argument("--samplelist", help='Path to samplelist (input)', required=True)
    parser.add_argument("--out_dir", help='Path to output directory', required=True)

    args = parser.parse_args()
    print("Running with the following argument values:")
    for k,v in sorted(vars(args).items()):
        print(f"\t{k}: {v}")

    PROBS_DIR = args.probs_dir
    SAMPLELIST = args.samplelist
    OUT_DIR = args.out_dir

    # create output directories
    prob_per_sample_dir = f"{OUT_DIR}/plots/prob_per_sample"
    prob_per_repeat_dir = f"{OUT_DIR}/plots/prob_per_repeat"
    os.makedirs(prob_per_sample_dir, exist_ok=True)
    os.makedirs(prob_per_repeat_dir, exist_ok=True)

    ############################################################################
    # Make plots for each OVO comparison
    cohorts_df = pd.read_csv(SAMPLELIST, sep='\t', header=0)
    print("cohort index")
    print(cohorts_df['library'])

    for csv in sorted(glob.glob(f"{PROBS_DIR}/*.csv")):
        # parse path to get "<c1>_v_<c2>" OR "<c1>_v_<c2>.n=<k>"
        name = csv.split("/")[-1].split(".csv")[0]
        # parse name string to get cohorts
        c1, c2 = name.split(".")[0].split("_v_")
        print(name, flush=True)

        # get probabilities df and cohorts
        probs_df = pd.read_csv(csv, index_col=0)
        print("probs df")
        print(probs_df)
        print("probs df index")
        print(probs_df.index)
        REPEAT_COLS = probs_df.columns
#        probs_df["cohort"] = cohorts_df.loc[probs_df.index]["cohort"]
        probs_df = probs_df.merge(cohorts_df[['library', 'cohort']], left_index=True, right_on='library', how='left')

        melted_df = pd.melt(probs_df.reset_index(),
                        id_vars=["library", "cohort"],
                        value_vars=REPEAT_COLS
                        )
        melted_df = melted_df.rename(
            columns={"index": "library",
                    "variable": "repeat", "value": "prob"}
        )
        print(melted_df)

        # plot probabilities per sample
        box_strip_plot(
            df=melted_df,
            x="prob",
            y="library",
            labels_order=probs_df.index.to_list().sort(),
            hue="cohort",
            hue_order=[c1, c2],
            plot_name=name,
            out_dir=prob_per_sample_dir,
            figsize=(6, 15)
        )
        # plot probabilities per repeat
        box_strip_plot(
            df=melted_df,
            x="prob",
            y="repeat",
            labels_order=REPEAT_COLS,
            hue="cohort",
            hue_order=[c1, c2],
            plot_name=name,
            out_dir=prob_per_repeat_dir,
            figsize=(6, 8)
        )
