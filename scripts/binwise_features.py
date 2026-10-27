import argparse
import numpy as np
import os
import pandas as pd
import re
from sklearn.preprocessing import StandardScaler

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--armz_csv", help='Path to bin wise fragmentomic arm-level z-scores .csv (input)', required=True)
parser.add_argument("--ratio_tsv", help='Path to ratio .tsv (input)', required=True)
parser.add_argument("--out_dir", help='Path to output directory', required=True)

args = parser.parse_args()
BINWISE_ZSCORE_CSV = args.armz_csv
BINWISE_RATIO_TSV = args.ratio_tsv
OUT_DIR = args.out_dir

os.makedirs(OUT_DIR,exist_ok=True)

################################################################################
# Helper
################################################################################
def merge_index_columns(df, index_cols, join_str="_", remove_chr=False):
    # Takes a pandas dataframe with no index column
    # Merges specified columns with "_" to create a single index column
    # Drops merged columns and sets computed column as index
    df["index"] = df[index_cols].apply(lambda row: join_str.join(row.values.astype(str)), axis=1)   
    if remove_chr:
        # df["index"] = df["index"].str.removeprefix("chr")
        # "removeprefix()" is not supported prior to python 3.9+
        # XXX: this removes all instances of each char is "chr"
        # df["index"] = df["index"].str.lstrip("chr")
        df["index"] = df["index"].str.replace("chr","")

    df.drop(index_cols, axis=1, inplace=True)
    df.set_index("index", inplace=True)

def sort_by_chromosome(df, sort_cols):
    # Reorder df to be chr order (non-lexicographic)
    # sort_cols specifies secondary columns to sort on (e.g. arm or start coord)
    
    # need to handle M (mitochondrial), X, and Y specially
    df["chr"] = df["chr"].str.replace("M","97")
    df["chr"] = df["chr"].str.replace("X","98")
    df["chr"] = df["chr"].str.replace("Y","99")
    df["chr_int"] = pd.to_numeric(df["chr"].str.replace("chr",""))
    df.sort_values(by=["chr_int"] + sort_cols,inplace=True)
    df["chr"] = df["chr"].str.replace("97","M")
    df["chr"] = df["chr"].str.replace("98","X")
    df["chr"] = df["chr"].str.replace("99","Y")
    df.drop(columns=["chr_int"],inplace=True)

def clean_df(df):
    # drop completely null rows (libs)
    df = df.dropna(axis=0, how="all")
    # drop cols with any null entries (features)
    df = df.dropna(axis=1, how="any")
    # drop non-autosomal coords (X,Y,M)
    for chrom in ['X','Y','M']:
        df = df.loc[:,~df.columns.str.startswith(chrom)]
    return df

def standard_scaling_rows(df):
    x_array = np.array(df)
    print(x_array)
    sc = StandardScaler()
    x_mplx_sc = sc.fit_transform(x_array.T)
    return x_mplx_sc.T

def arm_key(arm_str):
    arm_int, arm_letter = arm_str.split("_")
    arm_int = int(arm_int[3:])
    return int(arm_int), arm_letter

################################################################################
# Bin wise fragmentomics z-score
################################################################################
# This is already mostly done from the binwise_armz.R script. Just need to order
# arms numerically and clean any empty rows (should be none)
print("Processing bin wise fragmentomic z-score data")

zscores_df = pd.read_csv(BINWISE_ZSCORE_CSV, index_col=0)
cols = zscores_df.columns
sorted_cols = sorted(cols, key=arm_key)
zscores_df = zscores_df[sorted_cols]
print(zscores_df)
zscores_df.to_csv(OUT_DIR + "/binwise_armz.sorted.csv")

zscores_df = clean_df(zscores_df)
print(zscores_df)
zscores_df.to_csv(OUT_DIR + "/binwise_armz.clean.csv")

print("Finished processing bin wise fragmentomic z-score data")

################################################################################
# Bin-wise fragmentomic ratio
################################################################################
print("Processing bin wise fragmentomic ratio data")

binwise_df = pd.read_csv(BINWISE_RATIO_TSV, sep="\t")
binwise_df.columns = ["library","chr","start","end","fract","ratio.centered"]
lib_ids = binwise_df["library"].unique()

# Treat ["chr","start","end"] as keys to merge each lib's fract as cols
index_cols = ["chr","start","end"]
ratios_df = pd.DataFrame(columns=index_cols)
for lib_id in lib_ids:
    print(lib_id)
    lib_df = binwise_df.loc[(binwise_df['library'] == lib_id)]
    ratios_df = pd.merge(ratios_df, lib_df[index_cols + ["fract"]], how="outer", on=index_cols)
    ratios_df.rename(columns={"fract":lib_id}, inplace=True)

sort_by_chromosome(ratios_df, ["start"])
merge_index_columns(ratios_df, index_cols, join_str="_", remove_chr=True)
ratios_df = ratios_df.T # set libs as rows, bins as cols
print(ratios_df)
ratios_df.to_csv(OUT_DIR + "/binwise_ratio.raw.csv")

# NOTE: this drops 9/322 columns.
ratios_df = clean_df(ratios_df)
print(ratios_df)
ratios_df.to_csv(OUT_DIR + "/binwise_ratio.clean.csv")

## checking scikit normalization - yields same result as StandardScaler
# from scipy import stats
# z_df = stats.zscore(np.array(ratios_df),axis=1)
# z_df = pd.DataFrame(z_df,                            
#                     columns=ratios_df.columns.to_list(),
#                     index=ratios_df.index.to_list())
# print(z_df)
# print(z_df.mean(axis=1))
# print(z_df.std(axis=1))

# z-score normalization over rows (i.e. over libs)
ratios_df = pd.DataFrame(standard_scaling_rows(ratios_df),
                            columns=ratios_df.columns.to_list(),
                            index=ratios_df.index.to_list())
print(ratios_df)
# NOTE: will not be exactly 0 mean 1 sd. Pretty sure this is a fp precision
# issue, but values will be pretty close.
print(ratios_df.mean(axis=1)) # check row mean is 0
print(ratios_df.std(axis=1)) # check row std-dev is 1
ratios_df.to_csv(OUT_DIR + "/binwise_ratio.norm_rows.clean.csv")

print("Finished processing bin-wise fragmentomic ratio data")
