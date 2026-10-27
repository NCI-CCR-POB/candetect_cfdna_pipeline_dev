import argparse
import glob
import os
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--hist_dir',
    help='Directory with picard InsertSizeMetrics histograms',
    required=True,
    )
parser.add_argument('-o','--output',
    help='Output dir',
    required=True,
    )
parser.add_argument('-v','--verbose',
    help='Verbose output',
    action='store_true',
    required=False,
    )

args = parser.parse_args()
print("Running with the following argument values:")
for k,v in sorted(vars(args).items()):
    print("\t{}: {}".format(k,v))
print("\n")

HIST_DIR = args.hist_dir
OUT_DIR  = args.output
VERBOSE  = args.verbose

if not os.path.exists(HIST_DIR):
    raise Exception("HIST_DIR was not found")
if not os.path.exists(OUT_DIR):
    raise Exception("OUT_DIR was not found")

# Set initial range
START = 30
END = 700

###############################################################################
# Expected format of filenames in histogram dir is "<lib>_metrics.txt"
#
# picard InsertSizeMetrics files have a 10 line header and then 
# two columns of data: "insert_size" and "All_Reads.fr_count"
#
# Individual reports are merged into a single dataframe where rows are
# libs and columns are insert sizes.

metrics = sorted(glob.glob(HIST_DIR + "/*_metrics.txt"))
df = pd.DataFrame(columns=['insert_size'])

print("Files to process:",len(metrics))
i = 0
for metric in metrics:
    if VERBOSE:
        print("Reading {}".format(metric))
    #name = metric.split("/")[-1].split("_")[0]
    name = metric.split("/")[-1].replace("_metrics.txt", "")
    print(f"Processing file: {metric}")
    print(f"Extracted name: {name}")
    metrics_df = pd.read_csv(metric, sep="\t", skiprows=10)
    metrics_df = metrics_df.rename(columns={"All_Reads.fr_count":name})
    df = pd.merge(df,metrics_df, how="outer", on="insert_size")
    
    i += 1
    if i % 20 == 0 and VERBOSE:
        print(i, "files completed")

df = df.set_index("insert_size")
df = df.sort_values(by=["insert_size"])
df = df.fillna(0)
print(df)

# Get the actual min and max row indices
min_row = df.index.min()
max_row = df.index.max()

# Adjust START and END based on actual data
actual_start = max(START, min_row)
actual_end = min(END, max_row)

# Create the range based on actual data
RANGE = [i for i in range(actual_start, actual_end + 1)]

# Write the full histogram
df.T.to_csv(f"{OUT_DIR}/count_histogram.csv")

# Write the ranged histogram
df_ranged = df.loc[RANGE]
print(df_ranged)

# Calculate and write frequencies
df_freq = df_ranged.apply(lambda x: x/x.sum(),axis=0)
print(df_freq)
print(df_freq.sum(axis=0))
df_freq.T.to_csv(f"{OUT_DIR}/freq_histogram.{START}_{END}.csv")
