import shutil
import pandas as pd
import os
import glob

SUMMARIES_DIR = "/data/panalc/mpnst_pipeline_data/results/features_summary"

SRC_DIR = "/data/panalc/ichor-data-raw/readDepth"
DST_DIR = "/data/panalc/ichor-data-6x-custom-PoN/readDepth"

for file in glob.glob(f"{SUMMARIES_DIR}/*.csv"):
    print(file)
    df = pd.read_csv(file)
    # print(df)
    libs = df["Library"].values
    print(libs)
    for lib in libs:
        src = f"{SRC_DIR}/{lib}.wig"
        dst = f"{DST_DIR}/{lib}.wig"
        # print(dst)
        shutil.copy(src, dst)
