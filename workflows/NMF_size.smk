#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#     NMF size Fragmentomic Analysis of Cell-free DNA Whole Genome Sequencing  #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8

import pandas as pd
import re
import numpy as np
import glob
import os
import itertools as it
import yaml


#######define software path in this part#########
current_dir = os.getcwd()
sys.path.append(f'{current_dir}/config')
import program

# -----------------------------
# Path expansion for strings in config (~, $VARS)
# -----------------------------
def resolve_config_paths(config_dict):
    for k, v in config_dict.items():
        if isinstance(v, str):
            config_dict[k] = os.path.expandvars(os.path.expanduser(v))
        elif isinstance(v, dict):
            resolve_config_paths(v)
        elif isinstance(v, list):
            config_dict[k] = [os.path.expandvars(os.path.expanduser(i)) if isinstance(i, str) else i for i in v]

resolve_config_paths(config)

print(config)


#########################################################
# READ CLUSTER PER-RULE REQUIREMENTS
#########################################################
CLUSTERYAML = os.path.join(current_dir,"config/cluster.yaml")
with open(CLUSTERYAML) as json_file:
    CLUSTER = yaml.safe_load(json_file)

## Create lambda functions to allow a way to insert read-in values
## as rule directives
getthreads=lambda rname:int(CLUSTER[rname]["threads"]) if rname in CLUSTER and "threads" in CLUSTER[rname] else int(CLUSTER["__default__"]["threads"])
getmemg=lambda rname:CLUSTER[rname]["mem"] if rname in CLUSTER and "mem" in CLUSTER[rname] else CLUSTER["__default__"]["mem"]
getmemG=lambda rname:getmemg(rname).replace("g","G")


# Values directly from configuration file
bwa_dir             = config["ref_dir"]
COHORTS             = config["cohorts"]
data_dir            = config["data_dir"]
fasta_base          = config["fasta_base"]
FRAG_DISTROS        = config["frag_distro"]
frag_repo           = config["frag_repo"]
frag_threads        = config["threads"]
genome_fasta        = config["genome_fasta"]
metadata_file       = config["metadata_file"]
normal_control      = config.get("normal_control", "healthy")  # Default to "healthy" if not defined
PE1                 = config.get("PE1", "need to define PE1")
PE2                 = config.get("PE2", "need to define PE2")
refdir              = config["ref_dir"]
samples             = config["samples"]
threads             = config["threads"]
ENV_FRAGENTOMICS    = config["env_fragentomics"]

# Directory values derived from configuration YAML
benchdir            = data_dir + "/benchmark"
frag                = data_dir + "/analysis/cfdna_wgs"
frag_bams           = data_dir  + "/bams"
frag_fastqs         = data_dir + "/fastq"
frag_scriptdir      = frag_repo +  "/scripts"
logdir              = data_dir + "/logs/NMF_size"
qcdir               = data_dir + "/qc"

# NMF specific inputs
min_components      = config["min_nmf_components"]
max_components      = config["max_nmf_components"]
nc                  = range(config["min_nmf_components"],config["max_nmf_components"]+1)

# Mode 3: apply existing model + existing cutpoints (no new AUC)
APPLY_NMF_CUTPOINTS = config.get("existing_nmf_scores_dir", "no") != "no"

# NMF specific directories
nmf_dir             =  frag + "/nmf"

nmf_agg_hist_dir    = nmf_dir + "/agg"
nmf_comps           = [f"{nmf_dir}/model/{c1}_v_{c2}.n={{nc}}.pkl" for (c1,c2) in it.combinations(COHORTS,2)]
nmf_probs           = [f"{nmf_dir}/probs/{c1}_v_{c2}.n={{nc}}.csv" for (c1,c2) in it.combinations(COHORTS,2)]
nmf_raw_hist_dir    = nmf_dir + "/raw"
nmf_scores          = [f"{nmf_dir}/scores/{c1}_v_{c2}.n={{nc}}.csv" for (c1,c2) in it.combinations(COHORTS,2)]

####------ Read samples ------####
frag_libs = pd.read_table(config["samples"])

# Ensure readable fastqs
readable = []
for x in frag_libs.file:
    readable.append(os.access(x, os.R_OK))
frag_libs['readable']=readable
frag_libs = frag_libs[frag_libs.readable == True]

# Make the dictionary
FRAG_LIBS = frag_libs["library"].tolist()
print("FRAG_LIBS:", FRAG_LIBS)
frag_libs_file_indict = frag_libs["file"].tolist()
frag_lib_dict = dict(zip(FRAG_LIBS, frag_libs_file_indict))

#Check that all BAM files exist
print("Checking BAM files:")
for lib in FRAG_LIBS:
    bam_path = f"{frag_bams}/{lib}_filt.bam"
    print(f"  {bam_path} exists: {os.path.exists(bam_path)}")

# Make  a list of healthy libraries
FRAG_HEALTHY_LIBRARIES = frag_libs[frag_libs['cohort'] == normal_control]['library'].tolist()

pairs = []
pairs = list(it.combinations(COHORTS, 2))
print(pairs)
print([pair[0] for pair in pairs])

# Mode 3: determine best n_components per cohort pair from the original run at parse time
if APPLY_NMF_CUTPOINTS:
    _existing_nmf_dir = config["existing_nmf_scores_dir"]
    _best_aucs_path = os.path.join(_existing_nmf_dir, "aucs", "best_med_aucs.csv")
    if not os.path.exists(_best_aucs_path):
        raise ValueError(
            f"Mode 3 (existing_nmf_scores_dir set) requires best_med_aucs.csv but not found: {_best_aucs_path}"
        )
    _best_aucs_df = pd.read_csv(_best_aucs_path)
    BEST_NC = {(row["cohort1"], row["cohort2"]): int(row["n_components"])
               for _, row in _best_aucs_df.iterrows()}
    BEST_NMF_FNAMES = [f"{c1}_v_{c2}" for c1, c2 in pairs]
    BEST_NMF_NC    = [BEST_NC.get((c1, c2), min_components) for c1, c2 in pairs]
    print(f"Mode 3 NMF — best n_components per pair: {dict(zip(BEST_NMF_FNAMES, BEST_NMF_NC))}")

#for c1, c2 in it.combinations(COHORTS, 2):
#    pairs.append((c1, c2))
#    print(f"{c1} vs {c2}") 

def get_cohort_pairs():
    return pairs

#Define cohort comparisons (one-versus-one)
#COHORT_1 = []
#COHORT_2 = []
#for c1, c2 in it.combinations(COHORTS,2):
#    COHORT_1.append(c1)
#    COHORT_2.append(c2)
#    print(c1,c2)

#def get_cohort_pairs():
#    pairs = []
#    for c1, c2 in it.combinations(COHORTS, 2):
#        pairs.append((c1, c2))
#    return pairs

# Shared histogram outputs (all modes)
_nmf_hist_outputs = [
    expand(nmf_raw_hist_dir + "/{lib}_metrics.txt", lib=FRAG_LIBS),
    expand(nmf_raw_hist_dir + "/{lib}_histogram.pdf", lib=FRAG_LIBS),
    nmf_agg_hist_dir + "/count_histogram.csv",
    nmf_agg_hist_dir + "/freq_histogram.30_700.csv",
]

# Mode 3: only best n_components per cohort pair, scores via cutpoint application
if APPLY_NMF_CUTPOINTS:
    _nmf_analysis_outputs = [
        expand(nmf_dir + "/probs/{fnames}.n={nc}.csv",  zip, fnames=BEST_NMF_FNAMES, nc=BEST_NMF_NC),
        expand(nmf_dir + "/model/{fnames}.n={nc}.pkl",  zip, fnames=BEST_NMF_FNAMES, nc=BEST_NMF_NC),
        #expand(nmf_dir + "/coeffs/{fnames}.n={nc}.csv", zip, fnames=BEST_NMF_FNAMES, nc=BEST_NMF_NC),
        expand(nmf_dir + "/scores/{fnames}.n={nc}.csv", zip, fnames=BEST_NMF_FNAMES, nc=BEST_NMF_NC),
    ]
# Modes 1/2: full component sweep + AUC + plots
else:
    _all_fnames = [f"{c1}_v_{c2}" for c1, c2 in pairs]
    _all_nc     = range(min_components, max_components + 1)
    _nmf_analysis_outputs = [
        expand(nmf_dir + "/probs/{fnames}.n={nc}.csv",           fnames=_all_fnames, nc=_all_nc),
        expand(nmf_dir + "/model/{fnames}.n={nc}.pkl",           fnames=_all_fnames, nc=_all_nc),
        expand(nmf_dir + "/coeffs/{fnames}.n={nc}.csv",          fnames=_all_fnames, nc=_all_nc),
        expand(nmf_dir + "/plots/components/{fnames}.n={nc}.pdf", fnames=_all_fnames, nc=_all_nc),
        nmf_dir + "/aucs/aucs.csv",
        nmf_dir + "/aucs/best_med_aucs.csv",
        expand(nmf_scores, nc=_all_nc),
    ]

rule all:
    input:
        _nmf_hist_outputs,
        _nmf_analysis_outputs,

rule raw_histograms:
    input:
        bam = frag_bams + "/{lib}_filt.bam",
    output:
        metrics   = nmf_raw_hist_dir + "/{lib}_metrics.txt",
        histogram = nmf_raw_hist_dir + "/{lib}_histogram.pdf",
    conda:
        ENV_FRAGENTOMICS,
    log:
        logdir + "/nmf/raw_histograms/{lib}.log",
    params:
        out_dir = nmf_raw_hist_dir,
        log_dir = logdir + "/nmf/raw_histograms",
        width   = config["histogram_width"],
        picardpath = config["picardpath"],
    shell:
        """
        mkdir -p {params.out_dir}
        mkdir -p {params.log_dir}

        #module load picard
        #module load R

        picard CollectInsertSizeMetrics \
        -I {input.bam} \
        -O {output.metrics} \
        -H {output.histogram} \
        -W {params.width} \
        &> {log}
        """
rule agg_histograms:
    input:
        metrics = expand(nmf_raw_hist_dir + "/{lib}_metrics.txt", lib=FRAG_LIBS),
        plots = expand(nmf_raw_hist_dir + "/{lib}_histogram.pdf", lib=FRAG_LIBS)
    output:
        counts = nmf_agg_hist_dir + "/count_histogram.csv",
        freqs = nmf_agg_hist_dir + "/freq_histogram.30_700.csv",
    conda:
        ENV_FRAGENTOMICS,
    log:
        logdir + "/nmf/agg_histograms.log",
    params:
        hist_dir = nmf_raw_hist_dir,
        out_dir = nmf_agg_hist_dir,
        logs_dir = logdir + "/nmf",
        script = frag_scriptdir + "/agg_histograms.py"
    shell:
        """
        mkdir -p {params.out_dir}
        mkdir -p {params.logs_dir}
        #module load python/3.9
        python3 {params.script} \
        --hist_dir {params.hist_dir} \
        --output {params.out_dir} \
        --verbose \
        &> {log}
        """

# Modes 1/2: compute AUC across all component numbers + component plots
if not APPLY_NMF_CUTPOINTS:
    rule nmf:
        input:
            freq_csv = nmf_agg_hist_dir + "/freq_histogram.30_700.csv",
            counts = nmf_agg_hist_dir + "/count_histogram.csv",

            raw_metrics = expand(nmf_raw_hist_dir + "/{lib}_metrics.txt", lib=FRAG_LIBS),
            raw_plots = expand(nmf_raw_hist_dir + "/{lib}_histogram.pdf", lib=FRAG_LIBS)
        output:
            probs = nmf_dir + "/probs/{fnames}.n={nc}.csv",
            coeffs = nmf_dir + "/coeffs/{fnames}.n={nc}.csv",
            model = nmf_dir + "/model/{fnames}.n={nc}.pkl",
        conda:
            ENV_FRAGENTOMICS,
        log:
            logdir + "/nmf/model/{fnames}.n={nc}.log",
        params:
            out_dir = nmf_dir,
            log_dir = logdir + "/nmf/model",
            script = frag_scriptdir + "/nmf_model.py",
            samplelist = samples,
            existing_model = config["existing_nmf_model"]
        shell:
            """
            mkdir -p {params.out_dir}/{{model,probs,coeffs}} {params.log_dir}
            #module load python/3.9
            python3 {params.script} \
                --freq_csv {input.freq_csv} \
                --samplelist {params.samplelist} \
                --n_components {wildcards.nc} \
                --cohort1 $(echo "{wildcards.fnames}" | sed 's/_v_.*//') \
                --cohort2 $(echo "{wildcards.fnames}" | sed 's/.*_v_//') \
                --out_dir {params.out_dir} \
                --existing_model {params.existing_model} \
                --name {wildcards.fnames}.n={wildcards.nc} \
                &> {log}
            """
    rule nmf_aucs:
        input:
            expand(nmf_probs,
                    nc = range(min_components,max_components+1)
            ),
            samplelist = samples,
        output:
            nmf_dir + "/aucs/aucs.csv",
            nmf_dir + "/aucs/best_med_aucs.csv",
            expand(nmf_scores,
                    nc = range(min_components,max_components+1)
            ),
        conda:
            ENV_FRAGENTOMICS,
        log:
            logdir + "/nmf/auc/aucs.log"
        params:
            out_dir = nmf_dir,
            log_dir = logdir + "/nmf/auc",
            script  = frag_scriptdir + "/nmf_aucs.py",
            probs_dir = nmf_dir + "/probs",
            min_components = min_components,
            max_components = max_components,
            cohorts = ",".join(config["cohorts"]),
        shell:
            """
            mkdir -p {params.out_dir}
            mkdir -p {params.log_dir}

            python3 {params.script} \
            --probs_dir {params.probs_dir} \
            --min_components {params.min_components} \
            --max_components {params.max_components} \
            --samplelist {input.samplelist} \
            --out_dir {params.out_dir} \
            --cohorts {params.cohorts} \
            &> {log}
            """

    rule plot_nmf_components:
        input:
            model = nmf_dir + "/model/{fnames}.n={nc}.pkl",
            freq_csv = nmf_agg_hist_dir + "/freq_histogram.30_700.csv",
        output:
            plot = nmf_dir + "/plots/components/{fnames}.n={nc}.pdf",
        conda:
            ENV_FRAGENTOMICS,
        log:
            logdir + "/nmf/plots/components/{fnames}.n={nc}.log"
        params:
            out_dir = nmf_dir + "/plots/components",
            log_dir = logdir + "/nmf/plots/components",
            script = frag_scriptdir + "/plot_nmf_components.py",
        shell:
            """
            mkdir -p {params.out_dir}
            mkdir -p {params.log_dir}

            python3 {params.script} \
                --model {input.model} \
                --freq_csv {input.freq_csv} \
                --output {output.plot} \
                --cohort1 $(echo "{wildcards.fnames}" | sed 's/_v_.*//') \
                --cohort2 $(echo "{wildcards.fnames}" | sed 's/.*_v_//') \
                --n_components {wildcards.nc} \
                &> {log}
            """

# Mode 3: apply existing cutpoints to best n_components probs, no new AUC
else:
    rule nmf:
        input:
            freq_csv = nmf_agg_hist_dir + "/freq_histogram.30_700.csv",
            counts = nmf_agg_hist_dir + "/count_histogram.csv",

            raw_metrics = expand(nmf_raw_hist_dir + "/{lib}_metrics.txt", lib=FRAG_LIBS),
            raw_plots = expand(nmf_raw_hist_dir + "/{lib}_histogram.pdf", lib=FRAG_LIBS)
        output:
            probs = nmf_dir + "/probs/{fnames}.n={nc}.csv",
            #coeffs = nmf_dir + "/coeffs/{fnames}.n={nc}.csv",
            model = nmf_dir + "/model/{fnames}.n={nc}.pkl",
        conda:
            ENV_FRAGENTOMICS,
        log:
            logdir + "/nmf/model/{fnames}.n={nc}.log",
        params:
            out_dir = nmf_dir,
            log_dir = logdir + "/nmf/model",
            script = frag_scriptdir + "/nmf_model_mod3.py",
            samplelist = samples,
            existing_model = config["existing_nmf_model"]
        shell:
            """
            mkdir -p {params.out_dir}/{{model,probs,coeffs}} {params.log_dir}
            #module load python/3.9
            python3 {params.script} \
                --freq_csv {input.freq_csv} \
                --samplelist {params.samplelist} \
                --n_components {wildcards.nc} \
                --cohort1 $(echo "{wildcards.fnames}" | sed 's/_v_.*//') \
                --cohort2 $(echo "{wildcards.fnames}" | sed 's/.*_v_//') \
                --out_dir {params.out_dir} \
                --existing_model {params.existing_model} \
                --name {wildcards.fnames}.n={wildcards.nc} \
                &> {log}
            """
    rule nmf_apply_cutpoints:
        input:
            probs = nmf_dir + "/probs/{fnames}.n={nc}.csv",
        output:
            scores = nmf_dir + "/scores/{fnames}.n={nc}.csv",
        conda:
            ENV_FRAGENTOMICS,
        log:
            logdir + "/nmf/apply_cutpoints/{fnames}.n={nc}.log"
        params:
            existing_nmf_dir = config.get("existing_nmf_scores_dir", "no"),
            script = frag_scriptdir + "/apply_cutpoints.py",
            samplelist = samples,
            log_dir = logdir + "/nmf/apply_cutpoints",
        shell:
            """
            mkdir -p {params.log_dir}
            mkdir -p $(dirname {output.scores})

            python3 {params.script} \
            --probs_csv {input.probs} \
            --scores_csv {params.existing_nmf_dir}/scores/{wildcards.fnames}.n={wildcards.nc}.csv \
            --samplelist {params.samplelist} \
            --cohort1 $(echo "{wildcards.fnames}" | sed 's/_v_.*//') \
            --cohort2 $(echo "{wildcards.fnames}" | sed 's/.*_v_//') \
            --out_csv {output.scores} \
            &> {log}
            """
