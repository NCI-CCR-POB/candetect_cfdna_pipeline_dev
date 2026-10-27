#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#     Fragmentomic Analysis of Cell-free DNA Whole Genome Sequencing           #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8

import os
import itertools as it
import pandas as pd
import re
import numpy as np

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

# Directory values derived from data_dir in configuration YAML
benchdir            = data_dir + "/benchmark"
frag                = data_dir + "/analysis/cfdna_wgs"
frag_bams           = data_dir + "/bams"
frag_fastqs         = data_dir  + "/fastq"
frag_scriptdir      = frag_repo +  "/scripts"
logdir              = data_dir + "/logs/end_motifs"
qcdir               = data_dir  + "/qc"

# End motif specific directories
motif_dir            = frag + "/motifs"

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
frag_libs_file_indict = frag_libs["file"].tolist()
frag_lib_dict = dict(zip(FRAG_LIBS, frag_libs_file_indict))

# Make  a list of healthy libraries
FRAG_HEALTHY_LIBRARIES = frag_libs[frag_libs['cohort'] == normal_control]['library'].tolist()

# Mode 3: apply existing model + use original fprofile_AUC.csv (no new ROC/AUC)
APPLY_MOTIF_CUTPOINTS = config.get("existing_motif_scores_dir", "no") != "no"

# Shared outputs (all modes)
_motif_base_outputs = [
    expand(motif_dir + "/{library}_motifs.tsv", library=frag_libs["library"]),
    data_dir + "/analysis/cfdna_wgs/features/combined_motif_counts_matrix.tsv",
    data_dir + "/analysis/cfdna_wgs/features/combined_motif_fractions_matrix.tsv",
    data_dir + "/analysis/cfdna_wgs/features/mds.RData",
    data_dir + "/analysis/cfdna_wgs/features/motif_diversity_scores.csv",
    data_dir + "/analysis/cfdna_wgs/features/fprofiles.tsv",
    data_dir + "/analysis/cfdna_wgs/features/motif_per_fprofile.tsv",
    data_dir + "/analysis/cfdna_wgs/features/fprofiles.pdf",
    data_dir + "/analysis/cfdna_wgs/features/fprofile_AUC.csv",
]
# Mode 3: fprofile_AUC.csv copied from existing dir; no ROC plot
if APPLY_MOTIF_CUTPOINTS:
    _motif_analysis_outputs = []
# Modes 1/2: compute full ROC + AUC
else:
    _motif_analysis_outputs = [
        data_dir + "/analysis/cfdna_wgs/features/fprofile_ROC.pdf",
    ]

rule all:
  input:
    _motif_base_outputs,
    _motif_analysis_outputs,

rule bam_to_motif:
    input: frag_bams + "/{library}_filt.bam",
    log: logdir + "/motifs/bam_to_motif/{library}_bam_to_motif.log",
    output: motif_dir + "/{library}_motifs.tsv",
    conda:
        ENV_FRAGENTOMICS,
    params:
        fasta = genome_fasta,
        script = frag_scriptdir + "/extract_5p_motifs.sh",
        threads = getthreads("bam_to_motif"),
        out_dir = motif_dir,
        log_dir = logdir + "/motifs/bam_to_motif",
    shell:
        """
        mkdir -p {params.log_dir}
        mkdir -p {params.out_dir} 
        #module load python/3.9 samtools bedtools
        unset PYTHONPATH
        bash {params.script} \
        {input} \
        {params.fasta} \
        {output} \
        {params.threads} \
        > {log} 2>&1
        """
# Merge short and long fragment counts by genomic poistion for all libraries
rule motif_merge:
    input: expand(motif_dir + "/{library}_motifs.tsv",  library = FRAG_LIBS),
    log: logdir + "/motifs/motif_merge/motif_merge.log",
    output:  
        combined_counts = data_dir + "/analysis/cfdna_wgs/features/combined_motif_counts_matrix.tsv",
        combined_fractions = data_dir + "/analysis/cfdna_wgs/features/combined_motif_fractions_matrix.tsv",
    conda:
        ENV_FRAGENTOMICS,
    params:
        motif_tsv_dir = motif_dir,
        script = frag_scriptdir + "/aggregate_motif_counts.py",
        out_dir = data_dir + "/analysis/cfdna_wgs/features",
        log_dir = logdir + "/motifs/motif_merge",
    shell:
        """
        mkdir -p {params.log_dir}
        mkdir -p {params.out_dir} 
        python {params.script} \
        {params.motif_tsv_dir} \
        {output.combined_counts} \
        {output.combined_fractions} \
        &> {log}
        """
rule mds:
    input: data_dir + "/analysis/cfdna_wgs/features/combined_motif_fractions_matrix.tsv",
    log: logdir + "/motifs/mds/motif_diversity_score.log",
    output:  
        mdsRDS = data_dir + "/analysis/cfdna_wgs/features/mds.RData",
        CSVout = data_dir + "/analysis/cfdna_wgs/features/motif_diversity_scores.csv",
    conda:
        ENV_FRAGENTOMICS,
    params:
        annotation_table = config["samples"],
        script = frag_scriptdir + "/motif_diversity_score.R",
        out_dir = data_dir + "/analysis/cfdna_wgs/features",
        log_dir = logdir + "/motifs/mds",
    shell:
        """
        mkdir -p {params.log_dir}
        mkdir -p {params.out_dir} 
        Rscript  {params.script} \
        {input} \
        {params.annotation_table} \
        {output.mdsRDS} \
        {output.CSVout} \
        &> {log}
        """
rule f_profile:
    input: data_dir + "/analysis/cfdna_wgs/features/combined_motif_fractions_matrix.tsv",
    log: logdir + "/motifs/f_profile/f-profile.log",
    output: 
        fprofile_per_lib = data_dir + "/analysis/cfdna_wgs/features/fprofiles.tsv",
        motif_per_fprofile = data_dir + "/analysis/cfdna_wgs/features/motif_per_fprofile.tsv",
        model_dir = directory(data_dir + "/analysis/cfdna_wgs/features/fprofile_model/"),
    conda:
        ENV_FRAGENTOMICS,
    params:
        script = frag_scriptdir + "/fprofiles_pkl.py",
        log_dir = logdir + "/motifs/f_profile",
        existing_model = config["existing_motif_model"],
    shell:
        r"""
            set -euo pipefail
            mkdir -p {params.log_dir} {output.model_dir}

            # If loading an existing model
            if [[ "{params.existing_model}" != "no" ]]; then
                model_file="{params.existing_model}"
                model_out="{output.model_dir}/fprofile_model.pkl"

                if [[ -f "$model_file" ]]; then
                    cp "$model_file" "$model_out"
                    python {params.script} \
                        --data_file {input} \
                        --output_fprof_per_lib {output.fprofile_per_lib} \
                        --output_motif_per_fprof {output.motif_per_fprofile} \
                        --load_model "$model_out" \
                        > {log} 2>&1
                else
                    echo "Specified model file not found: $model_file" >&2
                    exit 1
                fi
            else
                # Train a new model
                python {params.script} \
                    --data_file {input} \
                    --output_fprof_per_lib {output.fprofile_per_lib} \
                    --output_motif_per_fprof {output.motif_per_fprofile} \
                    --save_model {output.model_dir}/fprofile_model.pkl \
                    > {log} 2>&1
            fi
            """


rule f_profile_plot:
    input: data_dir + "/analysis/cfdna_wgs/features/fprofiles.tsv",
    log: logdir + "/f-profile_plot.log",
    output: 
        data_dir + "/analysis/cfdna_wgs/features/fprofiles.pdf",
    conda:
        ENV_FRAGENTOMICS,
    params:
        script = frag_scriptdir + "/fprofile_stackedbarplot.py",
        plot_title = config["fprofile_stacked_title"],
        metadata_file = metadata_file, #Note, if not using e.g. different cycles, set as "none" in config
        colorscheme = "viridis",
    shell:
        """
        #module load python/3.9
        python {params.script} --input_file {input} --metadata_file {params.metadata_file} \
        --output_file {output} --title "{params.plot_title}" --colormap {params.colorscheme} \
        > {log} 2>&1
	    """      
# Modes 1/2: compute ROC and AUC from new samples
if not APPLY_MOTIF_CUTPOINTS:
    rule fprofile_roc_plots:
        input: data_dir + "/analysis/cfdna_wgs/features/fprofiles.tsv",
        log: logdir + "/motifs/roc_plots/roc_f-profiles.log",
        output:
            plot_out = data_dir + "/analysis/cfdna_wgs/features/fprofile_ROC.pdf",
            csv_out  = data_dir + "/analysis/cfdna_wgs/features/fprofile_AUC.csv",
        conda:
            ENV_FRAGENTOMICS,
        params:
            libraries_tsv = config["samples"],
            script = frag_scriptdir + "/motif_ROC_plot.R",
            out_dir = data_dir + "/analysis/cfdna_wgs/features",
            log_dir = logdir + "/motifs/roc_plots",
        shell:
            """
            mkdir -p {params.log_dir}
            mkdir -p {params.out_dir}
            Rscript {params.script} \
            {input} \
            {params.libraries_tsv} \
            {output.plot_out} \
            {output.csv_out} \
            &> {log}
            """

# Mode 3: copy fprofile_AUC.csv from original run; no new ROC plot
else:
    rule motif_use_existing_auc:
        output:
            fprofile_auc = data_dir + "/analysis/cfdna_wgs/features/fprofile_AUC.csv",
        conda:
            ENV_FRAGENTOMICS,
        log:
            logdir + "/motifs/use_existing_auc.log"
        params:
            existing_auc = os.path.join(config.get("existing_motif_scores_dir", ""), "aucs/aucs.csv"),
            out_dir = data_dir + "/analysis/cfdna_wgs/features",
        shell:
            """
            mkdir -p {params.out_dir}
            if [ ! -f "{params.existing_auc}" ]; then
                echo "ERROR: fprofile_AUC.csv not found at {params.existing_auc}" | tee {log}
                exit 1
            fi
            cp {params.existing_auc} {output.fprofile_auc}
            echo "Copied existing fprofile_AUC.csv from {params.existing_auc}" | tee {log}
            """
