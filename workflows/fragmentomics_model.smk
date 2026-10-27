#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#     Fragmentomic Analysis of Cell-free DNA Whole Genome Sequencing           #
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
ENV_FRAG_READS      = config["env_frag_reads"]



# Directory values derived from data_dir in configuration YAML
benchdir            = data_dir + "/benchmark"
frag                = data_dir + "/analysis/cfdna_wgs"
frag_bams           = data_dir + "/bams"
frag_fastqs         = data_dir + "/fastq"
frag_scriptdir      = frag_repo +  "/scripts"
logdir_binwise            = data_dir + "/logs/bin_wise_fragmentomics"
logdir_nmf           = data_dir + "/logs/nmf_size_fragmentomics"
logdir_motif            = data_dir + "/logs/end_motif_fragmentomics"
qcdir               = data_dir + "/qc"

bench_dir           = data_dir + "/benchmark"
frag_analysis_dir   = data_dir + "/analysis/cfdna_wgs"
logdir              = data_dir + "/logs/preprocessing"
qcdir               = data_dir  + "/qc"
tmpdir              = data_dir + "/tmp"
print(frag_scriptdir)

log_dir_combine = data_dir + "/logs/combine_datasets"
plot_dir_combine = data_dir + "/analysis/cfdna_wgs/plots"


#Make sure logdir exists
#Make sure logdir exists
os.makedirs(logdir, exist_ok=True)
os.makedirs(logdir_binwise, exist_ok=True)
os.makedirs(logdir_nmf, exist_ok=True)
os.makedirs(logdir_motif, exist_ok=True)
# Create directories if they don't exist
for d in [log_dir_combine, plot_dir_combine]:
    os.makedirs(d, exist_ok=True)


os.makedirs(tmpdir, exist_ok=True)

size_selected_dir = config["data_dir"] + "/size_selected_bams"
wig_dir = config["data_dir"] + "/wig_files"
ichor_dir = config["data_dir"] + "/analysis/cfdna_wgs/ichorCNA"
logdir_ichor = config["data_dir"] + "/logs/ichor"

os.makedirs(size_selected_dir, exist_ok=True)
os.makedirs(wig_dir, exist_ok=True)
os.makedirs(ichor_dir, exist_ok=True)
os.makedirs(logdir_ichor, exist_ok=True)



# Bin-wise fragmentomic specific directories
bin_wise_frag            = frag + "/bin_wise_fragmentomics"

bin_wise_beds      = bin_wise_frag + "/beds"
bin_wise_counts     = bin_wise_frag + "/counts"
bin_wise_gc_distros = bin_wise_frag + "/distros"

healthy_med_rds = config.get("existing_health_gc_rds", "")
print("########existing healthy med rds#####################")
print(healthy_med_rds)
print("#####################################################")

# NMF specific inputs
min_components      = config["min_nmf_components"]
max_components      = config["max_nmf_components"]
nc                  = range(config["min_nmf_components"],config["max_nmf_components"]+1)

# NMF specific directories
nmf_dir             =  frag + "/nmf"

nmf_agg_hist_dir    = nmf_dir + "/agg"
nmf_comps           = [f"{nmf_dir}/model/{c1}_v_{c2}.n={{nc}}.pkl" for (c1,c2) in it.combinations(COHORTS,2)]
nmf_probs           = [f"{nmf_dir}/probs/{c1}_v_{c2}.n={{nc}}.csv" for (c1,c2) in it.combinations(COHORTS,2)]
nmf_raw_hist_dir    = nmf_dir + "/raw"
nmf_scores          = [f"{nmf_dir}/scores/{c1}_v_{c2}.n={{nc}}.csv" for (c1,c2) in it.combinations(COHORTS,2)]



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
FRAG_HEALTHY_LIBRARIES = frag_libs[frag_libs['cohort2'] == normal_control]['library'].tolist()

#Define cohort comparisons (one-versus-one)
COHORT_1 = []
COHORT_2 = []
for c1, c2 in it.combinations(COHORTS,2):
    COHORT_1.append(c1)
    COHORT_2.append(c2)
    print(c1,c2)
pairs = []
pairs = list(it.combinations(COHORTS, 2))
print(pairs)
print([pair[0] for pair in pairs])




# binwise output
_binwise_base_outputs = [
    expand(bin_wise_beds + "/{library}_filt.bed", library=frag_libs["library"]),
    expand(bin_wise_gc_distros + "/{library}_gc_distro.csv", library=frag_libs["library"]),
    bin_wise_gc_distros + "/healthy_med.rds",
    expand(bin_wise_beds + "/{library}_sampled_frag.bed", library=frag_libs["library"]),
    expand(bin_wise_beds + "/{library}_norm_short.bed", library=frag_libs["library"]),
    expand(bin_wise_beds + "/{library}_norm_long.bed", library=frag_libs["library"]),
    expand(bin_wise_counts + "/{library}_cnt_short.tmp", library=frag_libs["library"]),
    expand(bin_wise_counts + "/{library}_cnt_long.tmp", library=frag_libs["library"]),
    bin_wise_frag + "/frag_counts.tsv",
    bin_wise_frag + "/ratios.tsv",
    bin_wise_frag + "/binwise_armz.raw.csv",
    bin_wise_frag + "/binwise_armz.clean.csv",
    bin_wise_frag + "/binwise_ratio.norm_rows.clean.csv"
]
_binwise_analysis_outputs = [
    expand(bin_wise_frag + "/model/coeffs/{c1}_v_{c2}.csv", zip, c1=COHORT_1, c2=COHORT_2),
    bin_wise_frag + "/aucs/aucs/aucs.csv",
    expand(bin_wise_frag + "/scores/{c1}_v_{c2}.csv", zip, c1=COHORT_1, c2=COHORT_2),
    expand(bin_wise_frag + "/plots/prob_per_repeat/{c1}_v_{c2}.pdf", zip, c1=COHORT_1, c2=COHORT_2),
    expand(bin_wise_frag + "/plots/prob_per_sample/{c1}_v_{c2}.pdf", zip, c1=COHORT_1, c2=COHORT_2),
    bin_wise_frag + "/binwise/plots/coeff_heatmap/binwise_armz.clean.pdf",
]
#print (_binwise_analysis_outputs)

# outputs for nmf feature
_nmf_hist_outputs = [
    expand(nmf_raw_hist_dir + "/{lib}_metrics.txt", lib=FRAG_LIBS),
    expand(nmf_raw_hist_dir + "/{lib}_histogram.pdf", lib=FRAG_LIBS),
    nmf_agg_hist_dir + "/count_histogram.csv",
    nmf_agg_hist_dir + "/freq_histogram.30_700.csv",
]
##NMF specific
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

# end_motif_outputs (all modes)
_motif_base_outputs = [
    expand(motif_dir + "/{library}_motifs.tsv", library=frag_libs["library"]),
    data_dir + "/analysis/cfdna_wgs/features/combined_motif_counts_matrix.tsv",
    data_dir + "/analysis/cfdna_wgs/features/combined_motif_fractions_matrix.tsv",
    data_dir + "/analysis/cfdna_wgs/features/mds.RData",
    data_dir + "/analysis/cfdna_wgs/features/motif_diversity_scores.csv"
]
_motif_analysis_outputs = [
    data_dir + "/analysis/cfdna_wgs/features/fprofiles.tsv",
    data_dir + "/analysis/cfdna_wgs/features/motif_per_fprofile.tsv",
    data_dir + "/analysis/cfdna_wgs/features/fprofiles.pdf",
    expand(data_dir + "/analysis/cfdna_wgs/features/fprofile_AUCROC/fprofile_AUC_{c1}_v_{c2}.csv",zip, c1=COHORT_1, c2=COHORT_2),
    expand(data_dir + "/analysis/cfdna_wgs/features/fprofile_AUCROC/fprofile_ROC_{c1}_v_{c2}.pdf",zip, c1=COHORT_1, c2=COHORT_2)
]

#ichor output
_ichor_output = [
    expand(size_selected_dir + "/{samples}.ss.bam", samples=frag_libs["library"]),
    expand(wig_dir + "/{samples}.wig", samples=frag_libs["library"]),
    expand(ichor_dir + "/{samples}/{samples}.cna.seg", samples=frag_libs["library"]),
    ichor_dir + "/tf.csv"
]

_combine_output = [
    expand(data_dir + "/analysis/cfdna_wgs/combined_dataset_{c1}_v_{c2}.csv", zip,c1 = COHORT_1,c2 = COHORT_2),
    expand(plot_dir_combine + "/subject_trajectories_{c1}_v_{c2}.png",zip,c1 = COHORT_1,c2 = COHORT_2),
    expand(plot_dir_combine + "/auc_comparison_{c1}_v_{c2}.png",zip,c1 = COHORT_1,c2 = COHORT_2),
    expand(plot_dir_combine + "/roc_comparison_{c1}_v_{c2}.png",zip,c1 = COHORT_1,c2 = COHORT_2),
    expand(plot_dir_combine + "/distributions_{c1}_v_{c2}",zip,c1 = COHORT_1,c2 = COHORT_2),
    expand(data_dir + "/analysis/cfdna_wgs/plots/voting/roc_calls_heatmap_{c1}_v_{c2}.pdf",zip,c1 = COHORT_1,c2 = COHORT_2),
    expand(data_dir + "/analysis/cfdna_wgs/plots/voting/roc_calls_heatmap_{c1}_v_{c2}_majority_vote.csv",zip,c1 = COHORT_1,c2 = COHORT_2),
    data_dir + "/analysis/cfdna_wgs/report/cfdna_wgs_report.html"
]
rule all:
  input:
    _binwise_base_outputs,
    _nmf_hist_outputs,
    _motif_base_outputs,
    _binwise_analysis_outputs,
    _nmf_analysis_outputs,
    _motif_analysis_outputs,
    _ichor_output,
    _combine_output

include: f"{frag_repo}/workflows/subworkflow/binwise_frags_feature.smk"
include: f"{frag_repo}/workflows/subworkflow/NMF_size_feature.smk"
include: f"{frag_repo}/workflows/subworkflow/end_motifs_feature.smk"
include: f"{frag_repo}/workflows/subworkflow/binwise_frags_model.smk"
include: f"{frag_repo}/workflows/subworkflow/NMF_size_model.smk"
include: f"{frag_repo}/workflows/subworkflow/end_motifs_model.smk"
include: f"{frag_repo}/workflows/ichor.smk"
include: f"{frag_repo}/workflows/combine_datasets.smk"
