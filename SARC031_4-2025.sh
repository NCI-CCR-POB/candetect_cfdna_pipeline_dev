#! /bin/bash

## NOTE: Copy this template and edit per job. Be sure to make the changes below:
#    -Must change path to repo directory. 
#    -Must change name of .yaml in all lines. 
#    -Must change start script path to new .sh file name that this template is copied to generate.

module load snakemake

cd /data/sundbyrt/repos/candetect_cfdna_pipeline #Update path before running

#module load mamba_install
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)" # Necessary to run conda snakemake command in shell script

#### Dry-Run ####
#Un-comment line-by-line to test each stage
#snakemake -s workflows/frag_reads.smk --configfile config/SARC031.yaml --profile profile/ -np
#snakemake -s workflows/end_motifs.smk --configfile config/SARC031.yaml --profile profile/ -np
#snakemake -s workflows/binwise_frags.smk --nolock  --configfile config/SARC031.yaml --profile profile/ -np
#snakemake -s workflows/NMF_size.smk --configfile config/SARC031.yaml --profile profile/ -np
#snakemake -s workflows/ichor.smk --configfile config/SARC031.yaml --profile profile/ -np
#snakemake -s workflows/combine_datasets.smk --nolock --configfile config/SARC031.yaml --profile profile/ -np

#### Actual Run #### 
#As this repo is currently set up, each line must be uncommented individually and run one-at-a-time (note: for end_motifs, binwise_frags, ichor, NMF_size you could add --nolock and run simultaneously)

#snakemake -s workflows/frag_reads.smk --configfile config/SARC031.yaml --profile profile/ 
#snakemake -s workflows/end_motifs.smk --nolock --configfile config/SARC031.yaml --profile profile/ 
snakemake -s workflows/binwise_frags.smk --nolock --configfile config/SARC031.yaml --profile profile/ 
#snakemake -s workflows/NMF_size.smk --nolock --configfile config/SARC031.yaml --profile profile/
#snakemake -s workflows/ichor.smk --nolock --configfile config/SARC031.yaml --profile profile/
#snakemake -s workflows/combine_datasets.smk --nolock --configfile config/SARC031.yaml --profile profile/ 

## START SCRIPT: Copy and paste this line, updated with paths, to terminal on the HPC to launch the smk pipeline
#sbatch --mem=1g --cpus-per-task=2 --time=120:00:00 --mail-type=BEGIN,TIME_LIMIT_90,END /data/sundbyrt/repos/candetect_cfdna_pipeline/SARC031_4-2025.sh