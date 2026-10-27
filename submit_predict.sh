#!/bin/bash
#SBATCH --partition=norm
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=48g
#SBATCH --time=120:00:00
#SBATCH --no-requeue

mkdir -p logs
mkdir -p stats

snakemake -s workflows/fragmentomics_pred.smk \
    --configfile config/CHLA_set.yaml \
    --cores 4 \
    --latency-wait 600 \
    --use-conda --conda-frontend mamba \
    -k --stats stats/snakemake_$(date +"%Y%m%d%H%M%S").stats \
    --rerun-incomplete --restart-times 0 \
    -j 300 --printshellcmds \
    --cluster-config config/cluster.yaml \
    --cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error}" \
    >& logs/snakemake_$(date +"%Y%m%d%H%M%S").log
