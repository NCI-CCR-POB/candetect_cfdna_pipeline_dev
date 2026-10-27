#! /bin/bash

## NOTE: Copy this template and edit per job. Be sure to make the changes below:
#    -Must change SMK_DIR path to repo directory. 
#    -Must change CONFIG_FILE to the name of the .yaml
#    -Update sbatch line before copying and pasting this into command line

# Set paths that will be parsed for  {params}
SMK_DIR="/path/to/repos/candetect_cfdna_pipeline"
CONFIG_FILE="config/default.yaml"

module load snakemake
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)" # Necessary to run conda snakemake command in shell script
cd $SMK_DIR

#### Run the Snakemake master workflow with the specified configuration ####
#snakemake -s workflows/master_workflow.smk --configfile $CONFIG_FILE --profile profile/ -np #Dry run
snakemake -s workflows/master_workflow.smk --configfile $CONFIG_FILE --profile profile/   #Actual run

## START SCRIPT: Copy and paste the below line (without the '#'), updated with paths, to terminal on the HPC to launch the smk pipeline
# sbatch --mem=1g --cpus-per-task=2 --time=120:00:00 --mail-type=BEGIN,TIME_LIMIT_90,END /path/to/repos/candetect_cfdna_pipeline/Default_launch_master_workflow.sh 