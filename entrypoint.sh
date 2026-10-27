#!/bin/bash

# Activate conda environment
source /opt/conda/etc/profile.d/conda.sh
conda activate LIQUORICE

# Set up environment variables for tools
export PATH="/opt/conda/bin:${PATH}"
export PATH="/opt/conda/envs/LIQUORICE/bin:${PATH}"

# Execute the command passed to the container
exec "$@" 