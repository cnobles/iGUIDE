#!/usr/bin/env bash
set -e

# Input arguments
__IGUIDE_ENV=${1-iguide}
__CORES=${2-1}

# Clear test directory
rm -rf analysis/simulation

# Activate conda environment
source ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate ${__IGUIDE_ENV}

# Create test analysis directory
iguide setup configs/simulation.config.yml -- -np
iguide setup configs/simulation.config.yml -- --nolock

# Generate test DAG graph
iguide run configs/simulation.config.yml -- -np

#iguide run configs/simulation.config.yml -- --dag --nolock | dot -Tsvg > \
#    analysis/simulation/reports/simulation.dag.svg

iguide run configs/simulation.config.yml -- -p -w 30 --nolock --cores ${__CORES}

iguide report configs/simulation.config.yml \
    -o analysis/simulation/reports/report.simulation.test.html \
    -s sampleInfo/simulation.supp.csv \
    -t html

# Test for precise outputs
Rscript tools/rscripts/check_file_digests.R tests/simulation.digests.yml -v

# Cleanup
iguide clean configs/simulation.config.yml
iguide clean configs/simulation.config.yml --remove_proj

# Deactivate conda environment
conda deactivate
