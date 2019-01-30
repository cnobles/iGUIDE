#!/usr/bin/env bash
set -e

# Input arguments
__IGUIDE_ENV=${1-iguide}
__CORES=${2-1}

# Clear test directory
rm -rf analysis/simulation

# Test script
source ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate ${__IGUIDE_ENV}

# Create test analysis directory
iguide setup configs/simulation.config.yml -- -np
iguide setup configs/simulation.config.yml -- --nolock

# Generate test DAG graph
iguide run configs/simulation.config.yml -- -np

#iguide run configs/simulation.config.yml -- --dag --nolock | dot -Tsvg > \
#    analysis/simulation/reports/simulation.dag.svg

iguide run configs/simulation.config.yml -- -w 30 --nolock --cores ${__CORES}

iguide report analysis/simulation/output/edited_sites.simulation.rds \
    -c configs/simulation.config.yml \
    -o analysis/simulation/reports/report.simulation.pdf \
    -s sampleInfo/simulation.supp.csv \
    -t pdf
