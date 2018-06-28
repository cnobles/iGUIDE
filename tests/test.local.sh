#!/usr/bin/env bash
set -ev

# Test script
PREFIX=${HOME}/miniconda3
export PATH=${PATH}:${PREFIX}/bin
source activate iguide
CORES=${1-1}

# Create test analysis directory
snakemake analysis/simulation --configfile configs/simulation.config.yml -np
snakemake analysis/simulation --configfile configs/simulation.config.yml

# Move test sequence files to analysis directory
cp tests/Data/Undetermined_S0_L001_* analysis/simulation/input_data/

# Generate test DAG graph
snakemake --configfile configs/simulation.config.yml -np
snakemake --configfile configs/simulation.config.yml --dag | dot -Tsvg > analysis/simulation/reports/simulation.dag.svg
snakemake --configfile configs/simulation.config.yml --latency-wait 30 --notemp --cores ${CORES}
head analysis/simulation/output/unique_sites.simulation.csv
