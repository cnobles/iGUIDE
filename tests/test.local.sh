#!/usr/bin/env bash
set -ev

# Test script
PREFIX=${HOME}/miniconda3
export PATH=${PATH}:${PREFIX}/bin
source activate iguide
CORES=${1-1}

# Create test analysis directory
snakemake analysis/default --configfile configs/default.config.yml -np
snakemake analysis/default --configfile configs/default.config.yml

# Move test sequence files to analysis directory
cp tests/Data/Undetermined_S0_L001_* analysis/default/input_data/

# Generate test DAG graph
snakemake --configfile configs/default.config.yml -np
snakemake --configfile configs/default.config.yml --dag | dot -Tsvg > analysis/default/reports/default.dag.svg
snakemake --configfile configs/default.config.yml --latency-wait 30 --notemp --cores ${CORES}
head analysis/default/output/unique_sites.default.csv
