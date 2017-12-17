#!/bin/bash
set -e

# Test script
PREFIX=${HOME}/miniconda3
export PATH=${PATH}:${PREFIX}/bin
source activate idsb

# Create test analysis directory
snakemake analysis/test --configfile configs/test.config.yml -np
snakemake analysis/test --configfile configs/test.config.yml

# Move test sequence files to analysis directory
cp tests/Data/Undetermined_S0_L001_* analysis/test/input_data/

# Generate test DAG graph
snakemake --configfile configs/test.config.yml -np
snakemake --configfile configs/test.config.yml --dag | dot -Tsvg > test.dag.svg
snakemake --configfile configs/test.config.yml --latency-wait 30
cat analysis/test/process/processedData/unique_sites.test.csv
