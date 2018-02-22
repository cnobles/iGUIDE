#!/usr/bin/env bash
set -ev

# Test script
PREFIX=${HOME}/miniconda3
export PATH=${PATH}:${PREFIX}/bin
source activate iguide

# Create test analysis directory
snakemake analysis/travis --configfile configs/travis.config.yml -np
snakemake analysis/travis --configfile configs/travis.config.yml

# Move test sequence files to analysis directory
cp tests/Data/Undetermined_S0_L001_* analysis/travis/input_data/

# Generate test DAG graph
snakemake --configfile configs/travis.config.yml -np
snakemake --configfile configs/travis.config.yml --dag | dot -Tsvg > travis.dag.svg
snakemake --configfile configs/travis.config.yml --latency-wait 30
cat analysis/travis/output/unique_sites.travis.csv
