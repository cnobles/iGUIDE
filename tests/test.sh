#!/usr/bin/env bash
set -ev

# Test script
PREFIX=${HOME}/miniconda3
export PATH=${PATH}:${PREFIX}/bin
source activate iguide
CORES=${1-1}

# Create test analysis directory
iguide setup configs/simulation.config.yml -- -np
iguide setup configs/simulation.config.yml

# Generate test DAG graph
iguide run configs/simulation.config.yml -- -np
iguide run configs/simulation.config.yml -- --dag | dot -Tsvg > analysis/simulation/reports/simulation.dag.svg
iguide run configs/simulation.config.yml -- --latency-wait 30 --cores ${CORES}
head analysis/simulation/output/unique_sites.simulation.csv
