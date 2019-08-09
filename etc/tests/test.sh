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
iguide setup configs/simulation.config.yml

# Generate test DAG graph and run
iguide run configs/simulation.config.yml -- -np

iguide run configs/simulation.config.yml -- --dag --nolock | dot -Tsvg > \
    analysis/simulation/reports/simulation.dag.svg

iguide run configs/simulation.config.yml -- -p -w 30 --notemp --nolock --cores ${__CORES}

# Evaluate and report out using a different metadata set
iguide eval configs/simulation.config.yml \
    -o analysis/simulation/reports/iguide.eval.simulation.test.rds \
    -s sampleInfo/simulation.supp.csv 

iguide report -e analysis/simulation/reports/iguide.eval.simulation.test.rds \
    -o analysis/simulation/reports/report.simulation.test \
    -g

# Test for accuracy and retention
Rscript tools/rscripts/check_test_accuracy.R configs/simulation.config.yml \
    etc/tests/Data/truth.csv -v

# Test for precise outputs if using config::Aligner : "blat"
function __blat_aligner () {
    grep Aligner configs/simulation.config.yml | grep blat > /dev/null && echo true || echo false
}

if [[ $(__blat_aligner) = true ]]; then
    Rscript tools/rscripts/check_file_digests.R etc/tests/simulation.digests.yml -v
fi

# Cleanup
iguide clean configs/simulation.config.yml
iguide clean configs/simulation.config.yml --remove_proj

# Deactivate conda environment
conda deactivate
