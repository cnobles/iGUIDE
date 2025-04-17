#!/usr/bin/env bash
set -e

# Input arguments
__IGUIDE_ENV=${1-iguide}
__CORES=${2-1}
__INSTALL=${3-miniconda3}  # Options include: miniconda3 (default) or anaconda

# Clear test directory
rm -rf analysis/simulation

# Activate conda environment
if [[ ${__INSTALL} = "miniconda3" ]]; then
    source ${HOME}/miniconda3/etc/profile.d/conda.sh
elif [[ ${__INSTALL} = "anaconda" ]]; then
    source ${HOME}/anaconda/etc/profile.d/conda.sh
else
    echo "Please use either 'miniconda3' or 'anaconda' for managing conda environments."
    exit 1
fi

conda activate ${__IGUIDE_ENV}

# Create test analysis directory
iguide setup configs/simulationA.config.yml
iguide setup configs/simulationB.config.yml

# Generate test DAG graph and run
iguide run configs/simulationA.config.yml -- -np

iguide run configs/simulationA.config.yml -- --dag --nolock | dot -Tsvg > \
    analysis/simulationA/reports/simulationA.dag.svg

iguide run configs/simulationA.config.yml -- -p -w 30 --notemp --nolock --cores ${__CORES}

# Evaluate and report out using a different metadata set
iguide eval configs/simulationA.config.yml \
    -o analysis/simulationA/reports/iguide.eval.simulationA.test.rds \
    -s sampleInfo/simulationA.supp.csv 

iguide report -e analysis/simulationA/reports/iguide.eval.simulationA.test.rds \
    -o analysis/simulationA/reports/report.simulationA.test \
    -g

# Generate simulation B data
iguide run configs/simulationB.config.yml -- -p -w 30 --notemp --nolock --cores ${__CORES}

# Generate combination simulation data
iguide eval configs/simulationA.config.yml configs/simulationB.config.yml \
    -o analysis/simulationB/reports/iguide.eval.combination.test.rds \
    -s sampleInfo/simulationB.supp.csv \
    --stat analysis/simulationB/reports/iguide.stat.combination.test.csv

iguide report -e analysis/simulationB/reports/iguide.eval.combination.test.rds \
    -o analysis/simulationB/reports/report.combination.test \
    -g

iguide summary -e analysis/simulationB/reports/iguide.eval.combination.test.rds \
    -o analysis/simulationB/reports/summary.combination.test
    
# Test for accuracy and retention
Rscript tools/rscripts/check_test_accuracy.R configs/simulationA.config.yml \
    etc/tests/DataA/truth.csv -v

Rscript tools/rscripts/check_test_accuracy.R configs/simulationB.config.yml \
    etc/tests/DataB/truth.csv -v

# Test for precise outputs if using config::Aligner : "blat"
function __blat_aligner () {
    grep Aligner configs/simulationA.config.yml | grep blat > /dev/null && echo true || echo false
}

if [[ $(__blat_aligner) = true ]]; then
    Rscript tools/rscripts/check_file_digests.R etc/tests/simulation.digests.yml -v
fi

# Cleanup
#iguide clean configs/simulationA.config.yml
#iguide clean configs/simulationA.config.yml --remove_proj
#iguide clean configs/simulationB.config.yml
#iguide clean configs/simulationB.config.yml --remove_proj

# Deactivate conda environment
conda deactivate
