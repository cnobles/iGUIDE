#!/usr/bin/env bash
set -e

# Input arguments
__IGUIDE_ENV=${1-iguide}
__CORES=${2-1}
__INSTALL=${3-miniconda3}  # Options include: miniconda3 (default) or anaconda

# Clear test directory
rm -rf analysis/simulation*

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
echo "Test 01 - Setup analysis directory for simulation A."
iguide setup configs/simulationA.config.yml
echo "Test 01 - PASS."

echo "Test 02 - Setup analysis directory for simulation B."
iguide setup configs/simulationB.config.yml
echo "Test 02 - PASS."

# Generate test DAG graph and run
echo "Test 03 - Dry run of simulation A."
iguide run configs/simulationA.config.yml -- -np
echo "Test 03 - PASS."

echo "Test 04 - DAG generation of simulation A workflow."

iguide run configs/simulationA.config.yml -- --dag --nolock | dot -Tsvg > \
    analysis/simulationA/reports/simulationA.dag.svg

if [[ -f "analysis/simulationA/reports/simulationA.dag.svg" ]]; then
    echo "Test 04 - PASS."
fi

echo "Test 05 - Run simulation A workflow."
iguide run configs/simulationA.config.yml -- -p -w 30 --notemp --nolock --cores ${__CORES}
echo "Test 05 - PASS."

# Evaluate and report out using a different metadata set
echo "Test 06 - Run Evaluation workflow on simulation A."

iguide eval configs/simulationA.config.yml \
    -o analysis/simulationA/reports/iguide.eval.simulationA.test.rds \
    -s sampleInfo/simulationA.supp.csv 

if [[ -f "analysis/simulationA/reports/iguide.eval.simulationA.test.rds" ]]; then
    echo "Test 06 - PASS."
fi

echo "Test 07 - Generate report for simulation A."

iguide report -e analysis/simulationA/reports/iguide.eval.simulationA.test.rds \
    -o analysis/simulationA/reports/report.simulationA.test \
    -g

if [[ -f "analysis/simulationA/reports/report.simulationA.test.html" ]]; then
    echo "Test 07 - PASS."
fi

# Generate simulation B data
echo "Test 08 - Run Evaluation workflow on simulation B."
iguide run configs/simulationB.config.yml -- -p -w 30 --notemp --nolock --cores ${__CORES}
echo "Test 08 - PASS."

# Generate combination simulation data
echo "Test 09 - Run evaluation on combination of simulations."

iguide eval configs/simulationA.config.yml configs/simulationB.config.yml \
    -o analysis/simulationB/reports/iguide.eval.combination.test.rds \
    -s sampleInfo/simulationB.supp.csv \
    --stat analysis/simulationB/reports/iguide.stat.combination.test.csv

if [[ -f "analysis/simulationB/reports/iguide.eval.combination.test.rds" ]]; then
    echo "Test 09 - PASS."
fi

echo "Test 10 - Generate report for combination of simulations."

iguide report -e analysis/simulationB/reports/iguide.eval.combination.test.rds \
    -o analysis/simulationB/reports/report.combination.test \
    -g

if [[ -f "analysis/simulationB/reports/report.combination.test.html" ]]; then
    echo "Test 10 - PASS."
fi

echo "Test 11 - Generate summary for combination of simulations."

iguide summary -e analysis/simulationB/reports/iguide.eval.combination.test.rds \
    -o analysis/simulationB/reports/summary.combination.test

if [[ -f "analysis/simulationB/reports/summary.combination.test.txt" ]]; then
    echo "Test 11 - PASS."
fi

# Test for accuracy and retention
echo "Test 12 - Check test accuracy for simulation A."
Rscript tools/rscripts/check_test_accuracy.R configs/simulationA.config.yml \
    etc/tests/DataA/truth.csv -v
echo "Test 12 - PASS."

echo "Test 13 - Check test accuracy for simulation B."
Rscript tools/rscripts/check_test_accuracy.R configs/simulationB.config.yml \
    etc/tests/DataB/truth.csv -v
echo "Test 13 - PASS."

# Test for precise outputs if using config::Aligner : "blat"
function __blat_aligner () {
    grep Aligner configs/simulationA.config.yml | grep blat > /dev/null && echo true || echo false
}

if [[ $(__blat_aligner) = true ]]; then
    echo "Test 14 - Check file digests for simulations if testing with BLAT aligner."
    Rscript tools/rscripts/check_file_digests.R etc/tests/simulation.digests.yml -v
    echo "Test 14 - PASS."
fi

# Cleanup
echo "Test 15 - Clean up analysis directories."
iguide clean configs/simulationA.config.yml
iguide clean configs/simulationA.config.yml --remove_proj
iguide clean configs/simulationB.config.yml
iguide clean configs/simulationB.config.yml --remove_proj
echo "Test 15 - PASS."

# Deactivate conda environment
conda deactivate
