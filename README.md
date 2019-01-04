## iGUIDE - improved Genome-wide Unbiased Identification of Double-strand DNA break Events
[![Build Status](https://travis-ci.org/cnobles/iGUIDE.svg?branch=master)](https://travis-ci.org/cnobles/iGUIDE)
[![Snakemake](https://img.shields.io/badge/snakemake-≥3.5.2-brightgreen.svg?style=flat)](http://snakemake.bitbucket.org)
[![Documentation Status](https://readthedocs.org/projects/iguide/badge/?version=latest)](http://iguide.readthedocs.io/en/latest/?badge=latest)

Bioinformatic pipeline for processing iGUIDE and GUIDE-seq samples.

### Description
Software pipeline for processing and analyzing double-strand DNA break events. These events may be induced, such as by designer nucleases like Cas9, or spontaneous, as produced through DNA replication or ionizing radiation. A laboratory bench-side protocol accompanies this software pipeline, and can be found XXX. 

Below, this readme gives the reader a overview of the pipeline, including how to install and process a sample dataset. Processing a sample data set is broken into three parts: 

1) developing a configuration file and sample information
2) setting up a run directory and acquiring the sequence data
3) initializing the pipeline and understanding the output

More complete documentation can be found on [ReadTheDocs.io](https://iguide.readthedocs.io/en/latest/index.html).

### Install
To install iGUIDE, simply clone the repository to the desired destination:

```
git clone https://github.com/cnobles/iGUIDE.git
```

Then initiate the install using the install script. If you would like the installed environment to be named something other than 'iguide', the new conda environment name can be provided to the 'install.sh' script as provided below.

```
cd path/to/iGUIDE
bash install.sh

# Or

cd path/to/iGUIDE
bash install.sh -e {env_name}
```

### An Example Run
To perform a local test of running the iGUIDE informatic pipeline, run the below code after installing. This block first activates your conda environment, 'iguide' by default, and then creates a test directory within the analysis directory. The run information is stored in the run specific configuration file (config file). Using the '-np' flag with the snakemake call will perform a dry-run (won't actually process anything) and print the commands to the terminal, so you can see what snakemake is about to perform. Next, the test data is moved to the input directory underneath the new test run directory. Then the entirety of processing can start. Using the '--dag' flag and piping the output to 'dot -Tsvg' will generate a vector graphic of the directed acyclic graph (dag) workflow that snakemake will follow given the data provided. 

```
# If conda is not in your path ...
PREFIX=${HOME}/miniconda3
export PATH=${PATH}:${PREFIX}/bin

# Activate iguide environment
conda activate iguide

# Check the setup workflow
iguide setup configs/simulation.config.yml -- -np

# Create test analysis directory
iguide setup configs/simulation.config.yml

# Check the run workflow
iguide run configs/simulation.config.yml -- -np

# Generate test DAG graph
iguide run configs/simulation.config.yml -- --dag | dot -Tsvg > \
    analysis/simulation/reports/simulation.dag.svg

# Run iGUIDE to analyze simulation dataset
iguide run configs/simulation.config.yml -- --latency-wait 30

# Output some of the analysis data to check.
cat analysis/simulation/output/unique_sites.simulation.csv

# Deactivate the environment
conda deactivate
```

### Changelog:

**v0.9.0 (January 4th, 2019)**

* Initial release.
* Supports setup and analysis of GUIDE-seq and iGUIDE experiments.
* Documentation on [ReadTheDocs.io](https://iguide.readthedocs.io/en/latest/index.html).
