## iGUIDE - integrative Genome-wide Unbiased Identification of Double-strand DNA break Events
[![Build Status](https://travis-ci.org/cnobles/iGUIDE.svg?branch=master)](https://travis-ci.org/cnobles/iGUIDE)
[![Snakemake](https://img.shields.io/badge/snakemake-≥3.5.2-brightgreen.svg?style=flat)](http://snakemake.bitbucket.org)
Bioinformatic pipeline for processing iGUIDE and original GUIDE-seq samples.

### Description

### Install
To install iGUIDE, simply clone the repository to the desired destination:
```
git clone https://github.com/cnobles/iGUIDE.git
```

Then initiate the install using the install script. If you would like the installed environment to be named something other than 'iguide', the new conda environment name can be provided to the 'install.sh' script as provided below.
```
cd path/to/iGUIDE
bash bin/install.sh

# Or

cd path/to/iGUIDE
bash bin/install.sh {env_name}
```

### Usage

#### Configuration Files - Configs
Configuration files, or configs for short, contain both run-related and pipeline-related information. This is by design, for reproducibility it is easiest to have what was processed and how it was processed in the same location. There should be one config file for each sequencing run to be processed. Below is a brief summary of how to 'configure' your config file to your specific run.

##### Layout


#### General

#### Configuration Files

#### SampleInfo / Manifests

### Tests
To perform a local test of running the iGUIDE informatic pipeline, run the below code after installing. This block first activates your conda environment, 'iguide' by default, and then creates a test directory within the analysis directory. The run information is stored in the run specific configuration file (config file). Using the '-np' flag with the snakemake call will perform a dry-run (won't actually process anything) and print the commands to the terminal, so you can see what snakemake is about to perform. Next, the test data is moved to the input directory underneath the new test run directory. Then the entirety of processing can start. Using the '--dag' flag and piping the output to 'dot -Tsvg' will generate a vector graphic of the directed acyclic graph (dag) workflow that snakemake will follow given the data provided. 

```
# Test script
PREFIX=${HOME}/miniconda3
export PATH=${PATH}:${PREFIX}/bin
source activate iguide

# Create test analysis directory
snakemake analysis/test --configfile configs/test.config.yml -np
snakemake analysis/test --configfile configs/test.config.yml

# Move test sequence files to analysis directory
cp tests/Data/Undetermined_S0_L001_* analysis/test/input_data/

# Generate test DAG graph
snakemake --configfile configs/test.config.yml -np
snakemake --configfile configs/test.config.yml --dag | dot -Tsvg > test.dag.svg
snakemake --configfile configs/test.config.yml --latency-wait 30
cat analysis/test/output/unique_sites.test.csv
```
### Uninstall
To uninstall iGUIDE, the user will need to remove the environment and the directory.

To remove the environment and channels used with conda:
```
cd path/to/iGUIDE
bash bin/uninstall.sh

# Or

cd path/to/iGUIDE
bash bin/uninstall.sh {env_name}
```

If the user would rather remove the environment created for iGUIDE, it is recommended directly use conda. This will leave the channels within the conda config for use with other conda configurations.
```
conda env remove -n iguide

# Or

conda env remove -n {env_name}
```

To remove the iGUIDE directory and conda, the following two commands can be used:
```
# Remove iGUIDE directory and software
rm -r path/to/iGUIDE

# Remove conda
rm -r path/to/miniconda3
```
