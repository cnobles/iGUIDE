## iGUIDE - improved Genome-wide Unbiased Identification of Double-strand DNA break Events
[![Build Status](https://travis-ci.org/cnobles/iGUIDE.svg?branch=master)](https://travis-ci.org/cnobles/iGUIDE)
[![CircleCI](https://circleci.com/gh/cnobles/iGUIDE.svg?style=svg)](https://circleci.com/gh/cnobles/iGUIDE)
[![Documentation Status](https://readthedocs.org/projects/iguide/badge/?version=latest)](http://iguide.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/88088016.svg)](https://zenodo.org/badge/latestdoi/88088016)


Bioinformatic pipeline for processing iGUIDE and GUIDE-seq samples.

### Description
iGUIDE is a pipeline written in [snakemake](http://snakemake.readthedocs.io/) for processing and analyzing double-strand DNA break events. These events may be induced, such as by designer nucleases like Cas9, or spontaneous, as produced through DNA replication or ionizing radiation. A laboratory bench-side protocol accompanies this software pipeline, and can be found XXX. 

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

# Or include simulation test
cd path/to/iGUIDE
bash install.sh -t

# For help with install options:
cd path/to/iGUIDE
bash install.sh -h
```

### An Example Run
To perform a local test of running the iGUIDE informatic pipeline, run the below code after installing. This block first activates your conda environment, 'iguide' by default, and then creates a test directory within the analysis directory. The run information is stored in the run specific configuration file (config file). Using the '-np' flag with the snakemake call will perform a dry-run (won't actually process anything) and print the commands to the terminal, so you can see what snakemake is about to perform. Next, the test data can be moved to the input directory underneath the new test run directory or the path to the input data needs to be included in the config file. Then the entirety of processing can start. 

```
# If conda is not in your path ...
source ${HOME}/miniconda3/etc/profile.d/conda.sh

# Activate iguide environment
conda activate iguide

# After constructing the config file and having reference files (i.e. sampleinfo)
# You can check the samples associated with the run.
iguide list_samples configs/simulation.config.yml

# Create test analysis directory
# (The simulation configuration file is used by default and does not need to be specified)
iguide setup configs/simulation.config.yml -- -np
iguide setup configs/simulation.config.yml

# Process a simulation dataset
iguide run configs/simulation.config.yml -- -np
iguide run configs/simulation.config.yml -- --latency-wait 30
cat analysis/simulation/output/unique_sites.simulation.csv

# After run completion, generate a report in a different format than standard
iguide report analysis/simulation/output/edited_sites.simulation.rds \
  -c configs/simulation.config.yml \
  -o analysis/simulation/reports/report.simulation.pdf \
  -s sampleInfo/simulation.supp.csv \
  -t pdf

# When all finished and ready to archive / remove excess files, a minimal configuration
# can be achived with the 'clean' subcommand.
iguide clean configs/simulation.config.yml

# Or you realized you messed all input up and need to restart
iguide clean configs/simulation.config.yml --remove_proj

# Deactivate the environment
conda deactivate
```

### Changelog:

**v0.9.4 (January 30th, 2019)**

* Updated 'report' utility and formating
  + custom templates now accepted
  + included as subcommand, check with 'iguide report -h'
  + pdf and html options report 'nicely' even when printed from either
* Updated build to v0.9.1 to support new formating in report
* Included the 'clean' subcommand to reduce size of processed projects
  + after cleaning a project, only terminal data files will remain

**v0.9.3 (January 11th, 2019)**

* Added 'list_samples' subcommand to list samples within a project.
* Caught a few bugs and worked them out for smoother processing and reports.

**v0.9.2 (January 7th, 2019)**

* Modified test dataset to run tests quicker and implemented CirclCI checking.

**v0.9.1 (January 6th, 2019)**

* Fixed problematic install for first time conda installers.

**v0.9.0 (January 4th, 2019)**

* Initial release.
* Supports setup and analysis of GUIDE-seq and iGUIDE experiments.
* Documentation on [ReadTheDocs.io](https://iguide.readthedocs.io/en/latest/index.html).
