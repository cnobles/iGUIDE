## iGUIDE - improved Genome-wide Unbiased Identification of Double-strand DNA break Events
[![Build Status](https://travis-ci.org/cnobles/iGUIDE.svg?branch=master)](https://travis-ci.org/cnobles/iGUIDE)
[![CircleCI](https://circleci.com/gh/cnobles/iGUIDE.svg?style=svg)](https://circleci.com/gh/cnobles/iGUIDE)
[![Documentation Status](https://readthedocs.org/projects/iguide/badge/?version=latest)](http://iguide.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/88088016.svg)](https://zenodo.org/badge/latestdoi/88088016)


Bioinformatic pipeline for processing iGUIDE and GUIDE-seq samples.

### Description
iGUIDE is a pipeline written in [snakemake](http://snakemake.readthedocs.io/) for processing and analyzing double-strand DNA break events. These events may be induced, such as by designer nucleases like Cas9, or spontaneous, as produced through DNA replication or ionizing radiation. A laboratory bench-side protocol accompanies this software pipeline, and can be found [**https://doi.org/10.1186/s13059-019-1625-3**](https://doi.org/10.1186/s13059-019-1625-3). 

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

iguide setup configs/simulation.config.yml

# Process a simulation dataset

iguide run configs/simulation.config.yml -- -np
iguide run configs/simulation.config.yml -- --latency-wait 30

# Processing will complete with several reports, but if additional analyses are required,
# you can re-evaluate a run by its config file. Multiple runs can be evaluated together, 
# just include multiple config files.

iguide eval configs/simulation.config.yml \
  -o analysis/simulation/output/iguide.eval.simulation.test.rds \
  -s sampleInfo/simulation.supp.csv

# After evaluation, generate a report in a different format than standard.
# Additionally the evaluation and report generation step can be combined using 
# config file(s) as inputs for the 'report' subcommand (using the -c flag instead of -e).

iguide report -e analysis/simulation/output/iguide.eval.simulation.test.rds \
  -o analysis/simulation/reports/report.simulation.pdf \
  -s sampleInfo/simulation.supp.csv \
  -t pdf

# When you are all finished and ready to archive / remove excess files, a minimal configuration
# can be achieved with the 'clean' subcommand.

iguide clean configs/simulation.config.yml

# Or you realized you messed up all the input and need to restart

iguide clean configs/simulation.config.yml --remove_proj

# Deactivate the environment

conda deactivate
```

### Changelog:

**v0.9.9 (August 9th,2019) - Additional updates**

* Implemented support for BWA aligner
* Added tools (samqc) for working with other SAM/BAM output aligners as well
* Switched iguide support code to iguideSupport R-package and added unit tests
* Fixed bugs related to quoted table inputs (csv/tsv)
* Implemented a method to skip demultiplexing, see documentation for setup
* Resoved a number of issues identified, check GitHub for history!

**v0.9.9 (June 10th, 2019)**

* Modified the assimilate + evaluate workflow
  + Assimilate now only includes reference genome data, meaning a cleaner intermediate file
  + Evaluate will now handle ref. gene sets and further analysis
  + This increases the modularity and consistancy of the workflow
* Revised the iGUIDE Report format to be more informational and clearer
* Revised a bit of the workflow to make reprocessing smoother
* Updated BLAT coupling script to be more memory efficient
* Fixed TravisCI testing!
* Changed stat workflow, now restarting analysis won't init a total reproc.

**v0.9.8 (April 19th, 2019)**

* iGUIDE can now support non-Cas9 nucleases as well!
  + Implemented nuclease profiles into configs
  + Updated assimilation, evaluation, and reporting scripts
* Added default resources to allow simpler HPC processing
* Included flexible system for identifying on-target sites
  + Config can accept a range rather than a single site
  + Acceptable notation: chr4:+:397-416 and chr3:*:397
* Changed build nomenclature from v0.9.3 to b0.9.3
  + So as not to confuse with version
* Added 'summary' subcommand to generate a consise text-based report
  + Working in the same manner as 'report', can generate from config(s) or eval file
* Added short stats-based report to be produced at the end of processing
* Additional bugfixes.

**v0.9.7 (March 6th, 2019)**

* Hotfix to workflow.
* Changed 'setup' subcommand to python script based rather than snakemake.
* Changed file organization.

**v0.9.6 (March 5th, 2019)**

* Introduced process workflow steps: assimilate and evaluate
  + Assimilate aligned data and compare with targeting sequences
    + Core data object that can be combined across runs / projects
  + Evaluated data incorporates reference data and statistical models
    + A staple data object for reports and can be constructed from multiple runs
* Included new subcommands 'eval' and modified 'report'
  + report from either config(s) or eval dataset
* Cleaned up file structure
* Updated documentation in code and docs.
* Implemented accuracy and retention checks with simulation dataset.
* Updated simulation dataset with larger set to test analysis.

**v0.9.5 (February 19th, 2019)**

* Updated demultiplexing to be more efficient and better HPC compatible.
* Added RefSeq Extended* reference gene sets
  + 'ext' includes curated, predicted, and other RefSeq sets
  + 'ext.nomodel' includes only curated and other RefSeq sets
* Incorporated resource allocation for job dependent memory consumption
  + Works great with HPC to specify memory requirements
* Streamlined input for report generation by only requiring config(s)


**v0.9.4 (January 30th, 2019)**

* Updated 'report' utility and formating
  + custom templates now accepted
  + included as subcommand, check with 'iguide report -h'
  + pdf and html options report 'nicely' even when printed from either
* Updated build to v0.9.2 to support new formating in report
* Builds are constructed from spec files rather than yaml requirements
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
