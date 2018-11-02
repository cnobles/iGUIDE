.. iGUIDE documentation master file, created by
   sphinx-quickstart on Fri Nov  2 14:34:12 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to iGUIDE's documentation
==================================

===========
Description
===========

Software pipeline for processing and analyzing double-strand DNA break events. 
These events may be induced, such as by designer nucleases like Cas9, or 
spontaneous, as produced through DNA replication or ionizing radiation. A 
laboratory bench-side protocol accompanies this software pipeline, and can be 
found XXX. 

Below, this readme gives the reader a overview of the pipeline, including how 
to install and process a sample dataset. Processing a sample data set is broken 
into three parts: 

1. developing a configuration file and sample information
2. setting up a run directory and acquiring the sequence data
3. initializing the pipeline and understanding the output

.. toctree::
   :maxdepth: 3
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


=======
Install
=======

To install iGUIDE, simply clone the repository to the desired destination::
  
  git clone https://github.com/cnobles/iGUIDE.git

Then initiate the install using the install script. If you would like the 
installed environment to be named something other than 'iguide', the new conda 
environment name can be provided to the 'install.sh' script as provided below::

  cd path/to/iGUIDE
  bash bin/install.sh

Or::

  cd path/to/iGUIDE
  bash bin/install.sh {env_name}

=====
Usage
=====

-----------------------------
Configuration Files - Configs
-----------------------------

Configuration files, or configs for short, contain both run-related and 
pipeline-related information. This is by design, for reproducibility it is 
easiest to have what was processed and how it was processed in the same location. 
There should be one config file for each sequencing run to be processed. Below 
is a brief summary of how to 'configure' your config file to your specific run.

Config files need to be named in the format '{RunName}.config.yml', where 
{RunName} is a parameter set within the config file for the run. For example, 
the default run configuration file is named 'simulation.config.yml', so the 
run name is 'simulation'.

Config files can be deposited anywhere in the users directory, but a dediacted 
directory has been included in the release of iGUIDE. For convienence, config 
files can be placed in iGUIDE/configs/.

For sample specific information, input is more easily placed in a sampleInfo 
file. See the section below regarding sample info files.

The path to the iGUIDE software will need to be included in the config file. 
Therefore, if processing may be done on multiple machines or platforms, make 
sure to update the file path of the install.

^^^^^^^^^^^^^^^^^^
Config File Layout
^^^^^^^^^^^^^^^^^^

Config files are in a yaml format, but are broken into two parts. The first 
contains run specific information that should be filled out by an individual 
familiar with the sequence data used in the laboratory bench-side protocol. The 
second part (below the divide) should be filled out by an individual familiar 
with the bioinformatic processing. Explanations of the different portions can 
be found below.

"""""""""""""""""""""""""""""""""
Config - Run Specific Information
"""""""""""""""""""""""""""""""""

"""""""""""""""""""""""""""""""
Config - Processing Information
"""""""""""""""""""""""""""""""

------------------------
Sample Information Files
------------------------

Sample information files (or sampleInfo files) contain information that may 
change from specimen to specimen. They need to contain at lease 3 columns of 
information: sample names, barcode 1, and barcode 2 sequences. Additionally, 
other parameters defined in the config file can be defined in the sampleinfo 
file if they change from specimen to specimen. 

Run specific config file will need to point to the sampleInfo files. For 
convienence, a directory can be found at iGUIDE/sampleInfo/ for depositing 
sampleInfo files.

SampleInfo files also need to have a specific naming format 
'{RunName}.sampleinfo.csv'.

------------------
Initializing a Run
------------------

Once the config and sampleInfo files have been configured, a run directory can 
be created using the command below where {RunName} is your run name::

  cd path/to/iGUIDE
  snakemake analysis/{RunName} --configfile configs/{RunName}.config.yml

The directory should look like this::
  
  >tree analysis/{RunName}
  analysis/{RunName}/
  ├── config.yml -> ../../configs/{RunName}.config.yml
  ├── input_data
  ├── logs
  ├── output
  ├── processData
  └── reports

Components within the run directory:

* config.yml - This is a symbolic link to the config file for the run
* input_data - Directory where input fastq.gz files are deposited
* logs - Directory containing log files from processing steps
* output - Directory containing output data from the analysis
* processData - Directory containing intermediate processing files
* reports - Directory containing output reports and figures

As a current convention, all processing is done within the analysis directory. 
The above command will create a file directory under the analysis directory for 
the run specified in by the config ('/iGUIDE/analysis/{RunName}'). At the end of 
this process, iGUIDE will give the user a note to deposit the input sequence 
files into the /analysis/{RunName}/input_data directory. Copy the fastq.gz files 
from the sequencing instrument into this directory.

Currently, iGUIDE needs each of the sequencing files (R1, R2, I1, and I2) for 
processing. If I1 and I2 are concatenated into the read names of R1 and R2, it 
is recommended the you run ``bcl2fastq ... --create-fastq-for-index-reads`` on 
the machine output directory to generate the I1 and I2 files. 

As iGUIDE has its own demultiplexing, it is recommend to not use the Illumina 
machine demultiplexing through input of index sequences in the SampleSheet.csv. 
See SampleSheet example in XXX. If sequence files are demultiplexed, they can be 
concatenated together into one file for each type of read using 'zcat'.

----------------
Processing a Run
----------------

Once the input_data directory has the required sequencing files, the run can be 
processed using the following command::

  cd path/to/iGUIDE/
  snakemake --configfile configs/{RunName}.config.yml

Snakemake offers a great number of resources for managing the processing through 
the pipeline. I recommend familiarizing yourself with the utility (XXX). Some helpful flags:

* [--configfile X] associate a specific configuration for processing, essential for processing
* [--cores X] multicored processing, specified cores to use by X
* [--nolock] process multiple runs a the same time, from different sessions
* [--notemp] keep all temporary files, otherwise removed
* [--keep-going] will keep processing if one or more job error out
* [-w X, --latency-wait X] wait X seconds for the output files to appear before erroring out

--------------
An Example Run
--------------

To perform a local test of running the iGUIDE informatic pipeline, run the below 
code after installing. This block first activates your conda environment, 
``iguide`` by default, and then creates a test directory within the analysis 
directory. The run information is stored in the run specific configuration file 
(config file). Using the ``-np`` flag with the snakemake call will perform a 
dry-run (won't actually process anything) and print the commands to the 
terminal, so you can see what snakemake is about to perform. Next, the test data 
is moved to the input directory underneath the new test run directory. Then the 
entirety of processing can start. Using the ``--dag`` flag and piping the output 
to ``dot -Tsvg`` will generate a vector graphic of the directed acyclic graph 
(dag) workflow that snakemake will follow given the data provided::

  # Test script
  PREFIX=${HOME}/miniconda3
  export PATH=${PATH}:${PREFIX}/bin
  source activate iguide

  # Create test analysis directory
  snakemake analysis/simulation --configfile configs/simulation.config.yml -np
  snakemake analysis/simulation --configfile configs/simulation.config.yml

  # Move test sequence files to analysis directory
  cp tests/Data/Undetermined_S0_L001_* analysis/simulation/input_data/

  # Generate test DAG graph
  snakemake --configfile configs/simulation.config.yml -np
  snakemake --configfile configs/simulation.config.yml --dag | dot -Tsvg > analysis/simulation/reports/simulation.dag.svg
  snakemake --configfile configs/simulation.config.yml --latency-wait 30
  cat analysis/simulation/output/unique_sites.simulation.csv

---------
Uninstall
---------

To uninstall iGUIDE, the user will need to remove the environment and the 
directory.

To remove the environment and channels used with conda::

  cd path/to/iGUIDE
  bash bin/uninstall.sh

Or::

  cd path/to/iGUIDE
  bash bin/uninstall.sh {env_name}

If the user would rather remove the environment created for iGUIDE, it is 
recommended directly use conda. This will leave the channels within the conda 
config for use with other conda configurations::

  conda env remove -n iguide

Or::

  conda env remove -n {env_name}

To remove the iGUIDE directory and conda, the following two commands can be 
used::

  # Remove iGUIDE directory and software
  rm -r path/to/iGUIDE

  # Remove conda
  rm -r path/to/miniconda3
