.. _quickstart:

.. contents::
   :depth: 2

------------------
Initializing a Run
------------------

Once the config and sampleInfo files have been configured, a run directory can 
be created using the command below where {ConfigFile} is the path to your configuration file::

  cd path/to/iGUIDE
  iguide setup --configfile {ConfigFile}

The directory should look like this (RunName is specified in the ConfigFile}::
  
  >tree analysis/{RunName}
  analysis/{RunName}/
  ├── config.yml -> {path to ConfigFile}
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
  iguide run --configfile {ConfigFile}

Snakemake offers a great number of resources for managing the processing through 
the pipeline. I recommend familiarizing yourself with the utility (XXX).
Here are some helpful flags that can be passed to snakemake by appending to the iguide command:

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
  # (The simulation configuration file is used by default and does not need to be specified)
  iguide setup -np
  iguide setup

  # Generate test DAG graph
  iguide run -np
  iguide run --dag | dot -Tsvg > analysis/simulation/reports/simulation.dag.svg
  iguide run --latency-wait 30
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
