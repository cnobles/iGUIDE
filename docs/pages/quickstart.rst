.. _quickstart:

.. contents::
   :depth: 2



Initializing a Run
------------------

Once the config and sampleInfo files have been configured, a run directory can 
be created using the command below where {ConfigFile} is the path to your configuration file::

  cd path/to/iGUIDE
  iguide setup {ConfigFile}

The directory should look like this (RunName is specified in the ConfigFile}::
  
  > tree analysis/{RunName}
  analysis/{RunName}/
  ├── config.yml -> {path to ConfigFile}
  ├── input_data
  ├── logs
  ├── output
  ├── process_data
  └── reports

Components within the run directory:

* config.yml - This is a symbolic link to the config file for the run
* input_data - Directory where input fastq.gz files can be deposited
* logs - Directory containing log files from processing steps
* output - Directory containing output data from the analysis
* process_data - Directory containing intermediate processing files
* reports - Directory containing output reports and figures

As a current convention, all processing is done within the analysis directory. 
The above command will create a file directory under the analysis directory for 
the run specified in by the config ('/iGUIDE/analysis/{RunName}'). At the end of 
this process, iGUIDE will give the user a note to deposit the input sequence 
files into the /analysis/{RunName}/input_data directory. Copy the fastq.gz files 
from the sequencing instrument into this directory if you do not have paths to
the files specified in the config file.

Currently, iGUIDE needs each of the sequencing files (R1, R2, I1, and I2) for 
processing since it is based on a dual barcoding scheme. If I1 and I2 are 
concatenated into the read names of R1 and R2, it is recommended the you run 
``bcl2fastq ... --create-fastq-for-index-reads`` on the machine output 
directory to generate the I1 and I2 files. 

As iGUIDE has its own demultiplexing, it is recommend to not use the Illumina 
machine demultiplexing through input of index sequences in the SampleSheet.csv. 
See SampleSheet example in XXX. If sequence files are demultiplexed, they can be 
concatenated together into one file for each type of read using 'zcat'.


List Samples for a Run
----------------------

As long as the config and sampleInfo files are present and in their respective 
locations, you can get a quick view of what samples are related to the project.
Using the 'list_samples' subcommand will produce an overview table on the 
console or write the table to a file (specified by the output option). 
Additionally, if a supplemental information file is associated with the run, the
data will be combined with the listed table.::

  > iguide list_samples configs/simulation.config.yml
  
  Specimen Info for : simulation.

   specimen   replicates       gRNA        nuclease
  ---------- ------------ --------------- ----------
     iGXA         1            TRAC         Cas9v1
     iGXB         1        TRAC;TRBC;B2M    Cas9v1
     iGXD         1             NA            NA


Processing a Run
----------------

Once the input_data directory has the required sequencing files, the run can be 
processed using the following command::

  cd path/to/iGUIDE/
  iguide run {ConfigFile}

Snakemake offers a great number of resources for managing the processing through 
the pipeline. I recommend familiarizing yourself with the utility 
(https://snakemake.readthedocs.io/en/stable/). Here are some helpful snakemake
options that can be passed to iGUIDE by appending to the iguide command after ``--``:

* ``[--configfile X]`` associate a specific configuration for processing, essential for processing but already passed in by ``iguide``.
* ``[--cores X]`` multicored processing, specified cores to use by X.
* ``[--nolock]`` process multiple runs a the same time, from different sessions.
* ``[--notemp]`` keep all temporary files, otherwise removed.
* ``[--keep-going]`` will keep processing if one or more job error out.
* ``[-w X, --latency-wait X]`` wait X seconds for the output files to appear before erroring out.
* ``[--restart-times X]`` X is the number of time to restart a job if it fails. Defaults to 0, but is used in ``iguide`` to increase memory allocation.
* ``[--resources mem_mb=X]`` Defined resources, for ``iguide`` the mem_mb is the MB units to allow for memory allocation to the whole run. For HPC, this can be coupled with ``--cluster-config`` to request specific resources for each job.
* ``[--rerun-incomplete, --ri]`` Re-run all jobs that the output is recognized as incomplete, useful if your run gets terminated before finishing.
* ``[--cluster-config FILE]`` A JSON or YAML file that defines wildcards used for HPC.


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
entirety of processing can start.::

  # After constructing the config file and having reference files (i.e. sampleinfo)
  # You can check the samples associated with the run.
  
  iguide list_samples configs/simulation.config.yml

  # Create test analysis directory
  # (The simulation configuration file is used by default and does not need to be specified)
  
  iguide setup configs/simulation.config.yml

  # Process a simulation dataset

  iguide run configs/simulation.config.yml -- -np
  iguide run configs/simulation.config.yml -- --latency-wait 30
  zcat analysis/simulation/output/unique_sites.simulation.csv.gz

  # Processing will complete with a report, but if additional analyses are required,
  # you can reevaluate the 'incorp_sites' object. Multiple objects can be evaluated
  # together, just include the run files.

  iguide eval analysis/simulation/output/incorp_sites.simulation.rds \
    -o analysis/simulation/output/iguide.eval.simulation.test.rds \
    -s sampleInfo/simulation.supp.csv

  # After evaluation, generate a report in a different format than standard.
  # Additionally the evaluation and report generation step can be combined using 
  # config file(s) as inputs for the 'report' subcommand.

  iguide report -e analysis/simulation/output/iguide.eval.simulation.test.rds \
    -o analysis/simulation/reports/report.simulation.pdf \
    -s sampleInfo/simulation.supp.csv \
    -t pdf

  # When you are all finished and ready to archive / remove excess files, a minimal configuration
  # can be achived with the 'clean' subcommand.

  iguide clean configs/simulation.config.yml

  # Or you realized you messed up all the input and need to restart

  iguide clean configs/simulation.config.yml --remove_proj


Uninstall
---------

To uninstall iGUIDE, the user will need to remove the environment and the 
directory.

To remove the environment and channels used with conda::

  cd path/to/iGUIDE
  bash etc/uninstall.sh

Or::

  cd path/to/iGUIDE
  bash etc/uninstall.sh {env_name}

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
