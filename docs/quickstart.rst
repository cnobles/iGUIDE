.. _quickstart:

Quickstart Guide
================

.. contents::
   :depth: 2

Install
*******

To install iGUIDE, simply clone the repository to the desired destination.::

  git clone https://github.com/cnobles/iGUIDE.git

Then initiate the install using the install script. If you would like the 
installed environment to be named something other than 'iguide', the new conda 
environment name can be provided to the ``install.sh`` script as provided 
below.::

  cd path/to/iGUIDE
  bash install.sh

Or::

  cd path/to/iGUIDE
  bash install.sh -e {env_name}
  
Additionally, help information on how to use the ``install.sh`` can be accessed
by::

  bash install.sh -h


Setup a Run
***********

Once the config and sampleInfo files have been configured, a run directory 
can be created using the command below where {ConfigFile} is the path to your 
configuration file::

  cd path/to/iGUIDE
  iguide setup {ConfigFile}

The directory should look like this (RunName is specified in the ConfigFile)::
  
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

iGUIDE typically uses each of the sequencing files (R1, R2, I1, and I2) for 
processing since it is based on a dual barcoding scheme. If I1 and I2 are 
concatenated into the read names of R1 and R2, it is recommended the you run 
``bcl2fastq ... --create-fastq-for-index-reads`` on the machine output 
directory to generate the I1 and I2 files. 

As iGUIDE has its own demultiplexing, it is recommend to not use the Illumina 
machine demultiplexing through input of index sequences in the SampleSheet.csv.
If your sequence data has already been demultiplexed though, please see the 
:ref:`usage` for setup instructions.


List Samples in a Run
*********************

As long as the config and sampleInfo files are present and in their respective 
locations, you can get a quick view of what samples are related to the project.
Using the ``iguide list_samples`` command will produce an overview table on 
the console or write the table to a file (specified by the output option).
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
****************

Once the input_data directory has the required sequencing files, the run can be 
processed using the following command::

  cd path/to/iGUIDE/
  iguide run {ConfigFile}

Snakemake offers a great number of resources for managing the processing through 
the pipeline. I recommend familiarizing yourself with the utility 
(https://snakemake.readthedocs.io/en/stable/).


An Example Workflow
*******************

To perform a local test of running the iGUIDE informatic pipeline, run the below 
code after installing. This block first activates your conda environment, 
'iguide' by default, and then creates a test directory within the analysis 
directory. The run information is stored in the run specific configuration file 
(config file). Using the ``-np`` flag with the snakemake call will perform a 
dry-run (won't actually process anything) and print the commands to the 
terminal, so you can see what snakemake is about to perform. Then the entirety 
of processing can start.::

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

  # When you are all finished and ready to archive / remove excess files, a minimal structure
  # can be achieved with the 'clean' subcommand.

  iguide clean configs/simulation.config.yml

  # Or you realized you messed up all the input and need to restart

  iguide clean configs/simulation.config.yml --remove_proj

  # Deactivate the environment

  conda deactivate


Reviewing Results
*****************

The output reports from a run are deposited under 
``analysis/{RunName}/reports``. For more informtion on output files, see :ref:`usage`!
