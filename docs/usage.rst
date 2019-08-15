.. _usage:

User Guide
==========

.. contents::
   :depth: 3



Nomenclature / Semantics
************************

Walk through of an generalized experiment with consistent terminology

Specimen vs Sample

Naming restrictions on Specimens vs Samples

Subcommands
***********

Overview of subcommands: setup, run, eval, report, summary, clean

.. figure:: figures/iguide_subcmd_fig.pdf
   :figwidth: 50%
   :alt: iGUIDE Subcommands: setup, run, eval, report, summary, clean.
   :align: center


Workflow
********

A generalized workflow, example with associated figures, primary then auxiliary

.. figure:: figures/iguide_prime_workflow_fig.pdf
   :figwidth: 90%
   :alt: Primary iGUIDE Workflow: setup, transfer/link data, run.
   :align: center
   



.. figure:: figures/iguide_aux_workflow_fig.pdf
   :figwidth: 90%
   :alt: Auxiliary iGUIDE Workflow: eval, report, summary, clean.
   :align: center
   




Requirements
************

- A relatively-recent Linux computer with more than 2Gb of RAM

We do not currently support Windows or Mac. (You may be able to run this on
Windows using the [WSL](https://docs.microsoft.com/en-us/windows/wsl/about), but
it has not been tested.




Installing
**********

To install iGUIDE, simply clone the repository to the desired destination.

.. code-block:: shell
  
  git clone https://github.com/cnobles/iGUIDE.git

Then initiate the install using the install script. If you would like the 
installed environment to be named something other than 'iguide', the new conda 
environment name can be provided to the ``install.sh`` script as provided below.

.. code-block:: shell

  cd path/to/iGUIDE
  bash install.sh

Or:

.. code-block:: shell

  cd path/to/iGUIDE
  bash install.sh -e {env_name}
  
Additionally, help information on how to use the ``install.sh`` can be accessed
by:

.. code-block:: shell

  bash install.sh -h
  
  
Testing
-------

How to run the test case and what to expect
Simulation dataset builder



Updating
--------

how to use the "install.sh" script for updating components (env, lib, pkg)


Uninstalling
------------

To uninstall iGUIDE, the user will need to remove the environment and the 
directory.

To remove the environment and channels used with conda::

  cd path/to/iGUIDE
  bash etc/uninstall.sh

Or::

  cd path/to/iGUIDE
  bash etc/uninstall.sh {env_name}

If the user would rather remove the environment created for iGUIDE, it is 
recommended to use conda. This will leave the channels within the conda 
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
(https://snakemake.readthedocs.io/en/stable/). Here are some helpful snakemake
options that can be passed to iGUIDE by appending to the iguide command after 
``--``:

* ``[--configfile X]`` associate a specific configuration for processing, 
  essential for processing but already passed in by ``iguide``.
* ``[--cores X]`` multicored processing, specified cores to use by X.
* ``[--nolock]`` process multiple runs a the same time, from different sessions.
* ``[--notemp]`` keep all temporary files, otherwise removed.
* ``[--keep-going]`` will keep processing if one or more job error out.
* ``[-w X, --latency-wait X]`` wait X seconds for the output files to appear 
  before erroring out.
* ``[--restart-times X]`` X is the number of time to restart a job if it fails. 
  Defaults to 0, but is used in ``iguide`` to increase memory allocation.
* ``[--resources mem_mb=X]`` Defined resources, for ``iguide`` the mem_mb is the
  MB units to allow for memory allocation to the whole run. For HPC, this can be
  coupled with ``--cluster-config`` to request specific resources for each job.
* ``[--rerun-incomplete, --ri]`` Re-run all jobs that the output is recognized 
  as incomplete, useful if your run gets terminated before finishing.
* ``[--cluster-config FILE]`` A JSON or YAML file that defines wildcards used 
  for HPC.




Config Files
************

Configuration files, or configs for short, contain both run-related and 
pipeline-related information. This is by design. For reproducibility it is 
easiest to have what was processed and how it was processed in the same 
location. There should be one config file for each sequencing run to be 
processed. Below is a brief summary of how to 'configure' your config file to 
your specific run.

Config files need to be named in the format '{RunName}.config.yml', where 
``{RunName}`` is a parameter set within the config file for the run. For 
example, the default run configuration file is named ``simulation.config.yml``, 
so the run name is ``simulation``.

Config files can be deposited anywhere in the users directory, but a dediacted 
directory has been included in the release of iGUIDE. For convienence, config 
files can be placed in ``iGUIDE/configs/``.

For sample specific information, input is more easily placed in a sampleInfo 
file. See the included section regarding sample info files.

File Layout
-----------

Config files are in a ``yaml`` format, but are broken into two parts. The first 
contains run specific information that should be filled out by an individual 
familiar with the sequence data used in the laboratory bench-side protocol. 
Additionally, they should be aware of the biochemistry related to the enzymes 
and sequences they are using.

The second part (below the divide ``----``) should be filled out by an 
individual familiar with the bioinformatic processing. Explanations of the 
different portions can be found in the following pages.









Run Specific Information
------------------------

Run configuration
"""""""""""""""""

``Run_Name``
  This is the name of the sequencing run, and should only contain alpha-numeric
  characters. Underscores (``_``) and dashes (``-``) are also allowed within the
  run name parameters. Other symbols should not be included, such as a dot 
  (``.``). The run name is further used by the software to link files and 
  directories together, so it will need to be consistent whenever it is used.
  Examples include: iGUIDE_190201_B6V99, 181213_PD1_T-cell_exp.
  
``Sample_Info``
  This is a file path to the sample information file. It can either be an 
  absolute file path or relative file path. If the file path is relative though,
  it will need to be relative to the Snakefile used by the iGUIDE software. For
  more information about this file, please see the Sample Information page.
  
``Supplemental_Info``
  Similar to ``Sample_Info``, this is a file path to a supplementary file which
  can contain information related to experimental parameters or patient 
  information. This will be used during the report output, which will group
  samples with identical parameters. The format for this file is quite loose, 
  and it only requires a single column ``Specimen``, which should match the 
  names of specimens in the sample information file. For more information about
  this file, please see the Supplemental Information page. If no file is to be 
  used, set the value for this parameter to ``"."`` and make sure to set the 
  ``suppFile`` in the run protion of the config to ``FALSE``. 
  
``Ref_Genome``
  This is a designation for the reference genome to used during processing. The
  genome will need to be included in the R libraries through BioConductoR prior
  to running the software. The human genome draft ``hg38`` is included by 
  default. Please see information on the BioConductoR package 'BSgenome' for 
  installing alternative genomes.
  
``Ref_Genome_Path``
  This is the file path (following the same workflow as the ``Sample_Info`` 
  parameter) to a reference genome file, if one is already available in a fasta
  format.
  
``Aligner``
  Options include either 'blat' or 'bwa', though at this time, only 'blat' is 
  supported. Future versions of iGUIDE will support other alignment softwares.
  
``UMItags``
  This is a logical parameter indicating whether to use unique molecular indices
  (UMI) sequence tags ('TRUE') or to only use unique fragments lengths (see
  `SonicAbundance <https://doi.org/10.1093/bioinformatics/bts004>`) to quantify
  abundances of unique observations.
  
  
Sequence files
""""""""""""""

``Seq_Path``
  This is the file path to the sequence files. Rather than repeating the path
  for each below, just include the path to the directory containing the files.

``R1 / R2 / I1 / I2``
  These parameters should be the file names of the sequence files to be 
  analyzed by the iGUIDE software. It is recommened to pass complete sequencing
  files to iGUIDE rather than demultiplexing prior to analysis.

``Demulti_Dir``
  Path to the directory containing demultiplexed sequence data. This is still
  under development and may present with bugs.


SampleInfo formating
""""""""""""""""""""

``Sample_Name_Column``
  This is the name of the column in the sample information file which contains 
  identifiable information about samples. An appropriate format for the sample 
  names is "{specimen}-{rep}" where 'specimen' is an alpha-numeric designator 
  for the specimen and 'rep' is a numeric identifier for technical or biological 
  replicates, separated by a dash (``-``). Replicates will be pooled during the
  final analysis, so if you want them to be separate in the report, make sure
  you give each specimen a different identifier. For example, iGSP0002-1 and
  iGSP0002-2, will be pooled together for the report and analysis, but 
  iGSP0002-1 and iGSP0003-1 will not. These names will be used in naming files,
  so do not include any special characters that will confuse file managment. 
  Try to stick to common delimiters, such as "-", "_", ".". A good practice is
  to put specimen identifiers at the beginning, replicate identifiers at the end
  following a "-", and anything else descriptive in the middle. For example, 
  iGSP0002-neg-1, can specify the orientation the sample was processed with.


Sequence information
""""""""""""""""""""

``R{1/2}_Leading_Trim``
  Sequence to be removed from the 5' or beginning of the R1 or R2 sequences. 
  Commonly a linker or fixed sequence that is part of the priming scheme during
  amplification. If no sequence should be removed, just include ``"."``. If the
  sequence is sample or specimen specific, it can be included in the sample 
  information file and indicated in these fields as ``"sampleInfo:{column}"``, 
  where 'column' is the column name with the data in the sample information 
  file.

``R{1/2}_Overreading_Trim``
  Similar to the ``Leading_Trim`` parameters, these parameters indicate the 
  sequence that should be removed from the 3' or end of the reads if it is 
  present. Again, if no sequence should be removed, use a ``"."`` or if the data
  is present in the sample information file, ``"sampleInfo:{column}"``.

``R2_Leading_Trim_ODN``
  This is a key parameter difference between iGUIDE and the original GUIDEseq
  method. This parameter indicates the sequence that is part of the dsODN but is
  **not** primed against. This sequence should directly follow the 
  ``R2_Leading_Trim`` sequence and should be a reverse complement of the 
  beginning of the ``R1_Overreading_Trim`` sequence if the iGUIDE dsODN is being 
  used. For GUIDEseq, simply include ``"."``, or if you have multiple sequences,
  then specify in the sample information file as ``"sampleInfo:{column}"``. 


Target sequence information
"""""""""""""""""""""""""""

``Target_Sequences``
  This parameter specifies the target sequences, **not including** the PAM 
  sequences for guide RNAs. An acceptable input format would be 
  ``{target_name} : "{sequence}"`` (i.e. ``B2M.3 : "GAGTAGCGCGAGCACAGCTANGG"``) 
  and additional target sequences can be included, one per line, and each 
  indented at the same level. The input format of 
  ``{target_name} : {target_seq}`` needs to be maintained for proper function. 
  The 'target_name' in this situation will need to match the 'target_name' used 
  in the ``On_Target_Sites`` and ``Treatment`` parameters. 'target_name' should 
  follow a common format, and use standard delimiters, such as "-", "_", and 
  ".". For example: ``B2M.3``, ``TRAC.1.5``, ``TruCD33v5``.

``On_Target_Sites``
  This parameter indicates the specific location for editing by the target 
  enzyme. There should be one line for each on-target site, even if there are 
  more than one on-target sites for a given target sequence. Typically the input
  format should follow ``{target_name} : "{seqname}:{+/-}:{position}"``, where 
  'target_name' matches the name of the given target sequence, and if multiple 
  on-target sites exist, then the names can be expanded using a 
  ``{target_name}'#`` notation. Additionally, the notation can be expanded to
  ``{target_name} : "{seqname}:{+/-/*}:{min.position}-{max.position}"``, where
  '*' indicates either orientation and 'min.position' and 'max.position' 
  represent the numerical range for the on-target site. The value for each 
  on-target site specifies the location or genomic coordinates of nuclease 
  activity. The 'seqname' indicates the chromosome or sequence name, an 
  orientation of '+' or '-' is given to the location depending on the editing 
  orientation (in line with positional numbering is '+' and opposite is '-', 
  unknown or both is '*'), and the 'position' or 'min/max.position' indicates 
  the nucleotide(s) of editing. For Cas9, the position of editing is commonly 
  between the 3rd and 4th nucleotide from the 3' end of the targeting sequence 
  (not including the PAM). Being off by a nucleotide or so will not cause any 
  problems. Example below.
  
  .. code-block:: shell
  
    On_Target_Sites :
      TRAC.5 : "chr14:+:22547664"
      TRBC.4'1 : "chr7:+:142792020"
      TRBC.4'2 : "chr7:+:142801367"
      PD1.3 : "chr2:-:241858808"
      TRAC.3.4 : "chr14:-:22550616-22550625"
      B2M.3 : "chr15:*:44711569-44711570"
      CIITA.15.1 : "chr16:+:10916399"


Specimen target treatment
"""""""""""""""""""""""""

``Treatment``
  This parameter indicates how samples were treated. If samples were all treated
  differently, then this information can be included in the sample information
  file as ``all : "sampleInfo:{column}"`` where 'column' is the name of the 
  column with the information. If a single sample was treated with more than one
  target sequence, then delimit multiple target names by a semicolon (``;``), 
  i.e. ``all : "B2M;TRAC;TRBC"``. Additionally, each specimen can be indicated 
  individually on a new line. Only specimen names should be given here and  
  provided individually, not sample identifiers. This means that if your sample
  names follow the suggested format, "{specimen}-{replicate}", you would only 
  specify the "{specimen} : {treatment}" underneath this parameter.


Specimen nuclease treatment

``Nuclease``
  Similar to target treatment above, this parameter dictates which nuclease(s)
  where used on the specimens. This refers to the class of nuclease, such as
  Cas9 or Cpf1, which behave differently when they edit DNA. Notation can follow
  the same as above, if all specimens were treated with the same class of
  nuclease, then just specify 'all : "{nuclease_profile}"', or list out by
  specimen. Additionally you can specify the column in sampleInfo in the same
  format as above. Currently, iGUIDE does not support processing for specimens
  with multiple classes of nuclease profiles. Only one profile can be specified
  per specimen.
  
  
  
  
  
  
  
  
  
  
  
  
Processing Information
----------------------

Below are parameters that are used to process the large amount of data, such as
setting memory suggestions if resources are specified or parameters for sequence
alignments. While these figues may not be relevant to the bench scientist, they
are particulars for computational scientists. 

Resource management is not required, but it can help when using HPC or limiting
jobs. You are encouraged to spend some time optimizing if you would like, these
parameters work out well on the designers platform.


iGUIDE configuration
""""""""""""""""""""

``Read_Types``
  This parameter should include which read types will be used in the analysis,
  i.e. ``["R1", "R2", "I1", "I2"]``. This follows a list notation is Python. If
  only single barcoding or some other method is employed and a read type is not
  included, simply leave it out of the example.

``Genomic_Reads``
  This parameter is similar to the ``Read_Types`` but only indicates which reads
  contain genomic information rather than indexing.


Memory Management
"""""""""""""""""

``defaultMB / demultiMB / trimMB / filtMB / consolMB / alignMB / qualCtrlMB / assimilateMB / evaluateMB / reportMB``
  Controls the amount of memory allocated to each of these processes during 
  snakemake processing. While working on a server or multicored machine, these
  parameters will work internally to help schedule jobs. Each value will act as
  an upper limit for the amount of MB of RAM to expect the process to take, and 
  schedule jobs appropriately using the ``--resources mem_mb={limitMB}`` flag with
  snakemake. During HPC use, these parameters can be combined with the cluster config
  to schedule specific memory requirements for jobs. Additionally, if the 
  ``--restart-times {x}`` is used where "x" is the number of times to restart a job
  if it fails, then the amount of memory for the job will increase by a unit of the 
  parameter. For example, if a trimming job fails because it runs out of memory, then
  restarting the job will try to allocate 2 times the memory for the second attempt.
  All parameters should be in megabytes (MB).


Demultiplexing parameters
"""""""""""""""""""""""""

``skipDemultiplexing``
  Logical (either TRUE or FALSE) to indicate if demultiplexing should be carried
  out. If TRUE, sequence files (*.fastq.gz) need to be placed or linked in the 
  input_data directory of an existing project directory (as with ``iguide setup``),
  one sequence file for each type (R1, R2, I1, I2). These need to be identified
  in the "Run" portion of the config file. If FALSE, then demultiplexed files need
  to be included in the input_data directory of an existing project directory. The
  files need to be appropriately named, in the format of ``{sampleName}.{readtype}.fastq.gz``,
  where ``sampleName`` matches the 'sampleName' column found in the associated 'sampleInfo'
  file, and ``readtype`` is R1, R2, I1, or I2. If ``UMItags`` is ``FALSE``, then only R1 and R2
  file types are required for analysis, if ``UMItags`` is ``TRUE``, then I2 is a
  required file type as well.

``barcode{1/2}Length``
  Integer values indicating the number of nucleotides in the barcodes or 
  indexing sequences.

``barcode{1/2}``
  Character values (i.e. ``"I1"``) indicating which reads to find the associated
  indexing information for demultiplexing.

``bc{1/2}Mismatch``
  An integer value indicating the number of tolarated mismatches in the barcode
  sequences for either barcode 1 or 2.


Sequence trimming
"""""""""""""""""

``R{1/2}leadMismatch``
  Integer values indicating the number of allowed mismatches in either R1 or R2
  leading sequence trimming. Recommend to set to less than 10% error.

``R2odnMismatch``
  Integer value indicating the number of allowed mismatches in the unprimed 
  ODN sequence, typically should be set to 0.

``R{1/2}overMismatch``
  Integer values indicating the number of allowed mismatches in either R1 or R2
  overreading trimming. This is converted into a percent matching and should be
  thought of as a number of mismatches allowed out of the total length of the 
  overreading trim sequence. 

``R{1/2}overMaxLength``
  Searching for overread trimming in sequences can be time consuming while not
  producing different results. For this the total length of searched for 
  sequences can be limited here. For example, if ``ATGCGTCGATCGTACTGCGTTCGAC`` 
  is used as the overreading sequence, and 5 mismatches are allowed, then the 
  tolerance will be 5/25 or 80% matching, but only the first 20 nucleotides of
  the sequence will be aligned for overtrimming, ``ATGCGTCGATCGTACTGCGT``. With
  an 80% matching requirement, 16 out of 20 nucleotides will need to align for
  overread trimming to be initiated.

Binning
"""""""

``bins``
  A number of bins to separate filtered sequences into for higher parallel 
  processing. The increasing the number of bins can help spread out the work
  required for processing to keep memory requirements lower.
  
``level``
  A number indicating the number of reads that should be targeted for each bin.
  Bins will be filled to the level amount, leaving remaining bins empty if 
  previous bins contain all the reads. Additionally, if all bins will
  "overflow", then reads will be evenly distributed across the number of bins.

Reference Alignment
"""""""""""""""""""

``BLATparams``
  A character string to be included with the BLAT call. A suggested example has
  been provided in the simulation config file. For options, please see the BLAT 
  help options by typing ``blat`` into the commandline after activating 
  ``iguide``.

``BWAparams``
  A character string to be inclued with the BWA call. A suggested example has
  been provided in the simulation config file. For options, please see BWA help
  by typing ``bwa mem`` into the commandline after activating ``iguide``.


Post-alignment filtering
""""""""""""""""""""""""

``maxAlignStart``
  Integer value indicating the number of nucleotides at the beginning of the 
  alignment that will be allowed to not align. Another way of thinking of this
  is the maximum start position on the query rather than the target reference.
  A default value of 5 means that the alignment needs to start in the first 5 
  nucleotides or the alignment is discarded during quality control filtering.

``minPercentIdentity``
  This is a value between 0 and 100 indicating the minimum global percent 
  identity allow for an alignment. If an alignment has less, then it is 
  discarded during quality control filtering.

``{min/max}TempLength``
  Specify the minimum (min) and maximum (max) template length expected. Joined
  alignments between R1 and R2 the are outside of this range are considered
  artifacts and are discarded or classified as chimeras.


Post-processing
"""""""""""""""

``refGenes / oncoGeneList / specialGeneList``
  These are special reference files in either text or BioConductoR's 
  GenomicRanges objects. They can be in an '.rds' format or table format 
  ('.csv' or '.tsv'). The ``file`` parameter should indicate the file path to
  the file (relative paths should be relative to the SnakeFile), and the 
  ``symbolCol`` parameter should indicate the column in the data object which 
  contains the reference names to be used in the analysis.
  
``maxTargetMismatch``
  The maximum number of mismatches between the reference genome and target
  sequence allowed for consideration to be a target matched incorporation 
  site. This is an integer value and is compared to the target sequence(s). 

``upstreamDist``
  The distance upstream of the incorporation site to look for a target
  similar sequence within the criteria specified by ``maxTargetMismatch``.

``downstreamDist``
  The distance downstream of the incorporation site to look / include for a 
  target similar sequence within the criteria specified by 
  ``maxTargetMismatch``.

``pileUpMin``
  An integer value indicating the number of alignments required to overlap
  before being considered a 'pileUp'.

``recoverMultihits``
  While multihit alignments are often difficult to analyze, some information 
  can still be gleamed from the data given reasonable assumptions. Adjusting 
  this parameter to ``TRUE`` will still only focuses on sites that are uniquely 
  mapped, but if a multihit includes a unique site and other locations, 
  contributions are given to the unique site location. Further, reads and their 
  contributions, umitags and fragments, are not double counted but instead 
  evenly distributed to all included unique sites. **Note**, some sequencing 
  artifacts may arrise in "off-target" associated sites. Users should be careful
  to conclude anything from these alignment artifacts. Leaving this option as 
  ``FALSE`` is recommended if the user does not have a target sequence that 
  locates a repetitive sequence. 


Report
""""""

``suppFile``
  Logical (``TRUE`` or ``FALSE``), if the supplemental file provided in 
  ``Supplemental_Info`` should be used in the default report generated at the
  end of processing. If set to ``FALSE``, the ``Supplemental_Info`` parameter
  is not required for processing.

``{tables/figures}``
  Logicals indicating if tables and figures should be generated from the report. 
  Data will be included under the ``reports`` directory in the project run directory. 
  For figures, both PDF and PNG formats will be generated if set to ``TRUE`` at 300 dpi
  while tables will be generated in a comma-separated values (csv) format.

``reportData``
  Logical indicating if a RData object should be saved during the report 
  generation in the ``reports`` directory.

``infoGraphic``
  Logical indicating if an info graphic displaying the genomic distribution of 
  incorporations should be generated at the beginning of the report. While 
  aesthetically pleasing, the graphic gives the report a unique twist and can 
  provide the knowledgeable user with information about the report at the very
  beginning.

``signature``
  Character string included at the beginning of reports to denote the author,
  analyst, laboratory, etc. Make sure you change if you don't want Chris 
  getting credit for your work.




Nuclease Profiles
-----------------

An additional component to the first part of the config file, is the Nuclease
Profiles. The user can specify which nuclease they are using and include
and profile to help identify edit sites. Nuclease can range from Cas9 to Cpf1
or TALEN based nickases. 

**Note:** For TALEN and dual flanking nickases / nucleases, each side will need
to be input as a different target. Specify in `Target_Sequences` the sequence
and `On_Target_Sites` the actual editing site. Make sure you include two 
distinct identifiers for the sequences on-target sites, then specify the 
target treatment as `{target_seq1};{target_seq2}.

Any name can be given in the `Nuclease` section, but that name needs to match
the profile name as well. So if you want to call it "Cas9v2", then just make 
sure you have a profile named "Cas9v2".

Below is some ascii art that indicates the differences between nucleases. 
Additionally, below the art are example profiles for input into the iGUIDE 
software::

  Editing strategies by designer nucleases
  Cas9 :
                   ><   PAM
  ATGCATGCATGCATGCATGCA TGG (sense strand)
  
   TGCATGCATGCATGCATGCA NGG # gRNA
   |||||||||||||||||||| |||
  TACGTACGTACGTACGTACGT ACC (anti-sense strand)
                   ><       # Dominant cutpoint
  
  Cpf1 : Also known as Cas12a (similar nuclease structure for CasX)
                          ><         # Dominant cutpoint
  GTTTG ATGCATGCATGCATGCATGCATGCATGC (sense strand)
    PAM
   TTTV ATGCATGCATGCATGCATGCA        # gRNA, nuclease activity leave overhang
   |||| |||||||||||||||||||||
  CTAAC TACGTACGTACGTACGTACGTACGTACG (anti-sense strand)
                              ><     # Dominant cutpoint
  
  TALEN : Protin-DNA binding domain fused with FokI nickase
  ATATATATATATATATATAT GCATGCATGCATGCAT GCGCGCGCGCGCGCGCGCGC (sense strand)
  \\\\\\\\\\\\\\\\\\\\
                      |------->
                               <-------|
                                        \\\\\\\\\\\\\\\\\\\\
  TATATATATATATATATATA CGTACGTACGTACGTA CGCGCGCGCGCGCGCGCGCG (anti-sense strand)
  # Proteins bind flanking the cleavage site and cut in the "insert" sequence.
  
  CasCLOVER : Clo051 or another nickases with CRISPR-based binding domains
  ATCCT ATGCATGCATGCATGCATGC TTAACCGGTTAACCGG TACGTACGTACGTACGTACG CGGTC
    ||| ||||||||||||||||||||                              (sense strand)
    PAM    Target Sequence  \------->
                                     <-------\   Target Sequence   PAM
  (anti-sense strand)                         |||||||||||||||||||| |||
  TAGGA TACGTACGTACGTACGTACG AATTGGCCAATTGGCC ATGCATGCATGCATGCATGC GCCAG


Below are the example profiles::

  Nuclease_Profiles :
  Cas9 :
    PAM : "NGG"
    PAM_Loc : "3p"
    PAM_Tol : 1
    Cut_Offset : -4
    Insert_size : FALSE

  Cpf1 :
    PAM : "TTTV"
    PAM_Loc : "5p"
    PAM_Tol : 1
    Cut_Offset : 26     #(Anywhere between 23 and 28)
    Insert_size : FALSE

  CasX :
    PAM : "TTCN"
    PAM_Loc : "5p"
    PAM_Tol : 1
    Cut_Offset : 22     #(Anywhere between 16 and 29)
    Insert_size : FALSE

  TALEN :
    PAM : FALSE
    PAM_Loc : FALSE
    PAM_Tol : 0
    Cut_Offset : Mid_insert
    Insert_size : "15:21"

  CasCLOVER :
    PAM : "NGG"
    PAM_Loc : "3p"
    PAM_Tol : 1
    Cut_Offset : Mid_insert
    Insert_size : "10:30"


Profile parameters
""""""""""""""""""

``PAM``
  protospacer adjacent motif - should be specified here and can contain 
  ambiguous nucleotides. 
  
``PAM_Loc`` 
  indicates the location of the PAM with respect to the pattern, either '5p', 
  '3p' or FALSE.
  
``PAM_Tol`` 
  indicates the tolerance for mismatches in the PAM sequence (ignorned if PAM 
  is FALSE). 
  
``Cut_Offset`` 
  indicates the offset from the 5' nucleotide of the PAM sequence where the 
  nuclease creates a double strand break, unless PAM is FALSE, then the 5' 
  position of the target sequence (also accepts "mid_insert" to specify middle 
  of region between paired alignments).
  
``Insert_size`` 
  is used if target sequences are expected to flank each other for editing, 
  such as with TALENs, and indicates the expected size of the insert. To input 
  a range, delimit the min and max by a colon, ie. 15:21. All names of 
  nucleases used to treat specimens need to have a profile. Additional profiles
  should be added under the 'Nuclease_Profiles' parameter.














Sample Information Files
************************

Sample information files (or sampleInfo files) contain information that may 
change from specimen to specimen. They need to contain at lease 3 columns of 
information: sample names, barcode 1, and barcode 2 sequences. Additionally, 
other parameters defined in the config file can be defined in the sample 
information file if they change from specimen to specimen. 

Run specific config file will need to point to the sample information files. For 
convienence, a directory can be found at ``iGUIDE/sampleInfo/`` for depositing 
these files.

SampleInfo files need to have a specific naming format that follows 
'{RunName}.sampleinfo.csv'.

An appropriate format for the sample names is "{specimen}-{rep}" where 
'specimen' is an alpha-numeric designator for the specimen and 'rep' is a 
numeric identifier for technical or biological replicates, separated by a dash 
(``-``). Replicates will be pooled during the final analysis, so if you want 
them to be separate in the report, make sure you give each specimen a different 
identifier. 

For example, iGSP0002-1 and iGSP0002-2, will be pooled together for 
the report and analysis, but iGSP0002-1 and iGSP0003-1 will not. These names 
will be used in naming files, so do not include any special characters that will
confuse file managment. Try to stick to common delimiters, such as ``-`` and ``_``.
Using a dot, ``.``, as a delimiter is not currently supported. 

A good practice is to put specimen identifiers at the beginning, replicate 
identifiers at the end following a "-", and anything else descriptive in the 
middle. For example, iGSP0002-neg-1, can specify the orientation the sample was 
processed with.








Supplemental Information Files
******************************

Supplemental information files (or supp files) contain information that may 
change from specimen to specimen. They have only one required column, 
"Specimen", but subsequence columns will be used to define conditions. Let's use
the below supp file as an example.::

  # Supplemental csv file example, padding included for visualization
  Specimen, Nuclease, gRNA
  iGXA,     Cas9,     TRAC
  iGXB,     Cas9,     TRAC
  iGXC,     Cas9,     B2M
  iGXD,     Cas9,     B2M
  iGXE,     Mock,     Mock
  iGXF,     Mock,     Mock
  
This type of setup would indicate that there are 6 specimens to be analyzed 
(iGXA - iGXF). Each of these would correlate with their sampleName'd replicates,
so for iGXA, all samples with the format iGXA-{number} or iGXA-{info}-{number}
would be pooled into the iGXA specimen.

Additionally, there are three conditions, defined by the distinct data excluding
information in the "Specimen" column. So in this case, the conditions are 
"Cas9-TRAC", "Cas9-B2M", and "Mock-Mock". Within the report format, there are 
several analyses that are conditionally based rather than specimen based. This 
adds to the flexibility and utility of the reporting functions supplied with 
iGUIDE. 

If the user would rather ever specimen analyzed independently and reported in 
that manner, then they can either run a report without a supp file or in a supp
file include a column that distinguishes each specimen from each other.

Column names and formating are transferred directly into the report. 
Additionally, this files sets the order presented in the report. If "iGXC"
comes before "iGXB" in the supp file, the it will be orderd as so throughout the
report. Conditions, as well, follow this format. As presented above, the report
will order the conditions in the following order "Cas9-TRAC", "Cas9-B2M", and 
"Mock-Mock", which is the order of first observation.






Outputs and Reports
*******************


Overview of different reports, runstats, summary and iguide report














