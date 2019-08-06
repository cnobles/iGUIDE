.. _configinfo:

.. contents::
   :depth: 4

Configs - Processing Information
================================

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

``defaultMB / demultiMB / trimMB / filtMB / consolMB / alignMB / coupleMB / assimilateMB / evaluateMB / reportMB``
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


Reference Alignment
"""""""""""""""""""

``BLATparams``
  A character string to be included with the BLAT call. For options, please see
  the BLAT help options by typing ``blat`` into the commandline after 
  activating ``iguide``.

``BWAparams``
  A character string to be inclued with the BWA call. BWA is not currently 
  supported, so this parameter is currently silent.


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
