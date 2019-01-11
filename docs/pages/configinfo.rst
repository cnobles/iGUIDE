.. _configinfo:

.. contents::
   :depth: 4

Setting up a config file
========================

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

Config File Layout
------------------

Config files are in a ``yaml`` format, but are broken into two parts. The first 
contains run specific information that should be filled out by an individual 
familiar with the sequence data used in the laboratory bench-side protocol. The 
second part (below the divide ``----``) should be filled out by an individual 
familiar with the bioinformatic processing. Explanations of the different 
portions can be found below.

Config - Run Specific Information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Run configuration
"""""""""""""""""

``Run_Name``
  This is the name of the sequencing run, and should only contain alpha-numeric
  characters. Underscores (``_``) and dashes (``-``) are also allowed within the
  run name parameters. Other symbols should not be included, such as a dot 
  (``.``). The run name is further used by the software to link files and 
  directories together, so it will need to be consistent whenever it is used.
  
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
  this file, please see the Supplemental Information page.
  
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

``R1 / R2 / I1 / I2``
  These parameters should be the file names of the sequence files to be 
  analyzed by the iGUIDE software. It is recommened to pass complete sequencing
  files to iGUIDE rather than demultiplexing prior to analysis.

SampleInfo formating
""""""""""""""""""""

``Sample_Name_Column``
  This is the name of the column in the sample information file which contains 
  information about samples. An appropriate format for the sample names is 
  "{specimen}-{rep}" where 'specimen' is an alpha-numeric designator for the 
  specimen and 'rep' is a numeric identifier for technical or biological 
  replicates, separated by a dash (``-``).

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

Guide RNA information
"""""""""""""""""""""

``Guide_RNA_Sequences``
  This parameter specifies the guide RNA sequences, including the PAM sequences.
  An acceptable input format would be ``B2M : "GAGTAGCGCGAGCACAGCTANGG"``, and 
  additional guide RNA sequences can be included, one per line, and each 
  indented at the same level. The input format of ``{gRNA_name} : {gRNA_seq}``
  needs to be maintained for proper function. The 'gRNA_name' in this situation
  will need to match the 'gRNA_name' used in the ``On_Target_Sites`` and 
  ``Treatment`` parameters.

``PAM_Sequence``
  A sequence indicating the pattern acquisition motif (PAM) of the guide RNA 
  sequence(s). Multiple PAM sequences can be separated by a linebreak, similar 
  to ``Guide_RNA_Sequences`` but do not need a name. The sequence provided needs
  to be identical to the end of the ``Guide_RNA_Sequences``.
  
``On_Target_Sites``
  This parameter indicates the specific location for editing by the guide RNAs.
  There should be one line for each on-target site, even if there are more than
  one on-target sites for a given guide RNA. Typically the input format should 
  follow ``{gRNA_name} : "{seqname}:{+/-}:{position}"``, where 'gRNA_name' 
  matches the name of the given guide RNA, and if multiple on-target sites 
  exist, then the names can be expanded using a ``{gRNA_name}'#`` notation. The
  value for each on-target site specifies the location or genomic coordinates of
  nuclease activity. The 'seqname' indicates the chromosome or sequence name, an
  orientation of '+' or '-' is given to the location depending on the editing 
  orientation (in line with positional numbering is '+' and opposite is '-'),
  and the 'position' indicates the nucleotide of editing. For Cas9, the position
  of editing is commonly between the 3rd and 4th nucleotide from the 3' end of
  the targeting sequence (not including the PAM). Being off by a nucleotide or 
  so will not cause any problems.

Specimen treatment
""""""""""""""""""

``Treatment``
  This parameter indicates how samples were treated. If samples were all treated
  differently, then this information can be included in the sample information
  file as ``all : "sampleInfo:{column}"`` where 'column' is the name of the 
  column with the information. If a single sample was treated with more than one
  guide RNA, then delimit multiple guide RNA names by a semicolon (``;``), i.e.
  ``all : "B2M;TRAC;TRBC"``. Additionally, each specimen can be indicated 
  individually on a new line. Only specimen names should be given here is 
  provided individually, not sample identifiers.


Config - Processing Information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

Demultiplexing parameters
"""""""""""""""""""""""""

``barcode{1/2}Length``
  Integer values indicating the number of nucleotides in the barcodes or 
  indexing sequences.

``barcode{1/2}``
  Character values (i.e. ``"I1"``) indicating which reads to find the associated
  indexing information for demultiplexing.

``bc{1/2}Mismatch``
  An integer value indicating the number of tolarated mismatches in the barcode
  sequences for either barcode 1 or 2.

``demultiCores``
  The number of core to be requested during demultiplexing. This can be a 
  memory intensive process and therefore can be limited here by using a smaller
  value than given the the ``iguide run`` command.

Sequence trimming
"""""""""""""""""

``R{1/2}leadMismatch``
  Integer values indicating the number of allowed mismatches in either R1 or R2
  leading sequence trimming. 

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
  
``maxGuideMismatch``
  The maximum number of mismatches between the reference genome and guide RNA 
  sequence allowed for consideration to be a guide RNA matched incorporation 
  site. This is an integer value and is compared to the guide RNA sequence(s). 

``upstreamDist``
  The distance upstream of the incorporation site to look for a guide RNA 
  similar sequence within the criteria specified by ``maxGuideMismatch``.

``downstreamDist``
  The distance downstream of the incorporation site to look / include for a 
  guide RNA similar sequence within the criteria specified by 
  ``maxGuideMismatch``.

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
  ``FALSE`` is recommended if the user does not have a guide RNA that targets a 
  repetitive sequence. 

Report
""""""

``suppFile``
  Logical (``TRUE`` or ``FALSE``), if the supplemental file provided in 
  ``Supplemental_Info`` should be used in the default report generated at the
  end of processing.

``figures``
  Logical indicating if figures should be generated from the report. Figures
  will be included under the ``reports`` directory in the run directory. Both
  PDF and PNG formats will be generated if set to ``TRUE`` at 300 dpi.

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
  getting credit for all the work.
