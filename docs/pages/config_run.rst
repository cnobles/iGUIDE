.. _configinfo:

.. contents::
   :depth: 4


Config - Run Specific Information
=================================

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
  
  