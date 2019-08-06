.. _sampleinfo:

.. contents::
   :depth: 2


Sample Information Files
========================

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

