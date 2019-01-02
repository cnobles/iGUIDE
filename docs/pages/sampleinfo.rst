.. _sampleinfo:

.. contents::
   :depth: 2

========================
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

