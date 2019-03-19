.. _configinfo:

.. contents::
   :depth: 4

Config - Setup
==============

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
familiar with the sequence data used in the laboratory bench-side protocol. 
Additionally, they should be aware of the biochemistry related to the enzymes 
and sequences they are using.

The second part (below the divide ``----``) should be filled out by an 
individual familiar with the bioinformatic processing. Explanations of the 
different portions can be found in the following pages.

