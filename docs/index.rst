.. iGUIDE documentation master file, created by
   sphinx-quickstart on Fri Nov  2 14:34:12 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to iGUIDE's documentation
==================================

iGUIDE is a software pipeline for processing and analyzing double-strand DNA 
break events. These events may be induced, such as by designer nucleases like 
Cas9, or spontaneous, as produced through DNA replication or ionizing 
radiation. A laboratory bench-side protocol accompanies this software pipeline, 
and can be found (https://doi.org/10.1186/s13059-019-1625-3). 

This documentation gives the reader an overview of the pipeline, including
how to install and process a sample dataset. Processing a sample data set is 
broken into a few parts: 

#. Developing a configuration file and sample information
#. Setting up a run directory and acquiring the sequence data
#. Initializing the pipeline and understanding the output

After processing sequencing run(s), additional evaluation and analysis can be 
performed by supplying supplemental data and combining specimen outputs from
several runs.

To get started, see :ref:`quickstart`!

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Contents:

   quickstart.rst
   usage.rst
   changelog.rst
