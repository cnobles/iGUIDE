.. iGUIDE documentation master file, created by
   sphinx-quickstart on Fri Nov  2 14:34:12 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to iGUIDE's documentation
==================================

===========
Description
===========

Software pipeline for processing and analyzing double-strand DNA break events. 
These events may be induced, such as by designer nucleases like Cas9, or 
spontaneous, as produced through DNA replication or ionizing radiation. A 
laboratory bench-side protocol accompanies this software pipeline, and can be 
found XXX. 

Below, this documentation gives the reader a overview of the pipeline, including
how to install and process a sample dataset. Processing a sample data set is 
broken into a few parts: 

#. developing a configuration file and sample information
#. setting up a run directory and acquiring the sequence data
#. initializing the pipeline and understanding the output

.. toctree::
   :hidden:
   :maxdepth: 4
   :caption: Contents:

   pages/install.rst
   pages/quickstart.rst
   pages/configinfo.rst
   pages/sampleinfo.rst   
