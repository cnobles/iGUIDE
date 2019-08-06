.. _sampleinfo:

.. contents::
   :depth: 2


Supplemental Information Files
========================

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
