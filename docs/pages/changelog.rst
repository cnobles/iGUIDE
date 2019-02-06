.. _changelog:

.. contents::
   :depth: 2

ChangeLog 
========================

**v0.9.5 (February XXth, 2019)**

* Updated demultiplexing to be more efficient and better HPC compatible.
* Added RefSeq Extended* reference gene set
  + Includeds curated, predicted, and other RefSeq sets

**v0.9.4 (January 30th, 2019)**

* Updated 'report' utility and formating
  + custom templates now accepted
  + included as subcommand, check with 'iguide report -h'
  + pdf and html options report 'nicely' even when printed from either
* Updated build to v0.9.2 to support new formating in report
* Builds are constructed from spec files rather than yaml requirements
* Included the 'clean' subcommand to reduce size of processed projects
  + after cleaning a project, only terminal data files will remain

**v0.9.3 (January 11th, 2019)**

* Added 'list_samples' subcommand to list samples within a project.
* Caught a few bugs and worked them out for smoother processing and reports.

**v0.9.2 (January 7th, 2019)**

* Modified test dataset to run tests quicker and implemented CirclCI checking.

**v0.9.1 (January 6th, 2019)**

* Fixed problematic install for first time conda installers.

**v0.9.0 (January 4th, 2019)**

* Initial release.
* Supports setup and analysis of GUIDE-seq and iGUIDE experiments.
* Documentation on [ReadTheDocs.io](https://iguide.readthedocs.io/en/latest/index.html).
