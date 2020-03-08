.. _changelog:

ChangeLog
=========

**v1.1.0 (March 8th, 2020)**

- Modified how samples designated as Mock are treated during the analysis
- Mock samples can now be indicated by "None" or "Control" as well 
  (case-insensitive)
- Abundance can now be selected as [Read], [UMI], or [Fragment]{default} within 
  config parameters and this selection will identify the abundance method used
  for analysis
- Added support for alternative UMI method (dx.doi.org/10.17504/protocols.io.wikfccw)

**v1.0.2 (February 15th, 2020)**

- Bugfix: UMItags set to FALSE will now process through to completion
- Rebuild: Updated to build version 1.0.1

**v1.0.1 (December 3rd, 2019)**

- Bugfix: Updated Gene set enrichment test within report

**v1.0.0 (August 15th, 2019)**

- Complete support for BLAT and BWA aligners
- Included a binning system to distribute workload into smaller loads
- Implemented a version tracking system into the intermediate data files
  (incorp_sites)
- Updated CLI with "hints" for snakemake processing

**v0.9.9 (August 9th, 2019) - Additional updates**

- Implemented support for BWA aligner
- Added tools (samqc) for working with other SAM/BAM output aligners as well
- Switched iguide support code to iguideSupport R-package and added unit tests
- Fixed bugs related to quoted table inputs (csv/tsv)
- Implemented a method to skip demultiplexing, see documentation for setup
- Resoved a number of issues identified, check GitHub for history!

**v0.9.9 (June 10th, 2019)**

- Revised the iGUIDE Report format to be more informational and clearer
- Revised a bit of the workflow to make reprocessing smoother
- Updated BLAT coupling script to be more memory efficient
- Fixed TravisCI testing!
- Changed stat workflow, now restarting analysis won't initiate a total 
  reprocessing.
- Modified the assimilate + evaluate workflow
- Assimilate now only includes reference genome data, meaning a cleaner
  intermediate file
- Evaluate will now handle ref. gene sets and further analysis
- This increases the modularity and consistancy of the workflow


**v0.9.8 (April 19th, 2019)**

- iGUIDE can now support non-Cas9 nucleases as well!
- Implemented nuclease profiles into configs
- Updated assimilation, evaluation, and reporting scripts
- Added default resources to allow simpler HPC processing
- Included flexible system for identifying on-target sites
- Config can accept a range rather than a single site
- Acceptable notation: chr4:+:397-416 and chr3:\*:397
- Changed build nomenclature from v0.9.3 to b0.9.3, so as not to confuse with
  version
- Added 'summary' subcommand to generate a consise text-based report
- Added short stats-based report to be produced at the end of processing
- Additional bugfixes.

**v0.9.7 (March 6th, 2019)**

- Hotfix to workflow.
- Changed 'setup' subcommand to python script based rather than snakemake.
- Changed file organization.

**v0.9.6 (March 5th, 2019)**

- Introduced process workflow steps: assimilate and evaluate
- Assimilate aligned data and compare with targeting sequences
- Incorp_sites now a core data object that can be combined across runs
- Evaluated data incorporates reference data and statistical models
- A staple data object for reports and can be constructed from multiple runs
- Included new subcommands 'eval' and modified 'report', report from either
  config(s) or eval dataset
- Cleaned up file structure
- Updated documentation in code and docs.
- Implemented accuracy and retention checks with simulation dataset.
- Updated simulation dataset with larger set to test analysis.

**v0.9.5 (February 19th, 2019)**

- Updated demultiplexing to be more efficient and better HPC compatible.
- Added RefSeq Extended reference gene sets
- 'ext' includes curated, predicted, and other RefSeq sets
- 'ext.nomodel' includes only curated and other RefSeq sets
- Incorporated resource allocation for job dependent memory consumption, works
  great with HPC to specify memory requirements
- Streamlined input for report generation by only requiring config(s)

**v0.9.4 (January 30th, 2019)**

- Updated 'report' utility and formating. Custom templates now accepted. 
  Included as subcommand, check with 'iguide report -h'. PDF and HTML options
  report 'nicely' even when printed from either
- Updated build to v0.9.2 to support new formating in report
- Builds are constructed from spec files rather than yaml requirements
- Included the 'clean' subcommand to reduce size of processed projects. After
  cleaning a project, only terminal data files will remain

**v0.9.3 (January 11th, 2019)**

- Added 'list_samples' subcommand to list samples within a project.
- Caught a few bugs and worked them out for smoother processing and reports.

**v0.9.2 (January 7th, 2019)**

- Modified test dataset to run tests quicker and implemented CirclCI checking.

**v0.9.1 (January 6th, 2019)**

- Fixed problematic install for first time conda installers.

**v0.9.0 (January 4th, 2019)**

- Initial release.
- Supports setup and analysis of GUIDE-seq and iGUIDE experiments.
- Documentation on [ReadTheDocs.io](https://iguide.readthedocs.io/en/latest/index.html).
