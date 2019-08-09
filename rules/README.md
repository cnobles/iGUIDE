# Rules
This directory contains snakemake rule files. Rules are used by SnakeMake to construct a processing pipeline (based on a directed-acyclic graph). 

Rules wil behave differently given the input configurations, adjustable in the configuration, manifest, and sample info files.

## Workflow Rules
The rules within this directory are related to support and architechture.

## Demultiplexing Rules
Demultiplexing is required to identify which reads are associated with which multiplexed sample.

This can be accomplished through single- or dual-indexing. iGUIDE and the original GUIDE-seq method use dual-indexing.

## Sequence Trimming Rules
When DNA is sequenced, it typically contains target or synthetic sequences. These sequences can disrupt the alignment to the reference genome, giving false data. 

Reads are trimmed for leading and overreading sequences, given the data provided in sample manifests or sample info files.

## Sequence Filtering Rules
After nucleotide sequences have been trimmed and selected for based on contained sequence information individually for paired-end reads, the pairs of reads should be compared to filter for only reads which contain both an R1 and R2 sequence read. 

This directory contains the rules for filtering R1 and R2 sequence files.

## Consolidate Rules
Consolidation is the process of identifying all unique sequences and producing a key which tells the user how to expand the file again. 

This is useful for sequence processing, such as alignment, as it reduces the total number of processes, typically leading to an increase in performance.

## Alignment Rules
These rules direct how the sequences will be aligned. BLAT outputs to 'psl' format while BWA and many other aligners output to 'sam/bam' outputs. Other aligners that output to 'sam/bam' formats could be incorporated into the workflow, similar to BWA. Please contact the maintainer(s) if you would like to use a different aligner.

Supported sequence aligners:

* BLAT / blat (v35) - BLAST-like Alignment Tool
* BWA / bwa (v0.7.17) - Burrows-Wheeler Aligner

## Post-Alignment Rules
After alignment of consolidated sequences, alignment data needs to be assigned to each read.

For blat-based alignment, alignments from R1 will need to be coupled to to R2 alignments since both R1 and R2 sequences are individually aligned with blat. 

For bwa-based alignments, alignments are expanded back to the read level.

Both methods are quality filtered by parameters given in the configuration files.

## Processing Rules
Once genomic locations have been identified, the data is processed to identify incorporation sites associated with CRISPR treatment. 

Specific parameters and methods for identifing incorporation sites can be adjusted in the configuration files.

