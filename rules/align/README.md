# Alignment Rules
These rules direct how the consolidated sequences will be aligned. 

Options for alignment:

* blat (v35) - BLAST-like Alignment Tool
* bwa (v??) - Burrows-Wheeler Aligner

## Pros / Cons :
### BLAT
Pros:
* all possible alignments for a sequence above a certain threshold (unique alignment and multihits).
* prefered method for total inclusion and discovery
Cons:
* requires reads (R1 and R2) to be aligned independently
* requires more computational resources and can take much longer

### BWA
Pros:
* fast and resource efficient
* typically gives best match for paired alignment
Cons:
* Minimal success in detecting multiple alignments
* Lower inclusion may lead to missing potential alignments in certain regions of the genome.
