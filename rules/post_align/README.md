# Post Alignment Rules
After alignment of consolidated sequences, alignment data needs to be assigned to each read.

For blat-based alignment, alignments from R1 will need to be coupled to to R2 alignments since both R1 and R2 sequences are individually aligned with blat. 

For bwa-based alignments, alignments are expanded back to the read level.

Both methods are quality filtered by parameters given in the configuration files.
