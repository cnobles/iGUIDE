#' For those reviewing the code below, the following is a small style guide 
#' outlining the various formats for the code. 
#' 
#' Names with "_": objects, inlucding data.frames, GRanges, vectors, ...
#' Names in caMel format: functions or components of objects (i.e. columns 
#' within a data.frame).
#' Names with ".": arguments / options for functions

#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, scipen = 99, width = 999)

code_dir <- dirname(sub(
  pattern = "--file=", 
  replacement = "", 
  x = grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
))

desc <- yaml::yaml.load_file(
  file.path(code_dir, "descriptions/couple.desc.yml")
)

# Set up and gather command line arguments
parser <- argparse::ArgumentParser(
  description = desc$program_short_description,
  usage = "Rscript couple.R <anchorPSL> <adriftPSL> [-h/--help, -v/--version] [optional args]"
)

parser$add_argument(
  "anchorPSL", nargs = 1, type = "character", help = desc$anchorPSL
)

parser$add_argument(
  "adriftPSL", nargs = 1, type = "character", help = desc$adriftPSL
)

parser$add_argument(
  "-k", "--keys", nargs = "*", type = "character", help = desc$keys
)

parser$add_argument(
  "-o", "--uniqOutput", nargs = 1, type = "character", help = desc$uniqOutput
)

parser$add_argument(
  "--condSites", nargs = 1, type = "character", help = desc$condSites
)

parser$add_argument(
  "--chimeras", nargs = 1, type = "character", help = desc$chimeras
)

parser$add_argument(
  "--multihits", nargs = 1, type = "character", help = desc$multihits
)

parser$add_argument(
  "--stat", nargs = 1, type = "character", default = FALSE, help = desc$stat
)

parser$add_argument(
  "-g", "--refGenome", nargs = 1, type = "character", default = "hg38",
  help = desc$refGenome
)

parser$add_argument(
  "--maxAlignStart", nargs = 1, type = "integer", default = 5L,
  help = desc$maxAlignStart
)

parser$add_argument(
  "--minPercentIdentity", nargs = 1, type = "integer", default = 95L,
  help = desc$minPercentIdentity
)

parser$add_argument(
  "--minTempLength", nargs = 1, type = "integer", default = 30L,
  help = desc$minTempLength
)

parser$add_argument(
  "--maxTempLength", nargs = 1, type = "integer", default = 2500L,
  help = desc$maxTempLength
)

parser$add_argument(
  "--keepAltChr", action = "store_true", help = desc$keepAltChr
)

parser$add_argument(
  "--readNamePattern", nargs = 1, type = "character", 
  default = "[\\w\\:\\-\\+]+", help = desc$readNamePattern
)

parser$add_argument(
  "--saveImage", nargs = 1, type = "character", help = desc$saveImage
)


args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

# Argument Conditionals
if( is.null(args$anchorPSL) | is.null(args$adriftPSL) ){
  stop("Anchor and adrift PSL files not found. Please provide.")
}

if( is.null(args$uniqOutput) ){
  stop("Please provide an output file name.")
}

# Print Inputs to terminal
input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(seq_along(args), function(i){
    paste(args[[i]], collapse = ", ")
  })
)

input_table <- input_table[
  match(
    c("anchorPSL :", "adriftPSL :", "keys :", "uniqOutput :", "condSites :", 
      "chimeras :", "multihits :", "stat :", "refGenome :", 
      "maxAlignStart :", "minPercentIdentity :", "minTempLength :", 
      "maxTempLength :", "readNamePattern :"),
    input_table$Variables
  ),
]

cat("Coupler Inputs:")
print(
  data.frame(input_table),
  right = FALSE, 
  row.names = FALSE
)

# Load supporting scripts
source(file.path(code_dir, "supporting_scripts", "printHead.R"))

source(file.path(code_dir, "supporting_scripts", "readKeyFile.R"))

source(file.path(code_dir, "supporting_scripts", "readPSL.R"))

source(file.path(code_dir, "supporting_scripts", "qualityFilter.R"))

source(file.path(code_dir, "supporting_scripts", "processBLATData.R"))

source(file.path(code_dir, "supporting_scripts", "condenseSites.R"))

source(file.path(code_dir, "supporting_scripts", "writeOutputFile.R"))

if( !all(
  c("printHead", "readKeyFile", "readPSL", "qualityFilter", 
    "processBLATData", "condenseSites", "writeOutputFile") %in% ls())
){
  stop(
    "Cannot load supporting scripts. You may need to clone from github again."
  )
}

# Load reference genome
if( grepl(".fa", args$refGenome) ){
  
  if( !file.exists(args$refGenome) ){
    stop("Specified reference genome file not found.")
  }
  
  ref_file_type <- ifelse(grepl(".fastq", args$refGenome), "fastq", "fasta")
  
  ref_genome <- Biostrings::readDNAStringSet(
    args$refGenome, format = ref_file_type
  )
  
}else{
  
  genome <- grep(
    args$refGenome, 
    unique(BSgenome::installed.genomes()), 
    value = TRUE
  )
  
  if( length(genome) == 0 ){
    
    cat("\nInstalled genomes include:\n")
    print(paste(unique(BSgenome::installed.genomes()), collapse = "\n"))
    stop("\nSelected reference '", args$refGenome, "'genome not in list.")
    
  }else if( length(genome) > 1 ){
    
    cat("\nInstalled genomes include:\n")
    print(paste(unique(BSgenome::installed.genomes(), collapse = "\n")))
    stop(paste(
        "\nPlease be more specific about reference genome.", 
        "Multiple matches to input."
    ))
    
  }
  
  suppressMessages(library(genome, character.only = TRUE))
  ref_genome <- get(genome)
  
}

## Set up stat object
if( args$stat != FALSE ){
  
  sampleName <- unlist(strsplit(args$uniqOutput, "/"))
  
  sampleName <- unlist(
    strsplit(sampleName[length(sampleName)], ".", fixed = TRUE)
  )[1]
  
  stat <- data.frame(
    sampleName = vector("character"),
    metric = vector("character"),
    count = vector("character")
  )
  
}

## Load and process alignment data ##
# Create single key file if one for each alignment file.
if( length(args$keys) > 1 ){
  
  anchor_key_type <- stringr::str_extract(args$keys[1], "[\\w]+$")
  
  if( !anchor_key_type %in% c("csv", "tsv", "rds", "RData") ){
    stop(
      "Output key file type not supported. Please use csv, tsv, rds, or RData."
    )
  }
  
  anchor_keys <- readKeyFile(args$keys[1], format = anchor_key_type)
  
  adrift_key_type <- stringr::str_extract(args$key[2], "[\\w]+$")
  
  if( !adrift_key_type %in% c("csv", "tsv", "rds", "RData") ){
    stop(
      "Output key file type not supported. Please use csv, tsv, rds, or RData."
    )
  }
  
  adrift_keys <- readKeyFile(args$keys[2], format = adrift_key_type)
  
  stopifnot(all(c("readNames", "seqID") %in% names(anchor_keys)))
  stopifnot(all(c("readNames", "seqID") %in% names(adrift_keys)))
  
  # Check input for data, if none, write files and exit
  if( nrow(anchor_keys) == 0 | nrow(adrift_keys) == 0 ){
    
    message("No sequences identified in at least one key file.")
    writeNullOutput(args)
    q()
    
  }
  
  # Verify readNames are in the same format.
  anchor_keys$readNames <- stringr::str_extract(
    anchor_keys$readNames, args$readNamePattern
  )
  
  adrift_keys$readNames <- stringr::str_extract(
    adrift_keys$readNames, args$readNamePattern
  )
  
  # Only interested in reads in common between the two.
  common_names <- intersect(anchor_keys$readNames, adrift_keys$readNames)

  # Check intersection is not 0
  if( length(common_names) == 0 | is.null(common_names) ){
    
    message("No sequences in common between key files.")
    writeNullOutput(args)
    q()
    
  }
  
  # Filter names in key files.
  anchor_keys <- anchor_keys[anchor_keys$readNames %in% common_names,]
  adrift_keys <- adrift_keys[adrift_keys$readNames %in% common_names,]
  
  # Create a common key
  adrift_keys <- adrift_keys[
    match(anchor_keys$readNames, adrift_keys$readNames),
  ]
  
  keys <- data.frame(
    "readNames" = anchor_keys$readNames,
    "anchorSeqID" = factor(anchor_keys$seqID),
    "adriftSeqID" = factor(adrift_keys$seqID)
  )
  
  keys$anchorKey <- as.integer(keys$anchorSeqID)
  keys$adriftKey <- as.integer(keys$adriftSeqID)
  keys$readPairKey <- paste0(keys$anchorKey, ":", keys$adriftKey)
  
  # Print beginning of keys
  printHead(
    keys, 
    title = "Beginning of Key for relating reads to sequences.",
    caption = paste0(
      "\tReads: ", length(unique(keys$readNames)), 
      "\n\tUnique Pairings: ", length(unique(keys$readPairKey))
    )
  )
  
}else if( length(args$keys) == 1 ){
  
  key_type <- str_extract(args$keys, "[\\w]+$")
  
  if( !keys_type %in% c("csv", "tsv", "rds", "RData") ){
    stop(
      "Output key file type not supported. Please use csv, tsv, rds, or RData."
    )
  }
  
  keys <- readKeyFile(args$keys, format = key_type)
  stopifnot(all(c("readNames", "anchorSeqID", "adriftSeqID") %in% names(keys)))
  
  if( nrow(keys) == 0 ){
    
    message("No sequences identified in key file.")
    writeNullOutput(args)
    q()
    
  }
  
  keys$anchorSeqID <- factor(keys$anchorSeqID)
  keys$adriftSeqID <- factor(keys$adriftSeqID)
  keys$anchorKey <- as.integer(keys$anchorSeqID)
  keys$adriftKey <- as.integer(keys$adriftSeqID)
  keys$readPairKey <- paste0(keys$anchorKey, ":", keys$adriftKey)
  
  # Print beginning of keys
  printHead(
    keys, 
    title = "Beginning of Key for relating reads to sequences.",
    caption = paste0(
      "\tReads: ", length(unique(keys$readNames)), 
      "\n\tUnique Pairings: ", length(unique(keys$readPairKey))
    )
  )
  
}else if( length(args$keys) > 2 ){
  
  stop("Cannot have more key files than sequence alignment files.")
  
}

# Load psl files and filter reads based on inputs
anchor_hits <- readPSL(args$anchorPSL)
adrift_hits <- readPSL(args$adriftPSL)

# Remove alignments to alternate chromosomes
# This helps in identifying unique locations instead of alignments that appear
# for both the standard and alternate chromosomes, becoming a multihit.
if( !args$keepAltChr ){
  
  anchor_hits <- anchor_hits[
    !stringr::str_detect(anchor_hits$tName, stringr::fixed("_")),
  ]
  
  adrift_hits <- adrift_hits[
    !stringr::str_detect(adrift_hits$tName, stringr::fixed("_")),
  ]
  
}


# Create base key if no key was supplied
if( is.null(args$keys) ){
  
  anchor_SeqID <- stringr::str_extract(
    unique(anchor_hits$qName), args$readNamePattern
  )
  
  adrift_SeqID <- stringr::str_extract(
    unique(adrift_hits$qName), args$readNamePattern
  )
  
  intersect_SeqID <- intersect(anchor_SeqID, adrift_SeqID)
  
  keys <- data.frame(
    readNames = intersect_SeqID,
    anchorSeqID = factor(intersect_SeqID),
    adriftSeqID = factor(intersect_SeqID)
  )
  
  keys$anchorKey <- as.integer(keys$anchorSeqID)
  keys$adriftKey <- as.integer(keys$adriftSeqID)
  keys$readPairKey <- paste0(keys$anchorKey, ":", keys$adriftKey)

  # Print beginning of keys
  printHead(
    keys, 
    title = "Beginning of Key for relating reads to sequences.",
    caption = paste0(
      "\tReads: ", length(unique(keys$readNames)), 
      "\n\tUnique Pairings: ", length(unique(keys$readPairKey))
    )
  )
  
}


# Print out basic alignment info.
cat(sprintf(
  "Anchor Alignments: %1$s from %2$s sequences", 
  nrow(anchor_hits),
  length(unique(anchor_hits$qName))
))

cat(sprintf(
  "\nAdrift Alignments: %1$s from %2$s sequences\n\n", 
  nrow(adrift_hits),
  length(unique(adrift_hits$qName))
))

# Stop if there are no alignments to couple.
if( nrow(anchor_hits) == 0 | nrow(adrift_hits) == 0 ){
  
  message("No sequences aligned for at least one of the sequence pairs.")
  writeNullOutput(args)
  q()
  
}

# Remove alignments that do not appear in the keys (single reads filtered out)
anchor_hits <- anchor_hits[anchor_hits$qName %in% levels(keys$anchorSeqID),]
adrift_hits <- adrift_hits[adrift_hits$qName %in% levels(keys$adriftSeqID),]

# Quality filter and convert alignments from data.frame to GRanges
anchor_hits <- qualityFilter(
  alignments = anchor_hits, 
  q.start.max = args$maxAlignStart, 
  global.identity.min = args$minPercentIdentity
)                     

if( nrow(anchor_hits) == 0 ){
  
  message("No alignments remaining after quality filtering anchor reads.")
  writeNullOutput(args)
  q()
  
}

anchor_hits <- processBLATData(
  algns = anchor_hits, from = "anchor", ref.genome = ref_genome
)

anchor_hits$anchorKey <- match(anchor_hits$qName, levels(keys$anchorSeqID))

adrift_hits <- qualityFilter(
  alignments = adrift_hits, 
  q.start.max = args$maxAlignStart, 
  global.identity.min = args$minPercentIdentity
)

if( nrow(adrift_hits) == 0 ){
  
  message("No alignments remaining after quality filtering adrift reads.")
  writeNullOutput(args)
  q()
  
}

adrift_hits <- processBLATData(
  algns = adrift_hits, from = "adrift", ref.genome = ref_genome
)

adrift_hits$adriftKey <- match(adrift_hits$qName, levels(keys$adriftSeqID))

# Info after quality filtering individual alignments.
printHead(
  anchor_hits,
  title = "Head of filtered anchor alignments.",
  caption = sprintf(
    "Alignments: %1$s from %2$s reads", 
    length(anchor_hits), 
    length(unique(anchor_hits$qName))
  )
)

printHead(
  adrift_hits,
  title = "Head of filtered adrift alignments.",
  caption = sprintf(
    "Alignments: %1$s from %2$s reads", 
    length(adrift_hits), 
    length(unique(adrift_hits$qName))
  )
)

# Stop if no alignments passed filtering for individual sequences.
if( length(anchor_hits) == 0 | length(adrift_hits) == 0 ){
  
  message("No alignments remaining after quality filtering for at least one of the sequence pairs.")
  writeNullOutput(args)
  q()
  
}

# All alignments should be either "+" or "-" strand.  
stopifnot(all(strand(anchor_hits) == "+" | strand(anchor_hits) == "-"))
stopifnot(all(strand(adrift_hits) == "+" | strand(adrift_hits) == "-"))

# Identify all combinations of unique anchor and adrift sequences present in the 
# data
unique_key_pairs <- unique(keys[,c("anchorKey", "adriftKey", "readPairKey")])

#' Reduced alignments identify the distinct genomic locations present in the 
#' data for the adrift sequences (breakpoint positions) and anchor sequences 
#' (integration site position). 
#' Levels: Reads --> Unique Sequences --> Alignments --> Unique Genomic Loci
red_anchor_hits <- GenomicRanges::reduce(
  x = GenomicRanges::flank(anchor_hits, -1, start = TRUE), 
  min.gapwidth = 0L, 
  with.revmap = TRUE
)

red_adrift_hits <- GenomicRanges::reduce(
  x = GenomicRanges::flank(adrift_hits, -1, start = TRUE), 
  min.gapwidth = 0L, 
  with.revmap = TRUE
)

#' The following finds all posible combinations of anchor and adrift loci which 
#' meet criteria for pairing. These include: oneEach (each pairing must come 
#' from one anchor and one adrift loci), opposite strands (paired loci should be
#' present on opposite strands), and correct downstream orientation (if an 
#' anchor loci is on the "+" strand, then the start of the anchor loci should be
#' less than the paired adrift, and vice versa for "-" strand).
#' (Inherent check for oneEach with findOverlaps())
pairs <- GenomicRanges::findOverlaps(
  query = red_anchor_hits, 
  subject = red_adrift_hits,
  maxgap = args$maxTempLength,
  ignore.strand = TRUE
)

anchor_loci <- red_anchor_hits[S4Vectors::queryHits(pairs)]
adrift_loci <- red_adrift_hits[S4Vectors::subjectHits(pairs)]

#Stop if no alignments coupled based on criteria.
if( length(pairs) == 0 ){

  message("No alignments coupled based on input criteria.")
  writeNullOutput(args)
  q()
  
}

# Check isDownstream and isOppositeStrand
adrift_loci_starts <- GenomicRanges::start(adrift_loci)
anchor_loci_starts <- GenomicRanges::start(anchor_loci)

adrift_loci_strand <- GenomicRanges::strand(adrift_loci)
anchor_loci_strand <- GenomicRanges::strand(anchor_loci)

keep_loci <- ifelse(
  anchor_loci_strand == "+", 
  as.vector(
    (adrift_loci_starts > anchor_loci_starts) & 
      (adrift_loci_strand != anchor_loci_strand)
  ), 
  as.vector(
    (adrift_loci_starts < anchor_loci_starts) & 
      (adrift_loci_strand != anchor_loci_strand)
  )
)

keep_loci <- as.vector(
  (keep_loci & anchor_loci_strand != "*") & (adrift_loci_strand != "*")
)

anchor_loci <- anchor_loci[keep_loci]
adrift_loci <- adrift_loci[keep_loci]

# Stop if no loci were properly paired
if( length(anchor_loci) == 0 | length(adrift_loci) == 0 ){
  
  message("No genomic loci from alignments were properly paired.")  
  writeNullOutput(args)
  q()
  
}

#' Below, the code constructs a genomic loci key which links genomic loci to
#' the various anchor and adrift sequences that were aligned. The technique used
#' below first matches the unique loci back to multiple alignments, then uses 
#' the indices of the unique_key_pairs data.frame (which matches alignments to 
#' unique sequence identifiers) as a GRanges object to match many alignments to 
#' many read identifiers with findOverlaps. For some reason, this method did not
#' work as anticipated with IRanges, and therefore objects were moved to GRanges
#' and GRangesLists.
loci_key <- data.frame(
  "anchorLoci" = S4Vectors::queryHits(pairs)[keep_loci],
  "adriftLoci" = S4Vectors::subjectHits(pairs)[keep_loci]
)

loci_key$lociPairKey <- paste0(loci_key$anchorLoci, ":", loci_key$adriftLoci)

loci_key$anchorKey <- IRanges::IntegerList(
  split(
    anchor_hits$anchorKey[unlist(anchor_loci$revmap)],
    S4Vectors::Rle(
      values = seq_along(anchor_loci), 
      lengths = width(anchor_loci$revmap@partitioning)
    )
  )
)

loci_key$adriftKey <- IRanges::IntegerList(
  split(
    adrift_hits$adriftKey[unlist(adrift_loci$revmap)],
    S4Vectors::Rle(
      values = seq_along(adrift_loci), 
      lengths = width(adrift_loci$revmap@partitioning)
    )
  )
)

anchor_readPair_hits <- GenomicRanges::findOverlaps(
  query = split(
    x = GenomicRanges::GRanges(
      seqnames = "int", 
      ranges = IRanges::IRanges(start = unlist(loci_key$anchorKey), width = 1), 
      strand = "*"
    ), 
    f = S4Vectors::Rle(
      values = seq_len(nrow(loci_key)), 
      lengths = S4Vectors::width(loci_key$anchorKey@partitioning)
    )
  ),
  subject = GenomicRanges::GRanges(
    seqnames = "int",
    ranges = IRanges::IRanges(start = unique_key_pairs$anchorKey, width = 1),
    strand = "*")
)

loci_key$anchorReadPairs <- IRanges::IntegerList(
  split(
    x = S4Vectors::subjectHits(anchor_readPair_hits), 
    f = S4Vectors::queryHits(anchor_readPair_hits)
  )
)

adrift_readPair_hits <- GenomicRanges::findOverlaps(
  query = split(
    x = GenomicRanges::GRanges(
      seqnames = "int", 
      ranges = IRanges::IRanges(start = unlist(loci_key$adriftKey), width = 1), 
      strand = "*"
    ), 
    f = S4Vectors::Rle(
      values = seq_len(nrow(loci_key)), 
      lengths = width(loci_key$adriftKey@partitioning)
    )
  ),
  subject = GenomicRanges::GRanges(
    seqnames = "int",
    ranges = IRanges::IRanges(start = unique_key_pairs$adriftKey, width = 1),
    strand = "*")
)

loci_key$adriftReadPairs <- IRanges::IntegerList(
  split(
    x = S4Vectors::subjectHits(adrift_readPair_hits), 
    f = S4Vectors::queryHits(adrift_readPair_hits)
  )
)

#' Determine intersecting readPairs between the anchor and adrift loci.
anchor_readPair_indices <- GenomicRanges::GRanges(
  seqnames = S4Vectors::queryHits(anchor_readPair_hits), 
  ranges = IRanges::IRanges(
    start = S4Vectors::subjectHits(anchor_readPair_hits), width = 1
  ),
  strand = "*"
)

adrift_readPair_indices <- GenomicRanges::GRanges(
  seqnames = S4Vectors::queryHits(adrift_readPair_hits), 
  ranges = IRanges::IRanges(
    start = S4Vectors::subjectHits(adrift_readPair_hits), width = 1
  ),
  strand = "*"
)

intersecting_readPairs <- GenomicRanges::findOverlaps(
  query = anchor_readPair_indices, 
  subject = adrift_readPair_indices
)

#' Using the range information from the filtered paired alignments, the code
#' constructs a GRanges object from the anchor_loci and adrift_loci. anchor_loci
#' are the integration site positions while the adrift_loci are the various 
#' breakpoints. The strand of the range is set to the same strand as the 
#' anchor_loci since the direction of sequencing from the viral or vector genome
#' is from the U5-host junction found at the 3' end of the integrated element.
paired_loci_key <- loci_key[
  as.integer(unique(GenomicRanges::seqnames(
    anchor_readPair_indices[S4Vectors::queryHits(intersecting_readPairs)]
  ))),
]


paired_loci_key$readPairKeys <- IRanges::CharacterList(
  split(
    x = unique_key_pairs$readPairKey[
      GenomicRanges::start(anchor_readPair_indices[
        S4Vectors::queryHits(intersecting_readPairs)
      ])
    ],
    f = as.integer(GenomicRanges::seqnames(
      anchor_readPair_indices[S4Vectors::queryHits(intersecting_readPairs)]
    ))
  )
)

select_anchor_loci <- anchor_loci[
  as.integer(unique(GenomicRanges::seqnames(
    anchor_readPair_indices[S4Vectors::queryHits(intersecting_readPairs)]
  )))
]

select_adrift_loci <- adrift_loci[
  as.integer(unique(GenomicRanges::seqnames(
    adrift_readPair_indices[S4Vectors::subjectHits(intersecting_readPairs)]
  )))
]

paired_loci <- GenomicRanges::GRanges(
  seqnames = GenomicRanges::seqnames(select_anchor_loci), 
  ranges = IRanges::IRanges(
    start = ifelse(
      GenomicRanges::strand(select_anchor_loci) == "+", 
      GenomicRanges::start(select_anchor_loci), 
      GenomicRanges::start(select_adrift_loci)),
    end = ifelse(
      GenomicRanges::strand(select_anchor_loci) == "+", 
      GenomicRanges::end(select_adrift_loci), 
      GenomicRanges::end(select_anchor_loci))),
  strand = GenomicRanges::strand(select_anchor_loci),
  lociPairKey = paired_loci_key$lociPairKey
)

#' Information on valid paired alignments from all sequences present.
printHead(
  paired_loci,
  title = "Head of valid paired loci present in the data.",
  caption = sprintf("Genomic loci: %s", length(paired_loci))
)

#' Add readPairKeys to the paired loci for expansion
paired_loci$readPairKeys <- paired_loci_key$readPairKeys

#' Stop if there are no paired_loci
if( length(paired_loci) == 0 ){
  
  message("No valid paired genomic loci were found within the data given input criteria.")
  writeNullOutput(args)
  q()
  
}

#' Expand readPairKeys and lociPairKeys to make a single object that maps loci
#' to unique sequences. This is analogous to a sparse matrix, but in a 
#' data.frame object. The keys object is still needed to jump from readPairKey
#' to readName.
read_loci_mat <- data.frame(
  "lociPairKey" = S4Vectors::Rle(
    values = paired_loci$lociPairKey,
    lengths = S4Vectors::width(paired_loci$readPairKeys@partitioning)),
  "readPairKey" = unlist(paired_loci$readPairKeys)
)

#' Templates aligning to single loci are termed unique, while templates
#' aligning to multiple loci are termed multihits.
read_pair_counts <- table(read_loci_mat$readPairKey)
uniq_read_pairs <- names(read_pair_counts[read_pair_counts == 1])
multihit_read_pairs <- names(read_pair_counts[read_pair_counts > 1])

#' Bin reads that would map to different loci on the same read (chimeras)
#' All unique and multihit templates were mapped successfully to 
#' genomic loci, yet some templates were sequenced but did not make it through
#' the selection criteria. These template either do not have alignments to the
#' reference genome (anchor or adrift did not align) or map to two distant 
#' genomic loci. The latter are termed chimeras and are considered to be 
#' artifacts of PCR amplification.
if( !is.null(args$chimeras) ){
  
  failed_reads <- keys[!keys$readPairKey %in% read_loci_mat$readPairKey,]
  
  chimera_reads <- failed_reads[
    failed_reads$anchorKey %in% anchor_hits$anchorKey & 
      failed_reads$adriftKey %in% adrift_hits$adriftKey,
  ]
  
  chimera_alignments <- GenomicRanges::GRangesList()
  
  if( nrow(chimera_reads) > 0 ){
    
    chim_anchor <- anchor_hits[
      anchor_hits$anchorKey %in% chimera_reads$anchorKey,
    ]
    
    chim_anchor <- split(x = chim_anchor, f = chim_anchor$qName)
    chim_anchor <- chim_anchor[chimera_reads$anchorSeqID]
    names(chim_anchor) <- chimera_reads$readNames
    chim_anchor <- unlist(chim_anchor)
    
    chim_adrift <- adrift_hits[
      adrift_hits$adriftKey %in% chimera_reads$adriftKey,
    ]
    
    chim_adrift <- split(x = chim_adrift, f = chim_adrift$qName)
    chim_adrift <- chim_adrift[chimera_reads$adriftSeqID]
    names(chim_adrift) <- chimera_reads$readNames
    chim_adrift <- unlist(chim_adrift)
    
    keepCols <- c(
      "from", "qName", "matches", "repMatches", "misMatches", "qStart", "qEnd", 
      "qSize", "tBaseInsert"
    )
    
    mcols(chim_anchor) <- mcols(chim_anchor)[,keepCols]
    mcols(chim_adrift) <- mcols(chim_adrift)[,keepCols]
    
    chimera_alignments <- c(chim_anchor, chim_adrift)
    chimera_alignments <- split(chimera_alignments, names(chimera_alignments))
    
  }
  
  if( args$stat != FALSE ){
    
    add_stat <- data.frame(
      sampleName = sampleName,
      metric = "chimera.reads",
      count = length(unique(chimera_reads$readNames))
    )
    
    stat <- rbind(stat, add_stat)
    
  }
  
  chimeraData <- list(
    "read_info" = chimera_reads, 
    "alignments" = chimera_alignments
  )
  
  writeOutputFile(chimeraData, file = args$chimeras, format = "rds")
  
}

#' Expand out reads uniquely mapped reads or unique sites
#' Below, the paired_loci object is expanded to create the genomic alignments
#' for each read that mapped to a single genomic loci. This data is then 
#' recorded in two formats. "allSites" is a GRanges object where each row is a
#' single read, while "sites.final" is a condensed form of the data where each
#' row is a unique integration site with the width of the range refering to 
#' the longest template aligned to the reference genome. 
uniq_read_loci_mat <- read_loci_mat[
  read_loci_mat$readPairKey %in% uniq_read_pairs,
]

uniq_templates <- paired_loci[
  match(uniq_read_loci_mat$lociPairKey, paired_loci$lociPairKey)
]

uniq_templates$readPairKeys <- NULL
uniq_templates$readPairKey <- uniq_read_loci_mat$readPairKey

uniq_keys <- keys[keys$readPairKey %in% uniq_read_pairs,]

uniq_reads <- uniq_templates[
  match(uniq_keys$readPairKey, uniq_templates$readPairKey)
]

names(uniq_reads) <- as.character(uniq_keys$readNames)

uniq_reads$sampleName <- stringr::str_extract(
  string = as.character(keys$anchorSeqID[
    match(uniq_reads$readPairKey, keys$readPairKey)
  ]), 
  pattern = "^[\\w-]+"
)

uniq_reads$ID <- names(uniq_reads)

uniq_sites <- uniq_reads
names(uniq_sites) <- NULL
writeOutputFile(uniq_sites, file = args$uniqOutput)

# Print out head of uniq_sites for reference.
printHead(
  uniq_sites,
  title = "Head of uniquely mapped genomic loci.",
  caption = sprintf(
    "Alignments yeilded %1$s unique anchor sites from %2$s properly-paired and aligned reads.",
    length(reduce(flank(uniq_sites, -1, start = TRUE), min.gapwidth = 0L)),
    length(uniq_sites)
  )
)
 
if( args$stat != FALSE ){
  
  add_stat <- data.frame(
    sampleName = sampleName,
    metric = c("unique.reads", "unique.algns", "unique.loci"),
    count = c(
      length(unique(uniq_sites$ID)), 
      length(unique(uniq_sites)),
      length(GenomicRanges::reduce(
        x = GenomicRanges::flank(uniq_sites, width = -1, start = TRUE), 
        min.gapwidth = 0L
      ))
    )
  )
  
  stat <- rbind(stat, add_stat)
  
}

# Generate condensed sites
if( !is.null(args$condSites) ){
  
  cond_sites <- condenseSites(
    uniq_sites, keep.cols = "sampleName", list.bp.counts = TRUE
  )
  
  writeOutputFile(cond_sites, file = args$condSites)
  
  printHead(
    cond_sites,
    title = "Head of unique anchor sites.",
    caption = sprintf(
      "There were %1$s unique anchor sites identified with a total of %2$s unique template lengths and %3$s read counts.",
      length(cond_sites),
      sum(cond_sites$fragLengths),
      sum(cond_sites$counts)
    )
  )
  
}

#' Clean up environment for expansion and clustering of multihits
#rm(uniq_read_loci_mat, uniq_templates, uniq_keys, 
#   uniq_reads, uniq_sites, )
#gc()

#' Group and characterize multihits 
#' Multihits are reads that align to multiple locations in the reference 
#' genome. There are bound to always be a certain proportion of reads aligning
#' to repeated sequence due to the high level degree of repeated DNA elements
#' within genomes. The final object generated, "multihitData", is a list of 
#' three objects. "unclustered_multihits" is a GRanges object where every 
#' alignment for every multihit read is present in rows. 
#' "clustered_multihit_positions" returns all the possible integration site 
#' positions for the multihit. Lastly, "clustered_multihit_lengths" contains the
#' length of the templates mapping to the multihit clusters, used for
#' abundance calculations.
if( !is.null(args$multihits) ){
  
  unclustered_multihits <- GenomicRanges::GRanges()
  clustered_multihit_positions <- GenomicRanges::GRangesList()
  clustered_multihit_lengths <- list()
  
  if( length(multihit_read_pairs) > 0 ){
    
    #' Only consider readPairKeys that aligned to multiple genomic loci
    multi_read_loci_mat <- read_loci_mat[
      read_loci_mat$readPairKey %in% multihit_read_pairs,
    ]
    
    multihit_templates <- paired_loci[
      paired_loci$lociPairKey %in% multi_read_loci_mat$lociPairKey
    ]
    
    multihit_expansion_map <- multihit_templates$readPairKeys
    multihit_templates$readPairKeys <- NULL
    
    multihit_templates <- multihit_templates[S4Vectors::Rle(
      values = seq_along(multihit_templates),
      lengths = S4Vectors::lengths(multihit_expansion_map)
    )]
    
    multihit_templates$readPairKey <- unlist(multihit_expansion_map)
    
    #' As the loci are expanded from the paired_loci object, unique templates 
    #' and readPairKeys are present in the readPairKeys unlisted from the 
    #' paired_loci object.
    multihit_templates <- multihit_templates[
      multihit_templates$readPairKey %in% multi_read_loci_mat$readPairKey
    ]
    
    multihit_keys <- keys[keys$readPairKey %in% multihit_read_pairs,]
    
    multihit_keys$sampleName <- stringr::str_extract(
      string = as.character(multihit_keys$anchorSeqID), pattern = "^[\\w-]+"
    )
    
    multihit_keys$ID <- multihit_keys$readNames
    
    #' Medians are based on all the potential sites for a given read, which will
    #' be identical for all reads associated with a readPairKey.
    multihit_medians <- round(
      median(GenomicRanges::width(split(
        x = multihit_templates, 
        f = multihit_templates$readPairKey
      )))
    )
    
    multihit_keys$medians <- multihit_medians[multihit_keys$readPairKey]
    
    multihits_pos <- GenomicRanges::flank(
      x = multihit_templates, width = -1, start = TRUE
    )
    
    multihits_red <- GenomicRanges::reduce(
      x = multihits_pos, min.gapwidth = 5L, with.revmap = TRUE
    )  #! Should make min.gapwidth a option
    
    revmap <- multihits_red$revmap
    
    axil_nodes <- as.character(S4Vectors::Rle(
      values = multihit_templates$readPairKey[
        S4Vectors::start(revmap@partitioning)
      ], 
      lengths = S4Vectors::width(revmap@partitioning)
    ))
    
    nodes <- multihit_templates$readPairKey[unlist(revmap)]
    edgelist <- unique(matrix( c(axil_nodes, nodes), ncol = 2 ))
    
    multihits_cluster_data <- igraph::clusters(
      igraph::graph.edgelist(el = edgelist, directed = FALSE)
    )
    
    clus_key <- data.frame(
      row.names = unique(as.character(t(edgelist))),
      "clusID" = multihits_cluster_data$membership
    )
    
    multihits_pos$clusID <- clus_key[multihits_pos$readPairKey, "clusID"]
    multihits_pos <- multihits_pos[order(multihits_pos$clusID)]
    clustered_multihit_index <- as.data.frame(mcols(multihits_pos))
    
    multihit_loci_rle <- S4Vectors::Rle(factor(
      x = clustered_multihit_index$lociPairKey, 
      levels = unique(clustered_multihit_index$lociPairKey)
    ))
    
    multihit_loci_intL <- split(
      multihit_loci_rle, clustered_multihit_index$clusID
    )
    
    clustered_multihit_positions <- GenomicRanges::granges(
      x = multihits_pos[
        match(
          x = unlist(S4Vectors::runValue(multihit_loci_intL)), 
          table = clustered_multihit_index$lociPairKey)
      ]
    )
    
    clustered_multihit_positions <- split(
      x = clustered_multihit_positions,
      f = S4Vectors::Rle(
        values = seq_along(multihit_loci_intL), 
        lengths = S4Vectors::width(S4Vectors::runValue(
          multihit_loci_intL
        )@partitioning)
      )
    )
    
    readPairKey_cluster_index <- unique(
      clustered_multihit_index[,c("readPairKey", "clusID")]
    )
    
    multihit_keys$clusID <- readPairKey_cluster_index$clusID[
      match(multihit_keys$readPairKey, readPairKey_cluster_index$readPairKey)
    ]
    
    multihit_keys <- multihit_keys[order(multihit_keys$medians),]
    
    clustered_multihit_lengths <- split(
      x = S4Vectors::Rle(multihit_keys$medians), 
      f = multihit_keys$clusID
    )
    
    #' Expand the multihit_templates object from readPairKey specific to read
    #' specific.
    multihit_keys <- multihit_keys[order(multihit_keys$readPairKey),]
    
    multihit_readPair_read_exp <- IRanges::IntegerList(
      split(x = seq_len(nrow(multihit_keys)), f = multihit_keys$readPairKey)
    )
    
    unclustered_multihits <- multihit_templates
    
    multihit_readPair_read_exp <- multihit_readPair_read_exp[
      as.character(unclustered_multihits$readPairKey)
    ]
    
    unclustered_multihits <- unclustered_multihits[S4Vectors::Rle(
      values = seq_along(unclustered_multihits),
      lengths = S4Vectors::width(multihit_readPair_read_exp@partitioning)
    )]
    
    names(unclustered_multihits) <- multihit_keys$ID[
      unlist(multihit_readPair_read_exp)
    ]
    
    unclustered_multihits$ID <- multihit_keys$ID[
      unlist(multihit_readPair_read_exp)
    ]
    
    unclustered_multihits$sampleName <- multihit_keys$sampleName[
      unlist(multihit_readPair_read_exp)
    ]
    
  }
  
  stopifnot(
    length(clustered_multihit_positions) == length(clustered_multihit_lengths)
  )
  
  multihitData <- list(
    unclustered_multihits, 
    clustered_multihit_positions, 
    clustered_multihit_lengths
  )
  
  names(multihitData) <- c(
    "unclustered_multihits", 
    "clustered_multihit_positions", 
    "clustered_multihit_lengths"
  )
  
  writeOutputFile(multihitData, file = args$multihits, format = "rds")
  
  print(
    data.frame(
      "multihit_reads" = length(unique(names(unclustered_multihits))),
      "multihit_alignments" = length(unique(unclustered_multihits)),
      "multihit_clusters" = length(clustered_multihit_positions),
      "multihit_lengths" = sum(lengths(clustered_multihit_lengths))
    )
  )
  
  if( args$stat != FALSE ){
    
    add_stat <- data.frame(
      sampleName = sampleName,
      metric = c("multihit.reads", "multihit.lengths", "multihit.clusters"),
      count = c(
        length(unique(names(unclustered_multihits))), 
        sum(lengths(clustered_multihit_lengths)), 
        length(clustered_multihit_positions))
    )
    
    stat <- rbind(stat, add_stat)
    
  }
  
}

if( args$stat != FALSE ){

  write.table(
    x = stat, file = args$stat, 
    sep = ",", row.names = FALSE, 
    col.names = FALSE, quote = FALSE
  )
  
}

if( !is.null(args$saveImage) ) save.image(args$saveImage)

q()
