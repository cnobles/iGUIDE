#' For those reviewing the code below, the following is a small style guide 
#' outlining the various formats for the code. 
#' 
#' Names with "_": objects, inlucding data.frames, GRanges, vectors, ...
#' Names in caMel format: functions or components of objects (i.e. columns 
#' within a data.frame).
#' Names with ".": arguments / options for functions

#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, scipen = 99, width = 120)
suppressMessages(library("magrittr"))

code_dir <- dirname(sub(
  pattern = "--file=", 
  replacement = "", 
  x = grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
))

desc <- yaml::yaml.load_file(
  file.path(code_dir, "descriptions/samqc.desc.yml")
)

# Set up and gather command line arguments
parser <- argparse::ArgumentParser(
  description = desc$program_short_description,
  usage = "Rscript samqc.R <bam> <bai> [-h/--help, -v/--version] [optional args]"
)

parser$add_argument(
  "bam", nargs = 1, type = "character", help = desc$bam
)

parser$add_argument(
  "bai", nargs = 1, type = "character", help = desc$bai
)

parser$add_argument(
  "-o", "--uniqOutput", nargs = 1, type = "character", required = TRUE,
  help = desc$uniqOutput
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
  "--batches", nargs = 1, type = "integer", default = 25L,
  help = paste(
    "A tuning parameter to batch process the alignments, specifies how many", 
    "batches to do. Default: 500."
  )
)

parser$add_argument(
  "--readNamePattern", nargs = 1, type = "character", 
  default = "[\\w\\:\\-\\+]+", help = desc$readNamePattern
)

parser$add_argument(
  "--saveImage", nargs = 1, type = "character", help = desc$saveImage
)


args <- parser$parse_args(commandArgs(trailingOnly = TRUE))


# Print Inputs to terminal
input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(seq_along(args), function(i){
    paste(args[[i]], collapse = ", ")
  })
)

input_table <- input_table[
  match(
    c("bam :", "bai :", "uniqOutput :", "condSites :", "chimeras :", 
      "multihits :", "stat :", "refGenome :", "maxAlignStart :", 
      "minPercentIdentity :", "minTempLength :", "maxTempLength :", 
      "keepAltChr :", "readNamePattern :"
    ),
    input_table$Variables
  ),
  ]

cat("\nSAM QC Inputs:\n")
print(
  data.frame(input_table),
  right = FALSE, 
  row.names = FALSE
)


# Load supporting scripts
source(file.path(code_dir, "supporting_scripts", "printHead.R"))
source(file.path(code_dir, "supporting_scripts", "condenseSites.R"))
source(file.path(code_dir, "supporting_scripts", "writeOutputFile.R"))

if( !all(c("printHead", "condenseSites", "writeOutputFile") %in% ls()) ){
  stop(
    "\n  Cannot load supporting scripts. ",
    "You may need to clone from github again.\n"
  )
}

# Load reference genome
if( grepl(".fa", args$refGenome) ){
  
  if( !file.exists(args$refGenome) ){
    stop("\n  Specified reference genome file not found.\n")
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
    stop("\n  Selected reference '", args$refGenome, "'genome not in list.\n")
    
  }else if( length(genome) > 1 ){
    
    cat("\nInstalled genomes include:\n")
    print(paste(unique(BSgenome::installed.genomes(), collapse = "\n")))
    stop(
      "\n  Please be more specific about reference genome. ", 
      "Multiple matches to input.\n"
    )
    
  }
  
  suppressMessages(library(genome, character.only = TRUE))
  ref_genome <- get(genome)
  
}

## Determine an associated sample name
sampleName <- unlist(strsplit(args$uniqOutput, "/"))

sampleName <- unlist(
  strsplit(sampleName[length(sampleName)], ".", fixed = TRUE)
)[1]


## Set up stat object
if( args$stat != FALSE ){
  
  stat <- data.frame(
    sampleName = vector("character"),
    metric = vector("character"),
    count = vector("character")
  )
  
}

# Additional functions ----
#' Load sorted BAM file into a data.frame
#' @param bam path to sorted BAM file (*.bam).
#' @param bai path to BAM index file (*.bai).
#' @param params character vector indicating the fields to import. Refer to 
#' SAMtools or BWA manual for field names.
#' @param tags character vector indicating the additional tags to import. Again,
#' refer to the SAMtools or BWA manual for tag names.

loadBAM <- function(bam, bai, params, tags, onlyPairMapped = TRUE){

  algn <- unlist(Rsamtools::scanBam(
      file = bam, index = bai,
      param = Rsamtools::ScanBamParam(
        flag = Rsamtools::scanBamFlag(
          isPaired = ifelse(onlyPairMapped, TRUE, NA), 
          isUnmappedQuery = ifelse(onlyPairMapped, FALSE, NA),
          hasUnmappedMate = ifelse(onlyPairMapped, FALSE, NA)
        ), 
        what = params, 
        tag = tags
      )
    ),
    recursive = FALSE
  )
  
  df <- as.data.frame(algn[seq_along(params)])
  
  for(t in tags){
    df[,t] <- algn$tag[[t]]
  }
  
  df
  
}

#' Calculate the global percent identity for alignments from cigar and MD tags
#' @param cigar character vector of cigar strings.
#' @param MD character vector of MDz tags.
#' @description Both input parameters must be the same length of vectors and 
#' indexed accordingly. Function calculates the global percent identity of the 
#' alignment.

calcPctID <- function(cigar, MD){
  # Must have same length to calc pct ID
  stopifnot( length(cigar) == length(MD) )
  
  data.frame("cig" = cigar, "md" = MD, stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      mismatch = rowSums(matrix(
        stringr::str_extract_all(md, "[ATGC]", simplify = TRUE) %in% 
          c("A", "T", "G", "C"), 
        nrow = n()), na.rm = TRUE
      ),
      match = rowSums(matrix(as.numeric(gsub(
          "M", "",stringr::str_extract_all(cig, "[0-9]+M", simplify = TRUE)
        )), 
        nrow = n()), na.rm = TRUE
      ) - mismatch,
      length = rowSums(matrix(as.numeric(gsub(
          "[HSMIDX=]", "", stringr::str_extract_all(
            cig, "[0-9]+[HSMIDX=]", simplify = TRUE
          )
        )),
        nrow = n()), na.rm = TRUE
      ),
      pctID = round(100 * (match / length), digits = 1)
    ) %>%
    .$pctID
  
}

#' Count clipped bases at start or end of aligments from cigar strings
#' @param cigar character string with cigar information
#' @param type character indicating hard ("H", "h"), soft ("S", "s"), or both 
#' ("both", default) clipping to be counted.
#' @param end character indicating the 5-prime ("5p", default) or 3-prime ("3p")
#' end of the alignment.

cntClipped <- function(cigar, type = "both", end = "5p"){
  
  # Format inputs
  type <- tolower(type)
  end <- tolower(end)
  
  # Check inputs
  stopifnot( type %in% c("both", "h", "s") )
  stopifnot( end %in% c("5p", "3p") )
  
  # Assign query
  if( type == "both" ){
    query_pat <- "[0-9]+[HS]"
  }else if( type == "h" ){
    query_pat <- "[0-9]+[H]"
  }else{
    query_pat <- "[0-9]+[S]"
  }
  
  # Assign end
  if( end == "5p" ){
    query_pat <- paste0("^", query_pat)
  }else{
    query_pat <- paste0(query_pat, "$")
  }
  
  # Capture all patterns and return integer of clipped bases
  rowSums(matrix(as.numeric(
        gsub("[HS]", "", stringr::str_extract_all(
          cigar, query_pat, simplify = TRUE
        ))
      ), 
      nrow = length(cigar)
    ), 
    na.rm = TRUE
  )
  
}

#' Process alignment data to valid paired-end alignments representing the input
#' template DNA.
#' @param id character vector indicating grouping of alignments.
#' @param chr character vector of seqnames. If using reference genome, these 
#' will need to match seqnames present in the reference object passed to 
#' `refGen`.
#' @param strand character vector of strand or alignment orientation, must be 
#' either "+" or "-".
#' @param pos numeric or integer vector indicating the "start" of the alignment.
#' @param width numeric or integer vector indicating the width of the alignment.
#' @param type character vector indicating type of alignment 
#' ("anchor" or "adrift").
#' @param maxLen numeric or integer value indicating the minimum distance 
#' between the two alignments that should be considered.
#' @param maxLen numeric or integer value indicating the maximum distance 
#' between the two alignments that should be considered.
#' @param refGen BSgenome object or other object with GenomeInfoDb::seqinfo.
#' This method is currently depreciated for the latter method.

.processAlignments <- function(id, chr, strand, pos, width, type, minLen = 30L,
                              maxLen = 2500L, refGen = NULL){
  
  # Check inputs
  inputs <- list(
    "grp" = id, "chr" = chr, "strand" = strand, 
    "pos" = pos, "width" = width, "type" = type
  )
  
  stopifnot( length(unique(sapply(inputs, length))) == 1 ) # All same length
  
  # Combine into data.frame and build GenomicRanges
  input_df <- as.data.frame(inputs) %>%
    dplyr::mutate(
      grp = as.character(grp),
      start = pos,
      end = pos + width - 1,
      type = as.character(type),
      strand = as.character(strand),
      posid = paste0(type, ":", chr, strand, ifelse(strand == "+", start, end))
    ) %>%
    dplyr::select(grp, chr, strand, start, end, type, posid)

  input_gr <- GenomicRanges::GRanges(
    seqnames = as.character(input_df$chr),
    ranges = IRanges::IRanges(
      start = as.numeric(input_df$start), 
      end = as.numeric(input_df$end)
    ),
    strand = as.character(input_df$strand),
    seqinfo = if(!is.null(refGen)){ GenomeInfoDb::seqinfo(refGen) }else{ NULL },
    grp = as.character(input_df$grp),
    type = as.character(input_df$type),
    posid = as.character(input_df$posid)
  )
  
  # Find overlaps within maxLen for anchors and adrift
  grl <- GenomicRanges::split(
    GenomicRanges::flank(input_gr, width = -1, start = TRUE), 
    input_gr$type
  )
  
  # Flip strand of adrift hits to enforce opposite strand requirement
  GenomicRanges::strand(grl$adrift) <- ifelse(
    GenomicRanges::strand(grl$adrift) == "+", "-", "+"
  )
  
  # Reduce to unique locations to minimize work
  red_list <- lapply(
    grl, GenomicRanges::reduce, min.gapwidth = 0L, with.revmap = TRUE
  )
  
  # ID all anchor-to-adrift alignment pairs
  ovlp_hits <- GenomicRanges::findOverlaps(
    red_list$anchor, red_list$adrift, maxgap = maxLen
  )
  
  # Gather data for each type
  anchor_df <- as.data.frame(red_list$anchor) %>%
    dplyr::mutate(
      seqnames = as.character(seqnames),
      strand = as.character(strand),
      type = "anchor",
      anchorid = seq_len(n()),
      posid = paste0("anchor:", seqnames, strand, start)
    )
  
  adrift_df <- as.data.frame(
      red_list$adrift[S4Vectors::subjectHits(ovlp_hits)]
    ) %>%
    dplyr::mutate(
      seqnames = as.character(seqnames),
      strand = as.character(strand),
      type = "adrift",
      anchorid = S4Vectors::queryHits(ovlp_hits),
      posid = paste0(
        "adrift:", seqnames, ifelse(strand == "+", "-", "+"), start
      )
    )
  
  # Combine hits to form valid paired-end alignments
  combo_df <- dplyr::bind_rows(anchor_df, adrift_df) %>%
    dplyr::group_by(anchorid) %>%
    dplyr::mutate(
      anchor.dist = start[type == "anchor"] - start,
      anchor.upstream = ifelse(
        type == "adrift", 
        ifelse(strand == "+", anchor.dist < 0, anchor.dist > 0),
        TRUE
      ),
      right.size = ifelse(
        type == "adrift",
        abs(anchor.dist) >= minLen & abs(anchor.dist) <= maxLen,
        TRUE
      )
    ) %>%
    dplyr::filter(anchor.upstream & right.size) %>%
    dplyr::ungroup()
  
  cond_df <- combo_df %>%
    dplyr::group_by(anchorid) %>%
    dplyr::mutate(
      anchor.posid = posid[type == "anchor"],
      adrift.posid = posid,
      start = ifelse(strand == "+", start - abs(anchor.dist), start),
      end = ifelse(strand == "+", end, end + abs(anchor.dist))
    ) %>%
    dplyr::filter(type == "adrift") %>%
    dplyr::ungroup()
  
  adrift_revmap <- cond_df$revmap
  
  cond_df[rep(seq_len(nrow(cond_df)), lengths(adrift_revmap)),] %>%
    dplyr::mutate(
      grp = grl$adrift$grp[BiocGenerics::unlist(adrift_revmap)]
    ) %>%
    dplyr::filter(
      paste0(grp, ":", anchor.posid) %in% 
        paste0(input_df$grp, ":", input_df$posid)
    ) %>%
    dplyr::mutate(grp = factor(grp, levels = unique(id))) %>%
    dplyr::arrange(grp) %>%
    dplyr::mutate(grp = as.character(grp)) %>%
    dplyr::select("id" = grp, "chr" = seqnames, strand, start, end)

}

#' Process alignment data to valid paired-end alignments representing the input
#' template DNA.
#' @param id character vector indicating grouping of alignments.
#' @param chr character vector of seqnames. If using reference genome, these 
#' will need to match seqnames present in the reference object passed to 
#' `refGen`.
#' @param strand character vector of strand or alignment orientation, must be 
#' either "+" or "-".
#' @param pos numeric or integer vector indicating the "start" of the alignment.
#' @param width numeric or integer vector indicating the width of the alignment.
#' @param type character vector indicating type of alignment 
#' ("anchor" or "adrift").
#' @param maxLen numeric or integer value indicating the minimum distance 
#' between the two alignments that should be considered.
#' @param maxLen numeric or integer value indicating the maximum distance 
#' between the two alignments that should be considered.
#' @param refGen BSgenome object or other object with GenomeInfoDb::seqinfo.
#' @param batches integer indicating the number of batches to serialize the 
#' data processing with. The number of reads analyzed within a batch will be
#' the number of unique `id`'s divided by the `batches`.

processAlignments <- function(id, chr, strand, pos, width, type, minLen = 30L,
                              maxLen = 2500L, refGen = NULL, batches = 25L){
  
  # Check inputs
  inputs <- list(
    "grp" = id, "chr" = chr, "strand" = strand, 
    "pos" = pos, "width" = width, "type" = type
  )
  
  stopifnot( length(unique(sapply(inputs, length))) == 1 ) # All same length
  
  # Combine into data.frame and build GenomicRanges
  input_df <- as.data.frame(inputs) %>%
    dplyr::mutate(
      grp = as.character(grp),
      start = pos,
      end = pos + width - 1,
      type = as.character(type),
      strand = as.character(strand),
      pos = ifelse(strand == "+", start, end)
    ) %>%
    dplyr::select(grp, chr, strand, pos, type)
  
  idx_list <- IRanges::IntegerList(split(seq_len(nrow(input_df)), input_df$grp))
  
  anchor_idx_list <- idx_list[
    IRanges::LogicalList(split(input_df$type == "anchor", input_df$grp))
  ]
  
  adrift_idx_list <- idx_list[
    IRanges::LogicalList(split(input_df$type == "adrift", input_df$grp))
  ]
  
  batch_list <- split(
    seq_along(idx_list), 
    ceiling(seq_along(idx_list) / (length(idx_list) / batches))
  )

  dplyr::bind_rows(lapply(seq_along(batch_list), function(i){
    
    print(i)
    idxs <- batch_list[[i]]
    
    # Identify which reads to analyze
    x <- names(idx_list)[idxs]
    
    # Pull in all anchors associated with reads
    anchor_aligns <- input_df[unlist(anchor_idx_list[x]),]
    
    # Pull in all adrift alignments associated with reads
    adrift_aligns <- input_df[unlist(adrift_idx_list[x]),] %>%
      dplyr::select(grp, "chr.d" = chr, "strand.d" = strand, "pos.d" = pos)
    
    anc_idx <- IRanges::IntegerList(
      split(seq_len(nrow(anchor_aligns)), anchor_aligns$grp)
    )
    
    adr_idx <- IRanges::IntegerList(
      split(seq_len(nrow(adrift_aligns)), adrift_aligns$grp)
    )
    
    exp_anc_idxs <- unlist(lapply(
      seq_along(anc_idx), 
      function(i) rep(anc_idx[[i]], each = length(adr_idx[[i]]))
    ))
        
    adrift_aligns[
        unlist(unname(adr_idx[rep(names(anc_idx), lengths(anc_idx))])),
      ] %>%
      dplyr::mutate(
        chr.n = anchor_aligns$chr[exp_anc_idxs],
        strand.n = anchor_aligns$strand[exp_anc_idxs],
        pos.n = anchor_aligns$pos[exp_anc_idxs]
      ) %>%
      dplyr::filter(
        # Filter for opposite strands
        strand.n != strand.d,
        # Filter for correct size window
        ifelse(strand.n == "+", pos.d - pos.n, pos.n - pos.d) >= minLen,
        ifelse(strand.n == "+", pos.d - pos.n, pos.n - pos.d) <= maxLen,
        # Filter for same chromosome
        chr.n == chr.d
      ) %>%
      dplyr::mutate(
        start = ifelse(strand.n == "+", pos.n, pos.d),
        end = ifelse(strand.n == "+", pos.d, pos.n)
      ) %>%
      dplyr::select(
        "id" = grp, "chr" = chr.n, "strand" = strand.n, start, end
      )
    
  }))
  

}

#' Determine if pair of reads are mapped
#' @param flag numeric or integer vector of flag codes indicating mapping 
#' status. This integer will be converted into binary bits and decoded to 
#' determine if the flag indicates paired mapping.
#' @description Given flag integer codes, this function returns a logical to 
#' indicate if the pair of reads are both mapped. If one or both reads are 
#' unmapped, then the return is "FALSE".

pair_is_mapped <- function(flag){
  
  #Check if input is in correct format.
  stopifnot( all(is.numeric(flag) | is.integer(flag)) )
  
  # Switch flag codes to binary bit matrix
  x <- matrix(as.integer(intToBits(flag)), ncol = 32, byrow = TRUE)
  
  # Flag codes designate 3rd and 4th bits to indicate unmapped read or mate
  # As long as both are zero, then the pair of reads are both mapped
  rowSums(x[,c(3,4)]) == 0
  
}

#' Determine the alignment is for the read or mate
#' @param flag numeric or integer vector of flag codes indicating mapping 
#' status. This integer will be converted into binary bits and decoded to 
#' determine if the flag indicates read or mate maping.
#' @param output character vector of length 2, indicating the output designation
#' for if the alignment is for the read or the mate.
#' @description Given flag integer codes, this function returns a logical or 
#' character vector to indicate if the alignment is for the read or mate

read_or_mate <- function(flag, output = NULL){
  
  #Check if input is in correct format.
  stopifnot( all(is.numeric(flag) | is.integer(flag)) )
  
  # Switch flag codes to binary bit matrix
  x <- matrix(as.integer(intToBits(flag)), ncol = 32, byrow = TRUE)
  
  # Flag codes designate 7th bit to indicate 1st read (read) and the 8th for mate
  # As long as both are zero, then the pair of reads are both mapped
  if( is.null(output) ){
    return(x[,c(7)] == 1)
  }else{
    return(ifelse(x[,c(7)] == 1, output[1], output[2]))
  }
  
}

# Additional parameters ---- 
# BAM parameters to get from file
bam_params <- c(
  "qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "cigar"
)

# BAM Tags to get from files
bam_tags <- c("MD")

# Import read alignments and filter on input criteria ----
input_hits <- loadBAM(
  bam = args$bam, bai = args$bai, params = bam_params, tags = bam_tags
)

# Top of inputs from alignments
printHead(
  input_hits,
  title = "Head of input alignments",
  caption = sprintf(
    "%1$s total alignments from %2$s reads.", 
    nrow(input_hits), 
    length(unique(input_hits$qname))
  )
)

# Stop if there are no remaining alignments
if( nrow(input_hits) == 0 ){ 
  cat("\nNo alignments in input bam file.\n")
  writeNullOutput(args)
  q()
}

## Initial quality filtering: min percent ID, minimum size, max align start ----
read_hits <- input_hits %>%
  dplyr::mutate(
    pairMapped = pair_is_mapped(flag),
    type = read_or_mate(flag, c("anchor", "adrift"))
  ) %>% 
  dplyr::filter(pairMapped) %>%
  dplyr::mutate(
    clip5p = cntClipped(cigar),
    pctID = calcPctID(cigar, MD)
  ) %>%
  dplyr::filter(
    pctID >= args$minPercentIdentity,
    qwidth >= args$minTempLength,
    clip5p <= args$maxAlignStart
  )

read_wo_pairs_after_init_filter <- read_hits %>%
  dplyr::group_by(qname) %>%
  dplyr::summarise(
    anchors = sum(type == "anchor"), 
    adrifts = sum(type == "adrift")
  ) %>%
  dplyr::filter(anchors == 0 | adrifts == 0) %>%
  dplyr::pull(qname)
  
read_hits <- dplyr::filter(
  read_hits, !qname %in% read_wo_pairs_after_init_filter
)

# Stop if there are no remaining alignments
if( nrow(read_hits) == 0 | dplyr::n_distinct(read_hits$type) == 1 ){
  
  cat(
    "\nNo valid alignments were found within the data given input criteria.\n"
  )

  writeNullOutput(args)
  q()
  
}

## Additional quality filtering: orientation structure, min and max size ----
all_valid_aligns <- with(
    read_hits, 
    processAlignments(
      qname, rname, strand, pos, qwidth, type, 
      refGen = ref_genome, batches = args$batches
    )
  ) %>%
  dplyr::mutate(
    lociPairKey = paste0(
      as.integer(factor(
        paste0(chr, strand, ifelse(strand == "+", start, end))
      )), ":", 
      as.integer(factor(
        paste0(chr, strand, ifelse(strand == "+", end, start))
      ))
    ),
    readPairKey = as.integer(factor(id))
  )

## Remove alternative sequence alignments if requested during input ----
if( !args$keepAltChr ){
  all_valid_aligns <- dplyr::filter(
    all_valid_aligns, !stringr::str_detect(chr, stringr::fixed("_"))
  )
}

### Print out top of valid alignments
printHead(
  all_valid_aligns,
  title = "Head of valid alignments",
  caption = sprintf(
    "%1$s valid alignments from %2$s reads.", 
    nrow(all_valid_aligns), 
    length(unique(all_valid_aligns$id))
  )
)

# Stop if there are no remaining alignments
if( nrow(all_valid_aligns) == 0 ){
  
  cat("\nNo valid alignments were found after QC filtering.\n")
  writeNullOutput(args)
  q()
  
}


## Group alignments into unique and multihit alignments ----
uniq_aligns <- all_valid_aligns %>%
  dplyr::group_by(id) %>%
  dplyr::filter(n() == 1) %>%
  dplyr::ungroup()

multihits <- all_valid_aligns %>%
  dplyr::group_by(id) %>%
  dplyr::filter(n() > 1) %>%
  dplyr::ungroup()

## Recover any reads not captured in the two groups above ----
failed_reads <- input_hits %>%
  dplyr::filter(
    !qname %in% c(unique(uniq_aligns$id), unique(multihits$id))
  )

# Log allocated read counts
cat(
  "\nReads associated with types of alignments:\n",
  "  unique alignments  : ", 
  format(length(unique(uniq_aligns$id)), big.mark = ","), "\n",
  "  multihit alignments: ", 
  format(length(unique(multihits$id)), big.mark = ","), "\n",
  "  chimera artifacts  : ", 
  format(length(unique(failed_reads$qname)), big.mark = ","), "\n"
)


# Bin reads that would map to different loci on the same read (chimeras)
# All unique and multihit templates are mapped successfully to 
# genomic loci, yet some templates are sequenced but do not make it through
# the selection criteria. These templates either do not have alignments to the
# reference genome (anchor or adrift did not align) or map to two distant 
# genomic loci. The latter are termed chimeras and are considered to be 
# artifacts of PCR amplification.

if( !is.null(args$chimeras) ){

  if( args$stat != FALSE ){
    
    add_stat <- data.frame(
      sampleName = sampleName,
      metric = "chimera.reads",
      count = length(unique(failed_reads$qname))
    )
    
    stat <- rbind(stat, add_stat)
    
  }
  
  chimeraData <- list(
    "failed_reads" = failed_reads
  )
  
  writeOutputFile(chimeraData, file = args$chimeras, format = "rds")
  
}


## Write unique (and condensed) output and record stats ----
uniq_sites <- uniq_aligns %>%
  dplyr::mutate(
    width = as.integer(end - start + 1),
    sampleName = sampleName
  ) %>%
  dplyr::select(
    "seqnames" = chr, start, end, width, strand, 
    lociPairKey, readPairKey, sampleName, "ID" = id
  ) %>%
  GenomicRanges::makeGRangesFromDataFrame(
    keep.extra.columns = TRUE, 
    seqinfo = GenomeInfoDb::seqinfo(ref_genome)
  )

writeOutputFile(uniq_sites, file = args$uniqOutput)

### Print out head of uniq_sites for reference.
printHead(
  uniq_sites,
  title = "Head of uniquely mapped genomic loci",
  caption = sprintf(
    paste(
      "Alignments yeilded %1$s unique anchor sites from %2$s", 
      "properly-paired and aligned reads."
    ),
    length(GenomicRanges::reduce(
      GenomicRanges::flank(uniq_sites, -1, start = TRUE), min.gapwidth = 0L
    )),
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

## Generate condensed sites ----
if( !is.null(args$condSites) ){
  
  cond_sites <- condenseSites(
    uniq_sites, keep.cols = "sampleName", list.bp.counts = TRUE
  )
  
  writeOutputFile(cond_sites, file = args$condSites)
  
  printHead(
    cond_sites,
    title = "Head of unique anchor sites",
    caption = sprintf(
      paste(
        "There were %1$s unique anchor sites identified with a total", 
        "of %2$s unique template lengths and %3$s read counts."
      ),
      length(cond_sites),
      sum(cond_sites$fragLengths),
      sum(cond_sites$counts)
    )
  )
  
}


## Write multihits output and record stats ----
if( !is.null(args$multihits) ){
  
  unclustered_multihits <- GenomicRanges::GRanges()
  clustered_multihit_positions <- GenomicRanges::GRangesList()
  clustered_multihit_lengths <- list()
  
  if( nrow(multihits) > 0 ){
    
    #' As the loci are expanded from the coupled_loci object, unique templates 
    #' and readPairKeys are present in the readPairKeys unlisted from the 
    #' paired_loci object.
    multihit_templates <- multihits %>%
      dplyr::mutate(
        width = end - start + 1,
        sampleName = sampleName
      ) %>%
      dplyr::select(
        "seqnames" = chr, start, end, width, strand, 
        lociPairKey, readPairKey, "ID" = id, sampleName
      ) %>%
      GenomicRanges::makeGRangesFromDataFrame(
        keep.extra.columns = TRUE, 
        seqinfo = GenomeInfoDb::seqinfo(ref_genome)
      )
    
    multihit_keys <- multihits %>%
      dplyr::mutate(sampleName = sampleName) %>%
      dplyr::distinct(sampleName, "ID" = id, readPairKey) %>%
      dplyr::select(sampleName, ID, readPairKey)
    
    #' Medians are based on all the potential sites for a given read, which will
    #' be identical for all reads associated with a readPairKey.
    multihit_medians <- round(
      median(GenomicRanges::width(split(
        x = multihit_templates, 
        f = multihit_templates$readPairKey
      )))
    )
    
    multihit_keys$medians <- multihit_medians[
      as.character(multihit_keys$readPairKey)
    ]
    
    multihits_pos <- GenomicRanges::flank(
      x = multihit_templates, width = -1, start = TRUE
    )
    
    multihits_red <- GenomicRanges::reduce(
      x = multihits_pos, min.gapwidth = 5L, with.revmap = TRUE
    )  #! Should make min.gapwidth a option
    
    revmap <- multihits_red$revmap
    
    axil_nodes <- as.character(S4Vectors::Rle(
      values = multihit_templates$readPairKey[min(revmap)], 
      lengths = lengths(revmap)
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
    
    multihits_pos$clusID <- clus_key[
      as.character(multihits_pos$readPairKey), "clusID"
    ]
    
    multihits_pos <- multihits_pos[order(multihits_pos$clusID)]
    
    clustered_multihit_index <- as.data.frame(
      GenomicRanges::mcols(multihits_pos)
    )
    
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
          x = BiocGenerics::unlist(S4Vectors::runValue(multihit_loci_intL)), 
          table = clustered_multihit_index$lociPairKey
        )
      ]
    )
    
    clustered_multihit_positions <- GenomicRanges::split(
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
      match(
        as.character(multihit_keys$readPairKey), 
        readPairKey_cluster_index$readPairKey
      )
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
      BiocGenerics::unlist(multihit_readPair_read_exp)
    ]
    
    unclustered_multihits$ID <- multihit_keys$ID[
      BiocGenerics::unlist(multihit_readPair_read_exp)
    ]
    
    unclustered_multihits$sampleName <- multihit_keys$sampleName[
      BiocGenerics::unlist(multihit_readPair_read_exp)
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
  
  printHead(
    data.frame(
      "multihit_reads" = length(unique(names(unclustered_multihits))),
      "multihit_alignments" = length(unique(unclustered_multihits)),
      "multihit_clusters" = length(clustered_multihit_positions),
      "multihit_lengths" = sum(lengths(clustered_multihit_lengths))
    ),
    title = "Multihit metrics", 
    caption = "Metrics highlighting the observation of multiple aligning reads."
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
