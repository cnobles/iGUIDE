#' For anyone reviewing the code below, the following is a small style guide 
#' outlining the various formats for the code. 
#' 
#' Names with "_": objects, inlucding data.frames, GRanges, vectors, ...
#' Names in caMel format: functions or components of objects (i.e. columns 
#' within a data.frame).
#' Names with ".": arguments / options for functions

# Required / highly suggested option parameters and library ----
options(stringsAsFactors = FALSE, scipen = 99, width = 999)
suppressMessages(library("magrittr"))

# Set up and gather command line arguments ----
parser <- argparse::ArgumentParser(
  description = "Assimilate incorporation data from iGUIDE pipeline.",
  usage = paste(
    "Rscript assimilate_incorp_data.R <uniqSites> -o <output> -c <config>",
    "[-h/--help, -v/--version] [optional args]"
  )
)

parser$add_argument(
  "uniqSites", nargs = 1, type = "character",
  help = paste(
    "Unique sites output from blatCoupleR. The output from an entire run can",
    "be concatenated together as a single input."
  )
)

parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", required = TRUE,
  help = "Output file name in .rds format."
)

parser$add_argument(
  "-c", "--config", nargs = 1, type = "character", required = TRUE,
  help = "Run specific config file in yaml format."
)

parser$add_argument(
  "-u", "--umitags", nargs = "+", type = "character",
  help = paste(
    "Path(s) to associated fasta files containing read specific",
    "random captured sequences. Multiple file paths can be separated by",
    "a space."
  )
)

parser$add_argument(
  "-m", "--multihits", nargs = "+", type = "character",
  help = paste(
    "Path(s) to associated multihit files (.rds) as produced by coupling",
    "BLAT output files. Multiple file paths can be separated by a space."
  )
)

parser$add_argument(
  "--stat", nargs = 1, type = "character", default = FALSE, 
  help = paste(
    "File name to be written in output directory of read couts for each",
    "sample. CSV file format. ie. test.stat.csv."
  )
)

parser$add_argument(
  "--iguide_dir", nargs = 1, type = "character", default = "IGUIDE_DIR",
  help = "iGUIDE install directory path, do not change for normal applications."
)


# Set arguments with parser ----
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

if( !dir.exists(args$iguide_dir) ){
  root_dir <- Sys.getenv(args$iguide_dir)
}else{
  root_dir <- args$iguide_dir
}

if( !dir.exists(root_dir) ){
  stop(paste0("\n  Cannot find install path to iGUIDE: ", root_dir, ".\n"))
}else{
  args$iguide_dir <- root_dir
}


input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(
    seq_along(args), 
    function(i) paste(args[[i]], collapse = ", ")
  )
)

input_table <- input_table[
  match(c(
    "uniqSites :", "output :", "config :", "umitags :", "multihits :", "stat :",
    "iguide_dir :"
    ), 
    input_table$Variables
  ),
]

# Log inputs
cat("\nAssimilate Inputs:\n")
print(
  x = data.frame(input_table),
  right = FALSE, 
  row.names = FALSE
)


# Source supporting functions ----
# Set the code_dir from the commandline call to the script
code_dir <- dirname(sub(
  pattern = "--file=", 
  replacement = "", 
  x = grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
))

## Load in supporting functions for the analysis
source(file.path(code_dir, "supporting_scripts/iguide_support.R"))
source(file.path(code_dir, "supporting_scripts/nucleotideScoringMatrices.R"))

# Inputs and parameters ----
# Run parameters and sample parameters
config <- yaml::yaml.load_file(args$config)

sample_info_path <- file.path(args$iguide_dir, config$Sample_Info)

if( !(file.exists(sample_info_path) | file.exists(config$Sample_Info)) ){
  
  stop(
    "\n  Specified Sample Info file not found: ", config$Sample_Info, "\n"
  )
  
}else if( !file.exists(sample_info_path) ){
  
  sample_info_path <- file.path(config$Sample_Info)
  
}

sample_info <- data.table::fread(input = sample_info_path, data.table = FALSE)

submat <- banmat()


# Load reference genome ----
## Load a reference genome from a fasta file or a BSGenome reference.
## Script stops if ref genome is not available

if( grepl(".fa", config$Ref_Genome) ){
  
  if( !(
    file.exists(file.path(args$iguide_dir, config$Ref_Genome)) | 
      file.exists(config$Ref_Genome)
  ) ){
    stop(
      "\n  Specified reference genome file not found: ", config$Ref_Genome, "\n"
    )
  }
  
  ref_file_type <- ifelse(grepl(".fastq", config$Ref_Genome), "fastq", "fasta")
  
  
  if( file.exists(
    file.path(args$iguide_dir, config$Ref_Genome) 
  ) ){
    
    ref_genome <- Biostrings::readDNAStringSet(
      filepath = file.path(args$iguide_dir, config$Ref_Genome),
      format = ref_file_type
    )
    
  }else{
    
    ref_genome <- Biostrings::readDNAStringSet(
      filepath = config$Ref_Genome,
      format = ref_file_type
    )
    
  }
  
}else{
  
  genome <- grep(
    pattern = config$Ref_Genome, 
    x = unique(BSgenome::installed.genomes()), 
    value = TRUE
  )
  
  if( length(genome) == 0 ){
    
    cat("\nInstalled genomes include:\n")
    print(unique(BSgenome::installed.genomes()))
    cat("\n  Selected reference genome not in list.")
    stop("\n  Genome not available.\n")
    
  }else if( length(genome) > 1 ){
    
    cat("\nInstalled genomes include:\n")
    print(unique(BSgenome::installed.genomes()))
    cat(
      "\n  Please be more specific about reference genome.",
      "Multiple matches to input."
    )
    stop("\n  Multiple genomes requested.\n")
    
  }
  
  suppressMessages(library(genome, character.only = TRUE))
  
  ref_genome <- get(genome)
  
}


# Incorporation site parameters ----
## These parameters are pulled straight from the run config file and describe 
## how the following analysis will be conducted.

upstream_dist <- config$upstreamDist
downstream_dist <- config$downstreamDist
max_target_mismatch <- config$maxTargetMismatch
pile_up_min <- config$pileUpMin
on_target_sites <- config$On_Target_Sites 

# Load data related to how samples were processed ----
## The treatment object dictates how each sample was treated, or which guide 
## RNAs where used on which samples. This is important in the results 
## interpretation.

treatment <- config$Treatment

if( any(grepl("sampleInfo:", treatment[1])) ){
  
  info_col <- match(
    x = stringr::str_extract(string = treatment[1], pattern = "[\\w]+$"), 
    table = names(sample_info)
  )
  
  if( length(info_col) != 1 ){
    stop("\n  Cannot parse treatment data. Check config yaml and sampleInfo.\n")
  }
  
  treatment_df <- data.frame(
    sampleName = sample_info$sampleName, 
    treatment = sample_info[,info_col]
  )
  
  treatment_df$specimen <- stringr::str_extract(
    string = treatment_df$sampleName, pattern = "[\\w]+"
  )
  
  treatment_df <- unique(treatment_df[,c("specimen", "treatment")])
  treatment <- strsplit(as.character(treatment_df$treatment), ";")
  names(treatment) <- treatment_df$specimen
  
}else if( any(grepl("all", names(treatment))) ){
  
  treatment_df <- data.frame(
    sampleName = sample_info$sampleName, 
    treatment = unique(unlist(treatment))
  )
  
  treatment_df$specimen <- stringr::str_extract(
    string = treatment_df$sampleName, 
    pattern = "[\\w]+"
  )
  
  treatment_df <- unique(treatment_df[,c("specimen", "treatment")])
  treatment <- strsplit(as.character(treatment_df$treatment), ";")
  names(treatment) <- treatment_df$specimen
  
}else{
  
  treatment_df <- data.frame(
    "specimen" = names(treatment), 
    "treatment" = sapply(treatment, paste, collapse = ";")
  )
  
}


## Identify the treatment of nucleases as well

nuclease <- config$Nuclease

if( any(grepl("sampleInfo:", nuclease[1])) ){
  
  info_col <- match(
    x = stringr::str_extract(string = nuclease[1], pattern = "[\\w]+$"), 
    table = names(sample_info)
  )
  
  if( length(info_col) != 1 ){
    stop("\n  Cannot parse nuclease data. Check config yaml and sampleInfo.\n")
  }
  
  nuclease_df <- data.frame(
    sampleName = sample_info$sampleName, 
    nuclease = sample_info[,info_col]
  )
  
  nuclease_df$specimen <- stringr::str_extract(
    string = nuclease_df$sampleName, pattern = "[\\w]+"
  )
  
  nuclease_df <- unique(nuclease_df[,c("specimen", "nuclease")])
  nuclease <- strsplit(as.character(nuclease_df$nuclease), ";")
  names(nuclease) <- nuclease_df$specimen
  
}else if( any(grepl("all", names(nuclease))) ){
  
  nuclease_df <- data.frame(
    sampleName = sample_info$sampleName, 
    nuclease = unique(unlist(nuclease))
  )
  
  nuclease_df$specimen <- stringr::str_extract(
    string = nuclease_df$sampleName, 
    pattern = "[\\w]+"
  )
  
  nuclease_df <- unique(nuclease_df[,c("specimen", "nuclease")])
  nuclease <- strsplit(as.character(nuclease_df$nuclease), ";")
  names(nuclease) <- nuclease_df$specimen
  
}else{
  
  nuclease_df <- data.frame(
    "specimen" = names(nuclease), 
    "nuclease" = sapply(nuclease, paste, collapse = ";")
  )
  
}

nuclease_treaments <- dplyr::left_join(
  nuclease_df, treatment_df, by = "specimen"
)

target_combn <- structure(
  strsplit(nuclease_treaments$treatment, ";"), 
  names = nuclease_treaments$specimen
)

combn_tbl <- data.frame(
    nuclease = nuclease_treaments$nuclease[
      as.vector(match(
        S4Vectors::Rle(names(target_combn), lengths(target_combn)), 
        nuclease_treaments$specimen
      ))
    ],
    target = unlist(target_combn),
    row.names = NULL
  ) %>%
  dplyr::filter(target != "Mock") %>%
  dplyr::distinct()


# Load target sequences and sample metadata ----
## Identify the target sequences used for the analysis and build an object 
## use further on in processing to analyse the samples.

target_seqs <- lapply(config$Target_Sequences, toupper)

pam_seq <- lapply(
  unique(unlist(nuclease)), 
  function(x) toupper(config$Nuclease_Profiles[[x]]$PAM)
)

names(pam_seq) <- unique(unlist(nuclease))
  
combn_tbl <- combn_tbl %>%
  dplyr::mutate(
    sequence = target_seqs[target],
    PAM = pam_seq[nuclease]
  )

# Log combination treatment table
cat("\nTarget Sequence Table:\n")
print(combn_tbl, right = FALSE, row.names = FALSE)


# Log treatment table
cat("\nSpecimen Target Treatment:\n")
null <- lapply(seq_along(treatment), function(i){
  cat("  ", names(treatment)[i], ":\n", sep = "")
  cat("    nuclease : ", nuclease[[i]], "\n", sep = "")
  cat("    targets  : ", paste(treatment[[i]], collapse = ", "), "\n", sep = "")
})


# Load input data ----
# Unique sites ----
## This object is the alignment positions for the sequences / reads that only 
## aligned to a single location on the reference genome.

reads <- data.table::fread(
  input = args$uniqSites, data.table = FALSE, stringsAsFactors = FALSE
)


# Multihits if requested ----
## Multihits are alignments that legitimately appear in multiple locations
## across the reference genome. These can be more difficult to interpret but are
## an option for this software. The user should be familiar and cautious of 
## alignment artifacts if using multihit data.

if( all(!is.null(args$multihits)) ){
  
  uniq_reads <- GenomicRanges::makeGRangesFromDataFrame(
    df = reads, 
    keep.extra.columns = TRUE, 
    seqinfo = GenomeInfoDb::seqinfo(ref_genome)
  )
  
  multi_reads <- unlist(GRangesList(lapply(args$multihits, function(x){
    
    multi <- readRDS(x)
    GenomeInfoDb::seqinfo(multi$unclustered_multihits) <- 
      GenomeInfoDb::seqinfo(ref_genome)
    
    if( length(multi$unclustered_multihits) > 0 ){
      
      GenomicRanges::mcols(multi$unclustered_multihits) <- 
        GenomicRanges::mcols(multi$unclustered_multihits)[
          ,c(names(GenomicRanges::mcols(uniq_reads)))
      ]
      
    }else{
      
      GenomicRanges::mcols(multi$unclustered_multihits) <- 
        GenomicRanges::mcols(uniq_reads)[
          0, c(names(GenomicRanges::mcols(uniq_reads)))
      ]
      
    }
    
    multi$unclustered_multihits
    
  })))
  
  comb_reads <- c(uniq_reads, multi_reads)
  
  GenomicRanges::mcols(comb_reads)$type <- rep(
    c("uniq", "multi"), c(length(uniq_reads), length(multi_reads))
  )
  
  GenomicRanges::mcols(comb_reads)$clus.id <- pileupCluster(
    gr = comb_reads, 
    grouping = "sampleName", 
    return = "ID"
  )
  
  filt_multi_reads <- dplyr::bind_rows(lapply(
    split(comb_reads, comb_reads$sampleName), 
    function(x){
      
      uniq_id <- unique(x$clus.id[x$type == "uniq"])
      multi_id <- unique(x$clus.id[x$type == "multi"])
      y <- x[x$type == "multi" & x$clus.id %in% intersect(uniq_id, multi_id)]
      mcols(y)$clus.id <- NULL
      
      if( length(y) > 0 ){
        
        contrib_amt <- 1 / table(mcols(y)$ID)
        GenomicRanges::mcols(y)$contrib <- 
          as.numeric(contrib_amt[GenomicRanges::mcols(y)$ID])
        
      }
      
      as.data.frame(y, row.names = NULL) %>%
        dplyr::mutate(
          seqnames = as.character(seqnames), 
          strand = as.character(strand)
        )
      
    }
  ))
  
  reads <- dplyr::mutate(reads, type = "uniq", contrib = 1) %>%
    dplyr::bind_rows(., filt_multi_reads)
  
}else{
  
  reads <- dplyr::mutate(reads, type = "uniq", contrib = 1)
  
}

# Print out stats during analysis.
cat("\nTabulation of aligned reads per specimen:\n")
temp_table <- table(stringr::str_extract(reads$sampleName, "[\\w]+"))

print(
  data.frame(
    "Specimen" = names(temp_table), 
    "Aligned_Reads" = format(as.numeric(temp_table), big.mark = ",")
  ),
  right = FALSE,
  row.names = FALSE
)

rm(temp_table)

# Umitags or captured random sequences ----
## Unique molecular index tags, or UMItags, are random sequences appended to the
## index 2 read. They are 8 or so nucleotides and are combined with the terminal 
## breakpoint sequence to be potentially used for a higher dynamic range 
## abundance measure. While ideal in theory, practice has identified these 
## sequences skewing with read counts and an over abundance of sharing of the 
## random sequence between difference breakpoints. Interpretation of UMItag 
## based abundances should be interpreted with caution as they are prone / 
## susceptable to PCR artifacts.

if( all(!is.null(args$umitags)) ){
  
  umitags <- lapply(args$umitags, ShortRead::readFasta)
  umitags <- serialAppendS4(umitags)
  
  reads$umitag <- as.character(ShortRead::sread(umitags))[
    match(reads$ID, as.character(ShortRead::id(umitags)))
  ]
  
}

# Process input data ----
# Format input alignments ----
## All alignments

algnmts <- dplyr::mutate(
  reads, 
  specimen = stringr::str_extract(sampleName, "[\\w]+")
)

## Determine abundance metrics, with or without UMItags
if( config$UMItags & !is.null(args$umitags) ){
  
  algnmts <- dplyr::arrange(algnmts, desc(contrib)) %>%
    dplyr::group_by(seqnames, start, end, strand, specimen, sampleName) %>%
    dplyr::summarise(
      count = sum(contrib),
      umitag = sum(as.integer(!duplicated(umitag[!is.na(umitag)])) * contrib),
      contrib = max(contrib)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      seqnames, start, end, strand, specimen, 
      sampleName, count, umitag, contrib
    ) %>%
    as.data.frame()
  
}else{
  
  algnmts <- dplyr::group_by(
      algnmts, seqnames, start, end, strand, specimen, sampleName
    ) %>%
    dplyr::summarize(count = sum(contrib), contrib = max(contrib)) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      seqnames, start, end, strand, specimen, sampleName, count, contrib
    ) %>%
    as.data.frame()
  
}

## Generate a sample table of the data for log purposes
sample_index <- ifelse(nrow(algnmts) > 10, 10, nrow(algnmts))
sample_index <- sample(seq_len(nrow(algnmts)), sample_index, replace = FALSE)

cat("\nSample of aligned templates:\n")

print(
  data.frame(algnmts[sample_index,]),
  right = FALSE,
  row.names = FALSE
)

cat(paste0("\nNumber of templates: ", nrow(algnmts), "\n"))

rm(sample_index)

## Transform the data into a GRanges object
algnmts_gr <- GenomicRanges::GRanges(
  seqnames = algnmts$seqnames,
  ranges = IRanges::IRanges(start = algnmts$start, end = algnmts$end),
  strand = algnmts$strand,
  seqinfo = GenomeInfoDb::seqinfo(ref_genome)
)

if( config$UMItags & !is.null(args$umitags) ){
  
  GenomicRanges::mcols(algnmts_gr) <- dplyr::select(
    algnmts, specimen, sampleName, count, umitag, contrib
  )
  
}else{
  
  GenomicRanges::mcols(algnmts_gr) <- dplyr::select(
    algnmts, specimen, sampleName, count, contrib
  )
  
}

# Analyze alignments ----
## Identify groups of alignments or pileups of aligned fragments
## These pileups give strong experimental evidence of directed incorporation of
## the dsODN into a region. Initially, pileups are identified and then checked 
## for pairing, or if there is another pileup on the opposite strand in close 
## proximity.
algnmts_gr$clus.ori <- pileupCluster(
  gr = algnmts_gr, 
  grouping = "specimen", 
  maxgap = 0L, 
  return = "simple"
)

algnmts_gr$paired.algn <- identifyPairedAlgnmts(
  gr = algnmts_gr, 
  grouping = "specimen", 
  maxgap = upstream_dist*2
)

algnmts_grl <- split(algnmts_gr, unlist(nuclease)[algnmts_gr$specimen])

annot_clust_info <- dplyr::bind_rows(lapply(
  seq_along(algnmts_grl), 
  function(i, grl){
    
    gr <- grl[[i]]
    nuc <- names(grl)[i]
    
    if( !nuc %in% names(config$Nuclease_Profiles) ){
      nuc_profile <- NULL
    }else{
      nuc_profile <- config$Nuclease_Profiles[[nuc]]
    }
    
    ## Create a GRange with only the unique cluster origins
    split_clus_id <- stringr::str_split(
      string = unique(paste0(gr$specimen, ":", gr$clus.ori)), 
      pattern = ":", 
      simplify = TRUE
    )
    
    algn_clusters <- GenomicRanges::GRanges(
      seqnames = split_clus_id[,2],
      ranges = IRanges::IRanges(
        start = as.numeric(split_clus_id[,4]), width = 1
      ),
      strand = split_clus_id[,3],
      seqinfo = GenomeInfoDb::seqinfo(ref_genome)
    )
    
    algn_clusters$specimen <- split_clus_id[,1]
    algn_clusters$clus.ori <- vcollapse(split_clus_id[, 2:4], sep = ":")
    
    algn_clusters$clus.seq <- getSiteSeqs(
      gr = algn_clusters, 
      upstream.flank = upstream_dist, 
      downstream.flank = downstream_dist, 
      ref.genome = ref_genome
    )
    
    ## Identify which target sequences binding near clusters
    if( !is.null(nuc_profile) ){
      
      algn_clusters <- compareTargetSeqs(
        gr.with.sequences = algn_clusters, 
        seq.col = "clus.seq", 
        target.seqs = target_seqs,
        tolerance = max_target_mismatch,
        nuc.profile = nuc_profile,
        submat = submat, 
        upstream.flank = upstream_dist, 
        downstream.flank = downstream_dist
      )
    
    }else{
      
      algn_clusters$target.match <- "No_valid_match"
      algn_clusters$target.mismatch <- NA
      algn_clusters$target.score <- NA
      algn_clusters$aligned.sequence <- NA
      algn_clusters$edit.site <- NA
      
    }
    
    as.data.frame(GenomicRanges::mcols(algn_clusters))
    
  },
  grl = algnmts_grl
))
  

## Merge the target sequence alignment information from the clusters back to all
## unique alignments
algnmts <- as.data.frame(merge(
  x = as.data.frame(algnmts_gr), 
  y = dplyr::select(annot_clust_info, -clus.seq),
  by = c("specimen", "clus.ori")
))

## Change guideRNA.match to No_Valid_Match if an inappropriate gRNA is annotated
algnmts$target.match <- filterInappropriateComparisons(
  guideRNA.match = algnmts$target.match, 
  specimen = algnmts$specimen, 
  treatment = treatment
)

## Fragment pileups, paired clustering, and guideRNA alignments have been used 
## to characterize the incorporation sites analyzed here. Each metric will be 
## used to create a list of incorporation sites that may be nuclease cut sites. 
## The following identifies which alignments are associated with each of these 
## criteria.
tbl_clus_ori <- algnmts %>% 
  dplyr::group_by(specimen, clus.ori) %>%
  dplyr::filter(n() >= pile_up_min) %>%
  dplyr::ungroup() %$%
  table(clus.ori)

idx_clus_ori <- which(algnmts$clus.ori %in% names(tbl_clus_ori))

tbl_paired_algn <- algnmts %>%
  dplyr::filter(!is.na(paired.algn)) %$%
  table(paired.algn)

idx_paired_algn <- which(algnmts$paired.algn %in% names(tbl_paired_algn))

idx_matched <- which(algnmts$target.match != "No_valid_match")

idx_combined <- sort(unique(c(idx_clus_ori, idx_paired_algn, idx_matched)))

idx_df <- data.frame(
  "Type" = c("PileUp", "Paired", "gRNA_Matched", "Combined"),
  "Counts" = sapply(
    list(idx_clus_ori, idx_paired_algn, idx_matched, idx_combined), 
    length
  )
)

cat("\nTable of uniquely aligned template counts:\n")
print(idx_df, right = FALSE, row.names = FALSE) 
cat(paste0("\nTotal number of templates: ", nrow(algnmts), "\n"))

probable_algns <- algnmts[idx_combined,]

probable_algns$on.off.target <- ifelse(
  probable_algns$edit.site %in% expandPosStr(unlist(on_target_sites)), 
  "On-target", 
  "Off-target"
)

cat("\nOn / Off target alignment counts:\n")
print(table(probable_algns$on.off.target))


## Create summary and output formated object related to each of the criteria for
## edited site detection.

## Matched alignments
matched_algns <- probable_algns[
  probable_algns$target.match != "No_valid_match",
]

matched_summary <- matched_algns %>%
  dplyr::mutate(
    target.match = stringr::str_replace(
      string = target.match, 
      pattern = "\\:\\([\\w]+\\)$",
      replacement = ""
    )
  ) %>%
  dplyr::group_by(
    specimen, edit.site, aligned.sequence, target.match, target.mismatch
  )

if( config$UMItags & !is.null(args$umitags) ){
  
  matched_summary <- dplyr::summarise(
    matched_summary,
    on.off.target = paste(sort(unique(on.off.target)), collapse = ";"),
    paired.algn = paste(sort(unique(paired.algn)), collapse = ";"),
    count = sum(count), 
    umitag = sum(umitag),
    algns = sum(contrib),
    orient = paste(sort(unique(as.character(strand))), collapse = ";")
  )
  
}else{
  
  matched_summary <- dplyr::summarise(
    matched_summary,
    on.off.target = paste(sort(unique(on.off.target)), collapse = ";"),
    paired.algn = paste(sort(unique(paired.algn)), collapse = ";"),
    count = sum(count), 
    algns = sum(contrib),
    orient = paste(sort(unique(as.character(strand))), collapse = ";")
  )
  
}

matched_summary <- dplyr::ungroup(matched_summary) %>% 
  dplyr::arrange(specimen, target.match, desc(algns)) %>%
  as.data.frame()

## Paired alignments
paired_algns <- probable_algns[
  probable_algns$paired.algn %in% names(tbl_paired_algn),
]

paired_regions <- paired_algns %>%
  dplyr::group_by(specimen, paired.algn, strand) %>%
  dplyr::mutate(pos = ifelse(strand == "+", min(start), max(end))) %>%
  dplyr::group_by(specimen, paired.algn)

if( config$UMItags & !is.null(args$umitags) ){
  
  paired_regions <- dplyr::summarise(
      paired_regions,
      seqnames = unique(seqnames),
      start = min(pos), 
      end = max(pos), 
      mid = start + (end-start)/2,
      strand = "*", 
      width = end - start, 
      count = sum(count), 
      umitag = sum(umitag), 
      algns = sum(contrib)
    ) %>%
    dplyr::ungroup()
  
}else{
  
  paired_regions <- dplyr::summarise(
      paired_regions,
      seqnames = unique(seqnames),
      start = min(pos), 
      end = max(pos), 
      mid = start + (end-start)/2,
      strand = "*", 
      width = end - start, 
      count = sum(count), 
      algns = sum(contrib)
    ) %>%
    dplyr::ungroup()
  
}

if( nrow(paired_regions) > 0 ){
  
  paired_regions <- paired_regions %>%
    dplyr::group_by(specimen, paired.algn) %>%
    dplyr::mutate(
      on.off.target = ifelse(
        any(sapply(
          expandPosStr(unlist(on_target_sites[
            which(
              stringr::str_extract(
                names(on_target_sites), "[\\w\\-\\_\\.]+") %in% 
              treatment[[specimen]]
            )
          ])),
          function(x, seq, st, en){
            
            match_seq <- seq == stringr::str_extract(x, "[\\w]+")
            
            within_start <- st <= 
              as.numeric(stringr::str_extract(x, "[\\w]+$")) + downstream_dist
            
            within_end <- en >= 
              as.numeric(stringr::str_extract(x, "[\\w]+$")) - downstream_dist
            
            match_seq & within_start & within_end
            
          }, 
          seq = seqnames, 
          st = start, 
          en = end
        )), 
        "On-target", 
        "Off-target"
      )
    ) %>%
    dplyr::ungroup() %>% 
    as.data.frame()
  
}else{
  
  paired_regions <- dplyr::mutate(
    paired_regions,
    on.off.target = vector(mode = "character")
  )
  
}

## Pile up alignments
pile_up_algns <- probable_algns[
  probable_algns$clus.ori %in% names(tbl_clus_ori),
]

pile_up_summary <- pile_up_algns %>%
  dplyr::mutate(
    target.match = stringr::str_replace(
      string = target.match, 
      pattern = "\\:\\([\\w]+\\)$",
      replacement = ""
    )
  ) %>%
  dplyr::group_by(specimen, clus.ori)

if( config$UMItags & !is.null(args$umitags) ){
  
  pile_up_summary <- dplyr::summarise(
    pile_up_summary,
    on.off.target = paste(sort(unique(on.off.target)), collapse = ";"),
    paired.algn = paste(sort(unique(paired.algn)), collapse = ";"),
    count = sum(count), 
    umitag = sum(umitag),
    algns = sum(contrib)
  )
  
}else{
  
  pile_up_summary <- dplyr::summarise(
    pile_up_summary,
    on.off.target = paste(sort(unique(on.off.target)), collapse = ";"),
    paired.algn = paste(sort(unique(paired.algn)), collapse = ";"),
    count = sum(count), 
    algns = sum(contrib)
  )
  
}

pile_up_summary <- dplyr::ungroup(pile_up_summary) %>% 
  dplyr::arrange(specimen, desc(algns)) %>%
  as.data.frame()


# Generate stats if requested ----
## If requested, generate stats from the analysis for qc.

if( args$stat != FALSE ){
  
  stat_summary <- function(x, y){
    
    x %>%
      dplyr::mutate(metric = y) %>%
      dplyr::group_by(sampleName, metric) %>%
      dplyr::summarize(count = sum(contrib)) %>%
      dplyr::ungroup()
    
  }
  
  total_stat <- stat_summary(algnmts, "total.algns")
  combined_stat <- stat_summary(probable_algns, "combined.algns")
  pileup_stat <- stat_summary(pile_up_algns, "pileup.algns")
  paired_stat <- stat_summary(paired_algns, "paired.algns")
  matched_stat <- stat_summary(matched_algns, "matched.algns")
  
  on_tar_stat <- dplyr::filter(
      matched_algns, on.off.target == "On-target"
    ) %>%
    stat_summary("ontarget.algns")
  
  off_tar_stat <- dplyr::filter(
      matched_algns, on.off.target == "Off-target"
    ) %>%
    stat_summary("offtarget.algns")
  
  stat <- dplyr::bind_rows(
    total_stat, combined_stat, pileup_stat, paired_stat, 
    matched_stat, on_tar_stat, off_tar_stat)
  
  write.table(
    x = stat, file = args$stat, 
    sep = ",", row.names = FALSE, 
    col.names = FALSE, quote = FALSE
  )
  
}

# Output data composition ----
## rds file that can be read into reports or loaded into a database with some 
## additional scripting.
data_comp <- list(
  "algnmts" = algnmts, 
  "probable_algns" = probable_algns,
  "matched_algns" = matched_algns,
  "matched_summary" = matched_summary,
  "paired_algns" = paired_algns,
  "paired_regions" = paired_regions,
  "pile_up_algns" = pile_up_algns,
  "pile_up_summary" = pile_up_summary
)

saveRDS(data_comp, file = args$output)

if( file.exists(args$output) ){
  message("Successfully completed script.")
}else{
  stop("Check output, it not detected after assimilating.")
}

q(status = 0)
