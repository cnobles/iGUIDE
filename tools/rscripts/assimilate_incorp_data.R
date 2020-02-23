#' For anyone reviewing the code below, the following is a small style guide 
#' outlining the various formats for the code. 
#' 
#' Names with "_": objects, inlucding data.frames, GRanges, vectors, ...
#' Names in caMel format: functions or components of objects (i.e. columns 
#' within a data.frame).
#' Names with ".": arguments / options for functions

# Required / highly suggested option parameters and library ----
options(stringsAsFactors = FALSE, scipen = 99, width = 180)
suppressMessages(library("magrittr"))
suppressMessages(library("iguideSupport"))

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
  "-u", "--umitags", nargs = 1, type = "character",
  help = paste(
    "Path to directory with associated fasta files containing read specific",
    "random captured sequences. Directory should contain files with file names",
    "like *.umitags.fasta."
  )
)

parser$add_argument(
  "-m", "--multihits", nargs = 1, type = "character",
  help = paste(
    "Path to directory with associated multihit files (*.multihits.rds) as",
    "produced by coupling alignment output files."
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
  "--readNamePattern", nargs = 1, type = "character", 
  default = "[\\w\\:\\-\\+]+", 
  help = "Regular expresion capturing the read name for a given sequence."
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

## Remove output file(s) if existing
if( args$stat != FALSE ){
  output_files <- c(args$output, args$stat)
}else{
  output_files <- c(args$output)
}

if( any(sapply(output_files, file.exists)) ){
  null <- lapply(output_files, unlink)
}

# Log inputs
cat("\nAssimilate Inputs:\n")
print(
  x = data.frame(input_table),
  right = FALSE, 
  row.names = FALSE
)


# Get versioning ----
soft_version <- as.character(read.delim(
  file = file.path(root_dir, ".version"), header = FALSE))

build_version <- list.files(file.path(root_dir, "etc")) %>%
  grep(pattern = "build.b[0-9\\.]+.*", x = ., value = TRUE) %>%
  stringr::str_extract(pattern = "b[0-9]+\\.[0-9]+\\.[0-9]+")


# Inputs and parameters ----
# Run parameters and sample parameters
config <- yaml::yaml.load_file(args$config)

## These parameters are dictate part of the following analysis if multihit 
## alignments are to be considered in the analysis.

upstream_dist <- config$upstreamDist
downstream_dist <- config$downstreamDist
pile_up_min <- config$pileUpMin

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


# Load input data ----
## Unique sites ----
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
  
  multihit_files <- list.files(path = args$multihit, full.names = TRUE)
  
  mulithit_files <- multihit_files[
    stringr::str_detect(mulithit_files, ".multihits.rds")
  ]
  
  multi_reads <- unlist(GRangesList(lapply(mulithit_files, function(x){
    
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

if( !is.null(args$umitags) ){
  
  umitag_files <- list.files(path = args$umitags, full.names = TRUE)
  
  umitag_files <- umitag_files[
    stringr::str_detect(umitag_files, ".umitags.fasta")
  ]
  
  umitags <- lapply(umitag_files, ShortRead::readFasta)
  umitags <- serialAppendS4(umitags)
  umitag_read_ids <- stringr::str_extract(
    as.character(ShortRead::id(umitags)),
    args$readNamePattern
  )
  
  reads$umitag <- as.character(ShortRead::sread(umitags))[
    match(reads$ID, umitag_read_ids)
  ]
  
}

# Generate stats if requested ----
## If requested, generate stats from the analysis for qc.

if( args$stat != FALSE ){

  stat <- reads %>% 
    dplyr::group_by(sampleName) %>% 
    dplyr::summarise(
      reads = dplyr::n_distinct(ID),
      aligns = dplyr::n_distinct(seqnames, start, end, strand),
      loci = dplyr::n_distinct(
        seqnames, strand, ifelse(strand == "+", start, end)
      )
    ) %>% 
    tidyr::gather(key = "type", value = "value", -sampleName)
  
  write.table(
    x = stat, file = args$stat, 
    sep = ",", row.names = FALSE, 
    col.names = FALSE, quote = FALSE
  )
  
}

# Output data ----
## rds file that can be read into evaluation or reports or loaded into a 
## database with some additional scripting.
fmt_reads <- reads %>%
  dplyr::select(-lociPairKey, -readPairKey)

output_file <- list(
  "soft_version" = soft_version,
  "build_version" = build_version,
  "config" = config,
  "reads" = fmt_reads
)
  
saveRDS(output_file, file = args$output)

if( all(sapply(output_files, file.exists)) ){
  message("Successfully completed script.")
}else{
  stop("Check output, it not detected after assimilating.")
}

q(status = 0)
