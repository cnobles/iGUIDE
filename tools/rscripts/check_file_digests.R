# Import yaml file with test paths and md5 digests
# Required libraries stringr, digest, Biostrings, data.table, yaml, pander
options(stringsAsFactors = FALSE, scipen = 99, warn = -1, window = 999)
suppressMessages(library("magrittr"))

# Set up and gather command line arguments ----
parser <- argparse::ArgumentParser(
  description = "Test checksums of files from a yaml input.",
  usage = "Rscript tools/rscripts/check_file_digests.R <yaml.input> <options>"
)

parser$add_argument(
  "yaml", nargs = 1, type = "character",
  help = "Yaml containing file paths and checksums (md5). ie. sim.test.yml"
)

parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", default = FALSE,
  help = "Output file name .csv, .tsv, or .rds format."
)

parser$add_argument(
  "-v", "--verbose", action = "store_true", 
  help = "Turns on diagnositc-based messages."
)

parser$add_argument(
  "--install_path", nargs = 1, type = "character", default = "IGUIDE_DIR",
  help = "iGUIDE install directory path, do not change for normal applications."
)

# Set arguments with parser ----
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

root_dir <- Sys.getenv("IGUIDE_DIR")
args$install_path <- root_dir

code_dir <- dirname(sub(
  pattern = "--file=", 
  replacement = "", 
  x = grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
))

input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(
    seq_along(args), 
    function(i) paste(args[[i]], collapse = ", ")
  )
)

input_table <- input_table[
  match(
    c("yaml :", "output :", "verbose :", "install_path :"), 
    input_table$Variables
  ),
]

## Log inputs
if( args$verbose ){
  
  cat("List Inputs")
  pander::pandoc.table(
    data.frame(input_table),
    justify = "left", 
    row.names = FALSE,
    style = "simple",
    split.table = Inf
  )
  
}


# Additional functions ----
readFile <- function(path, root){
  
  if( !file.exists(path) ){
    root_path <- file.path(root, path)
    if( !file.exists(root_path) ){
      stop("Cannot find file:", path)
    }else{
      path <- root_path
    }
  }
  
  # Read extension form path
  ext <- stringr::str_extract(path, "[\\w]+$")
  supported_ext <- c("tsv", "csv", "gz", "fasta", "fastq", "rds")
  
  stopifnot( ext %in% supported_ext )
  
  # Check additional extension if compressed
  if( ext == "gz" ){
    ext2 <- stringr::str_extract(path, "[\\w]+.gz")
    ext2 <- gsub(".gz", "", ext2)
    stopifnot( ext2 %in% supported_ext )
  }else{
    ext2 <- NA
  }
  
  exts <- c(ext, ext2)
  exts <- exts[!is.na(exts)]
  
  # Read in methods based on inputs.
  if( any(exts %in% c("tsv", "csv")) ){
    
    if( ext == "gz" ) path <- paste0("zcat ", path)
    return(data.table::fread(input = path, data.table = FALSE))
    
  }else if( any(stringr::str_detect(exts, "fast")) ){
    
    return(Biostrings::readDNAStringSet(path))
    
  }else{
    
    return(readRDS(path))
    
  }
  
}

# Load inputs ----
config <- yaml::yaml.load_file(args$yaml)
paths <- lapply(config, "[[", "path")
data_objs <- lapply(paths, readFile, root = args$install_path)

# Check digests ----
test_digests <- sapply(data_objs, digest::digest)
check_digests <- sapply(config, "[[", "md5")

df <- data.frame(
  "file_name" = sapply(config, "[[", "name"),
  "md5_standard" = check_digests,
  "md5_tested" = test_digests,
  "outcome" = ifelse(test_digests == check_digests, "pass", "FAIL")
)

# Log output if requested ----
if( args$verbose ){
  
  cat("\nList of Outcomes")
  pander::pandoc.table(
    df,
    justify = "left", 
    row.names = FALSE,
    style = "simple",
    split.table = Inf
  )
  
}

# Write output file if requested ----
if( args$output != FALSE ){
  if( stringr::str_detect(args$output, ".tsv$") ){
    write.table(df, file = args$output, quote = FALSE, row.names = FALSE)
  }else if( stringr::str_detect(args$output, ".csv$") ){
    write.table(df, file = args$output, quote = FALSE, row.names = FALSE)
  }else if( stringr::str_detect(args$output, ".rds$") ){
    saveRDS(df, file = args$output)
  }else if( stringr::str_detect(args$output, ".RData$") ){
    save(df, file = args$output)
  }
}

# Finish up and close out ----
if( all(df$outcome == "pass") ){
  q(save = "no", status = 0)
}else{
  q(save = "no", status = 1)
}

