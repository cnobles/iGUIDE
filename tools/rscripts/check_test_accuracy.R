# Import yaml file with test paths and md5 digests
# Required libraries stringr, digest, Biostrings, data.table, yaml, pander
options(stringsAsFactors = FALSE, scipen = 99, warn = -1, window = 999)
suppressMessages(library("magrittr"))

# Set up and gather command line arguments ----
parser <- argparse::ArgumentParser(
  description = "Check accuracy in processing test data set.",
  usage = paste(
    "Rscript tools/rscripts/check_test_accuracy.R <run.config> <test.truth>",
    "<options>"
  )
)

parser$add_argument(
  "run.config", nargs = 1, type = "character",
  help = paste(
    "Yaml config file used to process the run.", 
    "i.e. simulation.config.yml"
  )
)

parser$add_argument(
  "test.truth", nargs = 1, type = "character",
  help = paste(
    "CSV file with the original 'true' data for testing accuracy.",
    "i.e. truth.csv"
  )
)

parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", default = FALSE,
  help = "Output file name in an .rds format."
)

parser$add_argument(
  "-v", "--verbose", action = "store_true", 
  help = "Turns on diagnositc-based messages."
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
    c(
      "run.config :", "test.truth :", "output :", "verbose :", "iguide_dir :"
    ), 
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
  supported_ext <- c("tsv", "csv", "gz", "fasta", "fastq", "rds", "yaml", "yml")
  
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
    
    if( ext == "gz" ){
      return(read.table(gzfile(path), header = TRUE, sep = ","))
    }else{
      return(read.table(path, header = TRUE, sep = ","))
    }
    
  }else if( any(stringr::str_detect(exts, "fast")) ){
    
    return(Biostrings::readDNAStringSet(path))
    
  }else if( any(exts %in% c("yaml", "yml")) ){
    
    return(yaml::yaml.load_file(path))
    
  }else if( any(exts %in% c("rds")) ){
    
    rds_import <- readRDS(path)

    if( class(rds_import) == "list" ){

      if( any(sapply(rds_import, class) == "data.frame") ){

        idx <- which(sapply(rds_import, class) == "data.frame")
        return(rds_import[[idx[1]]])

      }else{
        return(as.data.frame(rds_import[[1]], row.names = NULL))
      }

    }else{
      return(rds_import)
    }
    
  }else{
    
    stop("\n  Unsupported input file time.\n")
    
  }
  
}

# Load inputs ----
run_config <- readFile(args$run.config, args$iguide_dir)
test_truth <- readFile(args$test.truth, args$iguide_dir)
sample_info <- readFile(run_config$Sample_Info, args$iguide_dir)

# Files to check ----
check_files <- paste0(
  "analysis/", run_config$Run_Name, "/output/incorp_sites.", 
  run_config$Run_Name, ".rds"
)

check_data <- lapply(check_files, readFile, root = args$iguide_dir)
names(check_data) <- c("uniq_sites")

check_data$multihits <- suppressMessages(dplyr::bind_rows(lapply(
  sample_info$sampleName, 
  function(x){

    readFile(
      paste0(
        "analysis/", run_config$Run_Name, "/process_data/multihits/", 
        x, ".multihits.rds"
      ), 
      args$iguide_dir
    )

  }), 
  .id = "specimen"
))

## Check for content ----
total_reads <- length(test_truth$read.name)

total_read_ids <- split(
  test_truth$read.name, 
  stringr::str_extract(test_truth$read.name, "[\\w]+\\-[\\w]+\\-[\\w]+")
)[
  unique(stringr::str_extract(test_truth$read.name, "[\\w]+\\-[\\w]+\\-[\\w]+"))
]

collected_stats <- dplyr::bind_rows(lapply(total_read_ids, function(x){
    
      ret <- c(
        "uniq" = sum(x %in% check_data$uniq_sites$ID), 
        "multi" = sum(x %in% check_data$multihits$ID), 
        "comb" = sum(x %in% check_data$uniq_sites$ID) + 
          sum(x %in% check_data$multihits$ID),
        "total" = length(x)
      )
      
      x <- x[x %in% check_data$uniq_sites$ID]
      spec_truth <- test_truth[match(x, test_truth$read.name),]
      uniq_sites <- check_data$uniq_sites[match(x, check_data$uniq_sites$ID),]
      seq_cnt <- sum(spec_truth$seqnames == uniq_sites$seqnames)
      std_cnt <- sum(spec_truth$strand == uniq_sites$strand)
      cum_dis <- sum(abs(spec_truth$start - uniq_sites$start)) + 
        sum(abs(spec_truth$end - uniq_sites$end))
      correct <- sum(
        spec_truth$seqnames == uniq_sites$seqnames & 
          spec_truth$strand == uniq_sites$strand &
          abs(spec_truth$start - uniq_sites$start) + 
          abs(spec_truth$end - uniq_sites$end) == 0
      )
      
      ret <- c(
        ret, 
        c(
          "seqs" = seq_cnt, "strand" = std_cnt, "dis" = cum_dis, 
          "cor" = correct
        )
      )
      
      data.frame(t(ret))
      
    }),
    .id = "type"
  ) %>%
  tidyr::separate(type, into = c("specimen", "target", "gRNA"), sep = "-")

missing_data <- dplyr::bind_rows(lapply(total_read_ids, function(x){
      
      x <- x[!x %in% check_data$uniq_sites$ID]
      x <- x[!x %in% check_data$multihits$ID]
      spec_truth <- test_truth[match(x, test_truth$read.name),]
      
      dist <- abs(
        as.numeric(stringr::str_extract(spec_truth$posid, "[0-9]+$")) - 
          as.numeric(stringr::str_extract(spec_truth$incorp, "[0-9]+$")))
      
      ret <- c(
        "count" = nrow(spec_truth), 
        "min_dist" = min(dist), "max_dist" = max(dist), 
        "min_width" = min(spec_truth$width), "max_width" = max(spec_truth$width)
      )
      
      data.frame(t(ret))
      
    }),
    .id = "type"
  ) %>%
  tidyr::separate(type, into = c("specimen", "target", "gRNA"), sep = "-") %>%
  dplyr::mutate(
    min_dist = ifelse(count == 0, 0, min_dist),
    max_dist = ifelse(count == 0, 0, max_dist),
    min_width = ifelse(count == 0, 0, min_width),
    max_width = ifelse(count == 0, 0, max_width)
  )

pct_retention <- 100 * sum(collected_stats$comb) / sum(collected_stats$total)
uniq_accuracy <- 100 * sum(collected_stats$cor) / sum(collected_stats$uniq)

# Log output if requested ----
if( args$verbose ){
  
  cat("\nCollected Stats:")
  pander::pandoc.table(
    collected_stats,
    justify = "left", 
    row.names = FALSE,
    style = "simple",
    split.table = Inf
  )

  cat("\nRead retention:", round(pct_retention, digits = 1), "%\n")
  cat("Unique accuracy:", round(uniq_accuracy, digits = 1), "%\n")
  
  cat("\nMissing data:")
  pander::pandoc.table(
    missing_data,
    justify = "left", 
    row.names = FALSE,
    style = "simple",
    split.table = Inf
  )
  
}

# Write output file if requested ----
if( args$output != FALSE ){
  if( stringr::str_detect(args$output, ".rds$") ){
    saveRDS(
      list(
        "collected_stats" = collected_stats,
        "missing_data" = missing_data,
        "test_truth" = test_truth,
        "checked_data" = check_data
      ), 
      file = args$output
    )
  }else{
    stop("\n  Output data object must be a .rds format.\n")
  }
}

# Finish up and close out ----
if( pct_retention >= 95 & uniq_accuracy >= 99 ){
  q(save = "no", status = 0)
}else{
  q(save = "no", status = 1)
}
