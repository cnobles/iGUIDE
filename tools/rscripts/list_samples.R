#' list_samples.R
#' 
#' usage: Rscript list_samples.R <config>
#' 
#' This script lists the samples and their supplementary data (if provided) to 
#' the consol. This can be a useful feature if you've systemtically named run 
#' directories and would like to know which of the processed (or unprocessed) 
#' directories contains a specific set of samples.
#' 
#' For anyone reviewing the code below, the following is a small style guide 
#' outlining the various formats for the code. 
#' 
#' Names with "_": objects, inlucding data.frames, GRanges, vectors, ...
#' Names in caMel format: functions or components of objects (i.e. columns 
#' within a data.frame).
#' Names with ".": arguments / options for functions

# Required / highly suggested option parameters and library ----
options(stringsAsFactors = FALSE, scipen = 99, warn = -1)
suppressMessages(library("magrittr"))


# Set up and gather command line arguments ----
parser <- argparse::ArgumentParser(
  description = "List samples associated with a config file for iGUIDE.",
  usage = "iguide list_samples <path/to/config.file> <options>"
)

parser$add_argument(
  "config", nargs = 1, type = "character",
  help = "Run specific config file in yaml format."
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

## Set arguments with parser
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

root_dir <- Sys.getenv("IGUIDE_DIR")

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
    c("config :", "output :", "verbose :", "install_path :"), 
    input_table$Variables
  ),
]

## Log inputs
if( args$verbose ){
  
  cat("List Sample Inputs\n")
  print(
    x = data.frame(input_table),
    right = FALSE, 
    row.names = FALSE
  )

}


# Load files ----
## Config
if( file.exists(args$config) ){
  config <- yaml::yaml.load_file(args$config)
}else{
  stop("\nCannot find config file: ", args$config, ".\n")
}

## Sample Info
if( file.exists(config$Sample_Info) ){
  
  sample_info <- data.table::fread(config$Sample_Info, data.table = FALSE)
  
}else if( file.exists(file.path(root_dir, config$Sample_Info)) ){
  
  sample_info <- data.table::fread(
    input = file.path(root_dir, config$Sample_Info), 
    data.table = FALSE
  )
  
}else{
  
  stop("\nCannot find associated Sample Info file: ", configs$Sample_Info, ".\n")
  
}

## Supplemental Info
if( config$Supplemental_Info != "." ){
  if( file.exists(config$Supplemental_Info) ){
    
    supp_info <- data.table::fread(config$Supplemental_Info)
    
  }else if( file.exists(file.path(root_dir, config$Supplemental_Info)) ){
    
    supp_info <- data.table::fread(
      input = file.path(root_dir, config$Supplemental_Info), 
      data.table = FALSE
    )
    
  }else{
    
    warning(
      "Cannot find Supplemental Info file: ", configs$Supplemental_Info, ".\n"
    )
    
  }
}

# Join appropriate tables together and / or format for output ----

sample_col <- match(config$Sample_Name_Column, names(sample_info))

if( is.na(sample_col) ){
  stop("\nCannot isolate sampleName column: ", config$Sample_Name_Column, ".\n")
}

names(sample_info)[sample_col] <- "sampleName"

sample_info <- sample_info %>%
  dplyr::mutate(
    specimen = stringr::str_extract(sample_info$sampleName, "[\\w]+")
  ) %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(replicates = n()) %>%
  dplyr::ungroup()

if( exists("supp_info") ){
  sample_info <- dplyr::left_join(sample_info, supp_info, by = "specimen")
}


# Output consolidated information ----
if( args$output != FALSE ){
  
  source(file.path(code_dir, "supporting_scripts/writeOutputFile.R"))
  writeOutputFile(as.data.frame(sample_info), args$output)
  
}else{
  
  run_name <- stringr::str_extract(args$config, "[\\w]+.config.yml$") %>%
    stringr::str_extract("[\\w]+")
  
  cat(paste0("\nSpecimen Info for : ", run_name, "."))
  pander::pandoc.table(sample_info, style = "simple", split.table = Inf)
  
}

q()
