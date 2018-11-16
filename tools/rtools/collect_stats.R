# Consolidate and write out stats for iGUIDE processing
#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, scipen = 99)
suppressMessages(library("argparse"))
suppressMessages(library("magrittr"))

# Capture commandline files
parser <- ArgumentParser(description = "Script to consolidate .stat files.")
parser$add_argument(
  "paths", nargs = "+", type = "character", help = "List of *.stat file paths.")
parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", 
  help = "Output file path and name, csv format. ie. path/to/file.csv")

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

# Manipulate file paths to determine stat types
all_paths <- args$paths
file_names <- stringr::str_extract(all_paths, "[\\w\\.\\-\\_]+$")
file_types <- stringr::str_remove(file_names, "[\\w\\-\\_]+.") %>%
  stringr::str_remove(".stat")

# Read in data in a long format
long_data <- dplyr::bind_rows(
  lapply(structure(all_paths, names = file_types), function(path){
    x <- read.csv(path, header = FALSE)
    names(x) <- c("sampleName", "metric", "count")
    x %>%
      dplyr::mutate(
        sampleName = stringr::str_extract(sampleName, "[\\w\\-\\_]+")) }),
  .id = "type")

# Transform data into a wide format
wide_data <- dplyr::mutate(
    long_data, 
    type = paste0(type, ".", metric),
    type = factor(type, levels = unique(type))) %>%
  dplyr::select(-metric) %>%
  dplyr::distinct() %>%
  tidyr::spread(type, count)

# Write data to output
write.csv(wide_data, file = args$output, quote = FALSE, row.names = FALSE)

q()