# Consolidate and write out stats for iGUIDE processing
#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, scipen = 99)
suppressMessages(library("magrittr"))

# Capture commandline files
parser <- argparse::ArgumentParser(
  description = "Script to consolidate .stat files."
)

parser$add_argument(
  "-f", "--file", nargs = "+", type = "character", default = "NA",
  help = "Path to files with *.stat files (long, csv format). "
)

parser$add_argument(
  "-d", "--dir", nargs = "+", type = "character", default = "NA",
  help = "Path to directory with *.stat files (long, csv format). "
)

parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", 
  help = "Output file path and name, csv format. ie. path/to/file.csv"
)

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

stopifnot(! all(c(args$file, args$dir) == "NA") )

# Manipulate file paths to determine stat types
if( args$file != "NA" ){
  is_present <- file.exists(args$file)
}else{
  is_present <- file.exists(args$dir)
}

if( !is_present ){
  stop(
    "\n  Cannot find the following file(s) or directory: ", 
    c(args$file, args$dir)[c(args$file, args$dir) != "NA"]
  )
}

if( args$file != "NA"){
  file_names <- stringr::str_extract(args$file, "[\\w\\.\\-\\_]+$")
  file_paths <- args$file
}else{
  file_names <- list.files(path = args$dir, pattern = "\\.stat$")
  file_paths <- file.path(args$dir, file_names)
}

file_types <- sub("[\\w\\-\\_]+.", "", file_names, perl = TRUE)
file_types <- sub(".stat", "", file_types)

# Read in data in a long format
long_data <- dplyr::bind_rows(
  lapply(
    structure(file_paths, names = file_types), 
    function(file){
      
      x <- try(
        expr = read.csv(file = file, header = FALSE), 
        silent = TRUE
      )
      
      if( class(x) == "try-error"){

        return(data.frame(
          sampleName = vector(mode = "character"), 
          metric = vector(mode = "character"), 
          count = vector("numeric")
        ))
        
      }else{
        
        names(x) <- c("sampleName", "metric", "count")
      
        return(dplyr::mutate(
          x,
          sampleName = stringr::str_extract(sampleName, "[\\w\\-\\_]+")
        )) 
        
      }
      
    }
  ),
  .id = "type"
)

fmt_long_data <- long_data %>%
  dplyr::distinct(type, sampleName, metric, count) %>%
  dplyr::mutate(
    bin = stringr::str_extract(type, "bin[0-9]+"),
    read = ifelse(
      stringr::str_detect(type, "R[12]."),
      ifelse(stringr::str_detect(type, "R1."), "R1", "R2"),
      NA
    ),
    type = stringr::str_remove(type, "bin[0-9]+.")
  ) %>%
  dplyr::group_by(sampleName, type, metric, read) %>%
  dplyr::summarise(count = sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(
    (stringr::str_detect(metric, "multihit") & 
      stringr::str_detect(type, "multihits")) | 
      !stringr::str_detect(metric, "multihit")
  ) %>%
  dplyr::mutate(type = ifelse(type == "multihits", "align", type)) %>%
  dplyr::ungroup()

# Transform data into a wide format
wide_data <- dplyr::mutate(
    fmt_long_data, 
    type = paste0(type, ".", metric),
    type = factor(type, levels = unique(type))
  ) %>%
  dplyr::select(-metric, -read) %>%
  tidyr::spread(type, count)

wide_cols <- names(wide_data)

wide_data <- wide_data[
  ,c("sampleName", sort(wide_cols[-match("sampleName", wide_cols)]))
]

# Write data to output
write.csv(wide_data, file = args$output, quote = FALSE, row.names = FALSE)

q()