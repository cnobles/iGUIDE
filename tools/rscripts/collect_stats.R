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

# Transform data into a wide format
wide_data <- dplyr::mutate(
    long_data, 
    type = paste0(type, ".", metric),
    type = factor(type, levels = unique(type))
  ) %>%
  dplyr::select(-metric) %>%
  dplyr::distinct() %>%
  tidyr::spread(type, count)

# Write data to output
write.csv(wide_data, file = args$output, quote = FALSE, row.names = FALSE)

q()