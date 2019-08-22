#' For anyone reviewing the code below, the following is a small style guide 
#' outlining the various formats for the code. 
#' 
#' Names with "_": objects, inlucding data.frames, GRanges, vectors, ...
#' Names in caMel format: functions or components of objects (i.e. columns 
#' within a data.frame).
#' Names with ".": arguments / options for functions

options(stringsAsFactors = FALSE, scipen = 99, width = 120)

# Set up and gather command line arguments ----
parser <- argparse::ArgumentParser(
  description = "Separate sequence files into bins of appropriate size.",
  usage = paste(
    "Rscript bin_seqs.R <seqs> -o <outputDir> [optional args] [-h/--help]"
  )
)

parser$add_argument(
  "seqs", nargs = "+", type = "character",
  help = paste(
    "Path(s) to sequence files to separate into bins. Only read names in first",
    "file will be used for indexing and splitting. Make sure all files have",
    "the same content! Read order will be set by first file. Fasta or Fastq",
    "formats allowed, as well as gzipped compression."
  )
)

parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", default = ".",
  help = "Directory for output files to be written. Default: '.'"
)

parser$add_argument(
  "-b", "--bins", nargs = 1, type = "integer", default = 2L,
  help = "The number of bins to separate files into, default is 2."
)

parser$add_argument(
  "-l", "--level", nargs = 1, type = "integer", default = 0L, 
  help = paste(
    "Fill level for each bin. If specified, then script will fill files to the",
    "specified level with reads before filling the next file, sequentially.",
    "If the total number of reads would fill all bins to their level, then",
    "reads will be evenly distributed across all bins, which is the default",
    "behavior. Default value: 0."
  )
)

parser$add_argument(
  "--compress", action = "store_true", 
  help = paste(
    "Output sequence file(s) in gzipped format. Otherwise this relies on the",
    "input format."
  )
)

parser$add_argument(
  "--readNamePattern", nargs = 1, type = "character", 
  default = "[\\w\\:\\-\\+]+", 
  help = paste(
    "Regular expression for pattern matching read names. Should not contain", 
    "R1/R2/I1/I2 specific components. Default is [\\w\\:\\-\\+]+"
  )
)


args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

# Create output directory if not currently available ----
if( !dir.exists(args$output) ){

  dir.create(args$output)
  if(!dir.exists(args$output)) stop("Cannot create output folder.\n")
  args$output <- normalizePath(args$output)
  
}

# Load sequence files
seq_files <- lapply(args$seqs, function(x){
  
  if( stringr::str_detect(x, ".fastq") | stringr::str_detect(x, ".fq") ){
    return(ShortRead::readFastq(x))
  }else{
    return(ShortRead::readFasta(x))
  }
  
})


# Score indices from first sequence for binning input sequences
if( length(seq_files[[1]]) <= args$bins * args$level ){

  seq_idx <- split(
      seq_along(seq_files[[1]]),
      ceiling(seq_along(seq_files[[1]]) / args$level)
  )
  
  if( length(seq_idx) < args$bins ){
    seq_idx <- c(
      seq_idx, 
      split(integer(), seq(length(seq_idx)+1, args$bins, 1))
    )
  }
  
}else{

  seq_idx <- split(
    seq_along(seq_files[[1]]), 
    ceiling(
      seq_along(seq_files[[1]])/(length(seq_files[[1]])/as.numeric(args$bins))
    )
  )
  
}

seq_names <- stringr::str_extract(
  as.character(ShortRead::id(seq_files[[1]])),
  args$readNamePattern
)

seq_name_list <- lapply(seq_idx, function(i) seq_names[i])

# Split and write sequences to output directory
output_files <- strsplit(args$seqs, "/")

output_files <- unlist(mapply(
  function(i, j) output_files[[i]][j],
  i = seq_along(output_files), 
  j = lengths(output_files),
  SIMPLIFY = FALSE
))

if( any(stringr::str_detect(output_files, ".gz$")) | args$compress ){
  args$compress <- TRUE
}else{
  args$compress <- FALSE
}

expanded_output_file_names <- lapply(output_files, function(x){
  
  x <- stringr::str_remove(x, ".gz$")
  
  ext <- unlist(strsplit(x, "\\."))
  lead <- paste(ext[-length(ext)], collapse = ".")
  ext <- ext[length(ext)]
  
  bins <- stringr::str_pad(seq_len(args$bins), nchar(args$bins), pad = 0)
  exp_names <- paste0(lead, ".bin", bins, ".", ext)
  
  if( args$compress ){
    exp_names <- paste0(exp_names, ".gz")
  }
  
  exp_names
  
})

# Write output files
null <- mapply(
  function(seqs, outputs, idx_names){
    
    null <- mapply(
      function(idx, outfile){
        
        matched_idx <- match(idx, stringr::str_extract(
          as.character(ShortRead::id(seqs)), args$readNamePattern
        ))
        
        if( any(table(matched_idx)) > 1 ){
          stop("\n  ReadNamePattern is ambiguous, please refine.")
        }
        
        if( file.exists(file.path(args$output, outfile)) ){
          unlink(file.path(args$output, outfile))
        }
  
        if( stringr::str_detect(outfile, ".fastq") | 
            stringr::str_detect(outfile, ".fq") ){
              
          ShortRead::writeFastq(
            object = seqs[matched_idx],
            file = file.path(args$output, outfile),
            compress = args$compress
          )
          
        }else{
          
          ShortRead::writeFasta(
            object = seqs[matched_idx],
            file = file.path(args$output, outfile),
            compress = args$compress
          )
          
        }
        
      },
      idx = idx_names,
      outfile = outputs
    )
    
  },
  seqs = seq_files,
  outputs = expanded_output_file_names,
  MoreArgs = list(idx_names = seq_name_list)
)

# Check for output files
if( 
  all(file.exists(file.path(args$output, unlist(expanded_output_file_names)))) 
){

  cat(
    "\nAll files written to output directory:\n ", 
    paste(
      file.path(args$output, unlist(expanded_output_file_names)), 
      collapse = "\n  "
    ),
    "\n"
  )
  
  q(save = "no", status = 0)
  
}else{
  
  stop("\n  Could not confirm existance of all output files.\n")
  
}

