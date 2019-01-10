#' For anyone reviewing the code below, the following is a small style guide 
#' outlining the various formats for the code. 
#' 
#' Names with "_": objects, including data.frames, GRanges, vectors, ...
#' Names in caMel format: functions or components of objects (i.e. columns 
#' within a data.frame).
#' Names with ".": arguments / options for functions
# Set Global options and load intiial packages ----
options(stringsAsFactors = FALSE, scipen = 99, width = 999)

code_dir <- dirname(sub(
  pattern = "--file=", 
  replacement = "", 
  x = grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
))

desc <- yaml::yaml.load_file(
  file.path(code_dir, "descriptions/demulti.desc.yml")
)

# Set up and gather command line arguments ----
## Argument parser ----
parser <- argparse::ArgumentParser(
  description = desc$program_short_description,
  usage = "Rscript demulti.R [-h/--help, -v/--version] [optional args]"
)

parser$add_argument(
  "-m", "--manifest", type = "character", 
  help = desc$manifest
)

parser$add_argument(
  "--read1", type = "character", default = "NA", 
  help = desc$read1
)

parser$add_argument(
  "--read2", type = "character", default = "NA", 
  help = desc$read2
)

parser$add_argument(
  "--idx1", type = "character", default = "NA", 
  help = desc$idx1
)

parser$add_argument(
  "--idx2", type = "character", default = "NA", 
  help = desc$idx2
)

parser$add_argument(
  "-o", "--outfolder", nargs = 1, type = "character", 
  help = desc$outfolder
)

parser$add_argument(
  "--bc1", nargs = 1, type = "character", default = "I1", 
  help = desc$bc1
)

parser$add_argument(
  "--bc2", nargs = 1, type = "character", default = "I2", 
  help = desc$bc2
)

parser$add_argument(
  "--bc1Man", nargs = 1, type = "character", default = "barcode1", 
  help = desc$bc1Man
)

parser$add_argument(
  "--bc2Man", nargs = 1, type = "character", default = "barcode2", 
  help = desc$bc2Man
)

parser$add_argument(
  "--bc1Len", nargs = 1, type = "integer", default = 8, 
  help = desc$bc1Len
)

parser$add_argument(
  "--bc2Len", nargs = 1, type = "integer", default = 8,
  help = desc$bc2Len
)

parser$add_argument(
  "--maxMis", nargs = 1, type = "integer", 
  help = desc$maxMis
)

parser$add_argument(
  "--bc1Mis", nargs = 1, type = "integer", default = 0, 
  help = desc$bc1Mis
)

parser$add_argument(
  "--bc2Mis", nargs = 1, type = "integer", default = 0,
  help = desc$bc2Mis
)

parser$add_argument(
  "--stat", nargs = 1, type = "character", default = FALSE, 
  help = desc$stat
)

parser$add_argument(
  "-c", "--cores", nargs = 1, default = 1, type = "integer", 
  help = desc$cores
)

parser$add_argument(
  "--compress", action = "store_true", 
  help = desc$compress
)

parser$add_argument(
  "-p", "--poolreps", action = "store_true", 
  help = desc$poolreps
)

parser$add_argument(
  "--singleBarcode", action = "store_true", 
  help = desc$singleBarcode
)

parser$add_argument(
  "--readNamePattern", nargs = 1, type = "character", 
  default = "[\\w\\:\\-\\+]+", 
  help = desc$readNamePattern
)


args <- parser$parse_args(commandArgs(trailingOnly = TRUE))


demulti <- data.frame(
  "readType" = c("R1", "R2", "I1", "I2"),
  "path" = c(args$read1, args$read2, args$idx1, args$idx2)
)

demulti$bc1 <- grepl(args$bc1, demulti$readType)
demulti$bc2 <- grepl(args$bc2, demulti$readType)


if( demulti$readType[demulti$bc1] == demulti$readType[demulti$bc2] ){
  stop("\nPlease select different read types for barcodes 1 and 2.")
}

if( demulti$readType[demulti$bc1] == "NA" ){
  stop("\nBarcode 1 is set to a read type that is not provided.")
}

if( demulti$readType[demulti$bc2] == "NA" ){
  stop("\nBarcode 2 is set to a read type that is not provided.")
}

if( args$singleBarcode ){
  demulti$bc2 <- FALSE
}

if( !is.null(args$maxMis) ){
  args$bc1Mis <- args$maxMis
  args$bc2Mis <- args$maxMis
}


input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply( 
    seq_along(args), 
    function(i){
      paste(args[[i]], collapse = ", ")
    }
  )
)

input_table <- input_table[
  match(
    c("manifest :", "idx1 :", "idx2 :", "read1 :", "read2 :", 
      "outfolder :", "stat :", "bc1 :", "bc2 :", "bc1Man :", "bc2Man :", 
      "bc1Len :", "bc2Len :", "bc1Mis :", "bc2Mis :", "cores :", "compress :", 
      "poolreps :", "singleBarcode :",  "readNamePattern :"
    ),
    input_table$Variables
  ),
]

cat("\nDemultiplex Inputs:\n")
print(
  x = data.frame(input_table, row.names = NULL), 
  right = FALSE, 
  row.names = FALSE
)

# Create output directory if not currently available ----
if( !file.exists(args$outfolder) ){
  
  attempt <- try(system(paste0("mkdir ", args$outfolder)))
  if(attempt == 1) stop("Cannot create output folder.")

}

# Check for required packages ----
required_packs <- c("stringr", "ShortRead", "Biostrings")
present_packs <- required_packs %in% row.names(installed.packages())

if( !all(present_packs) ){
  
  cat("\nMissing required r-packages:\n")
  print(
    data.frame(
      "R-Packages" = required_packs, 
      "Installed" = present_packs, 
      row.names = NULL
    ), right = FALSE, row.names = FALSE)

  stop("\nCheck dependancies.")

}

# Operating functions ----
parseIndexReads <- function(barcode, index.file.path, barcode.length, 
                            max.mismatch, read.name.pattern){
  
  # Load index file sequences and sequence names
  index <- ShortRead::readFastq(index.file.path)
  index <- ShortRead::narrow(index, start = 1, end = barcode.length)
  index@id <- Biostrings::BStringSet(
    stringr::str_extract(
      as.character(ShortRead::id(index)), 
      read.name.pattern)
  )
  
  # Trim barcode if necessary
  barcode <- as.character(
    Biostrings::DNAStringSet(barcode, start = 1, end = barcode.length)
  )
  
  # Identify read names with sequences above or equal to the minscore
  vmp <- Biostrings::vmatchPattern(
    barcode, ShortRead::sread(index), max.mismatch = max.mismatch, fixed = FALSE
  )
  
  return( which(lengths(vmp) == 1) )
  
}

writeDemultiplexedSequences <- function(read.file.path, type, multiplexed.data, 
                                        read.name.pattern, out.folder, compress){
  
  # Load read sequences and sequence names then write to file
  reads <- ShortRead::readFastq(read.file.path)
  reads@id <- Biostrings::BStringSet(
    stringr::str_extract(as.character(ShortRead::id(reads)), read.name.pattern)
  )
  
  reads <- reads[multiplexed.data$index]
  reads <- split(reads, multiplexed.data$sampleName)
  
  null <- lapply(
    seq_along(reads), 
    function(i, reads, type, out.folder, compress){
      
      if(compress){  
        
        filePath <- file.path(
          out.folder, paste0(names(reads[i]), ".", type, ".fastq.gz")
        )
        
      }else{
        
        filePath <- file.path(
          out.folder, paste0(names(reads[i]), ".", type, ".fastq")
        )
        
      }
      
      if( file.exists(filePath) ) unlink(filePath)
      
      ShortRead::writeFastq(reads[[i]], file = filePath, compress = compress)
      
      cat(
        paste0("\nWrote ", length(reads[[i]]), " reads to:\n", filePath, "."))
      
    },
    reads = reads, type = type, out.folder = out.folder, compress = compress
  )
  
  return(list(read.file.path, type, out.folder))
  
}

# Load manifest / sample mapping file ----
file_ext <- unlist(strsplit(args$manifest, "\\."))
file_ext <- file_ext[length(file_ext)]

if( file_ext %in% c("yaml", "yml") ){
  
  if( !"yaml" %in% row.names(installed.packages()) ){
    stop("\nPackage:yaml not loaded or installed.")
  }
  
  manifest <- yaml::yaml.load_file(args$manifest)
  
  if( args$singleBarcode ){
    
    samples_df <- data.frame(
      "sampleName" = names(manifest$samples),
      "bc1" = sapply( manifest$samples, function(x) x[args$bc1Man] ),
      row.names = NULL
    )
    
  }else{
    
    samples_df <- data.frame(
      "sampleName" = names(manifest$samples),
      "bc1" = sapply( manifest$samples, function(x) x[args$bc1Man] ),
      "bc2" = sapply( manifest$samples, function(x) x[args$bc2Man] ),
      row.names = NULL
    )
    
  }

}else{
  
  if( file_ext == "csv" ){
    manifest <- read.csv(args$manifest)
  }else if( file_ext == "tsv" ){
    manifest <- read.delim(args$manifest)
  }
  
  if( args$singleBarcode ){
    samples_df <- manifest[, c("sampleName", args$bc1Man)]
    names(samples_df) <- c("sampleName", "bc1")
  }else{
    samples_df <- manifest[, c("sampleName", args$bc1Man, args$bc2Man)]
    names(samples_df) <- c("sampleName", "bc1", "bc2")
  }
  
}


if( !args$singleBarcode ){
  
  unique_samples <- nrow(samples_df[,c("bc1", "bc2")]) == 
    nrow(unique(samples_df[,c("bc1", "bc2")]))
  
  if( !unique_samples ){
    stop("\nAmbiguous barcoding of samples. Please correct.")
  }

}else{

  unique_samples <- length(samples_df[,c("bc1")]) == 
    length(unique(samples_df[,"bc1"]))
  
  if( !unique_samples ){
    stop("\nAmbiguous barcoding of samples. Please correct.")
  }

}

# Read in barcode sequences ----
bc1_reads <- ShortRead::readFastq(demulti$path[demulti$bc1])
cat(paste("\nReads to demultiplex : ", length(bc1_reads), "\n"))

if( args$cores > 1 ){
  
  cluster <- parallel::makeCluster(min(c(parallel::detectCores(), args$cores)))
  
  BC1_parsed <-  parallel::parLapply(
    cluster,
    unique(samples_df$bc1), 
    parseIndexReads,
    index.file.path = demulti$path[demulti$bc1],
    barcode.length = args$bc1Len,
    max.mismatch = args$bc1Mis,
    read.name.pattern = args$readNamePattern
  )
  
  names(BC1_parsed) <- unique(samples_df$bc1)
  
  cat("\nbc1 breakdown:\n")
  print(
    data.frame(
      "bc1" = names(BC1_parsed),
      "Read Counts" = sapply( BC1_parsed, length )
    ),
    right = TRUE, 
    row.names = FALSE
  )
  
  if( !args$singleBarcode ){
    
    BC2_parsed <- parallel::parLapply(
      cluster, 
      unique(samples_df$bc2), 
      parseIndexReads,
      index.file.path = demulti$path[demulti$bc2],
      barcode.length = args$bc2Len,
      max.mismatch = args$bc2Mis,
      read.name.pattern = args$readNamePattern
    )
    
  }
  
  parallel::stopCluster(cluster)
  
}else{
  
  BC1_parsed <-  lapply(
    unique(samples_df$bc1), 
    parseIndexReads,
    index.file.path = demulti$path[demulti$bc1],
    barcode.length = args$bc1Len,
    max.mismatch = args$bc1Mis,
    read.name.pattern = args$readNamePattern
  )

  names(BC1_parsed) <- unique(samples_df$bc1)
  
  cat("\nbc1 breakdown:\n")
  print(
    data.frame(
      "bc1" = names(BC1_parsed),
      "Read Counts" = sapply( BC1_parsed, length )
    ),
    right = TRUE, 
    row.names = FALSE
  )
  
  if( !args$singleBarcode ){
    
    BC2_parsed <- lapply(
      unique(samples_df$bc2), 
      parseIndexReads,
      index.file.path = demulti$path[demulti$bc2],
      barcode.length = args$bc2Len,
      max.mismatch = args$bc2Mis,
      read.name.pattern = args$readNamePattern
    )
    
  }
  
}

if( !args$singleBarcode ){

  names(BC2_parsed) <- unique(samples_df$bc2)
  cat("\nbc2 breakdown:\n")
  print(
    data.frame(
      "bc2" = names(BC2_parsed),
      "Read Counts" = sapply(BC2_parsed, length)
    ),
    right = TRUE,
    row.names = FALSE
  )

}

if( !args$singleBarcode ){
  
  demultiplexed_indices <- mapply(
    function(bc1, bc2){
      Biostrings::intersect(BC1_parsed[[bc1]], BC2_parsed[[bc2]])
    },
    bc1 = samples_df$bc1,
    bc2 = samples_df$bc2,
    SIMPLIFY = FALSE
  )
  
  names(demultiplexed_indices) <- paste0(
    samples_df$bc1, samples_df$bc2
  )
  
}else{
  
  demultiplexed_indices <- BC1_parsed
  
}

# As there is some flexibility in the barcode matching, some reads may be 
# be assigned to multiple samples. These reads are ambiguous and will be 
# removed.
ambiguous_indices <- unique(
  unlist(demultiplexed_indices)[duplicated(unlist(demultiplexed_indices))]
)

demultiplexed_indices <- lapply(demultiplexed_indices, function(x, reads){
    x[!x %in% reads]
  },
  reads = ambiguous_indices
)

# Reads by sample
samples_df$read_counts <- sapply( demultiplexed_indices, length )
cat("\nRead counts for each sample.\n")
print(samples_df, split.tables = Inf)

# Ambiguous reads
cat(paste0("\nAmbiguous reads: ", length(ambiguous_indices), "\n"))

# Unassigned reads
all_indices <- seq_along(bc1_reads)

unassigned_indices <- all_indices[
  !all_indices %in% unlist(demultiplexed_indices, use.names = FALSE)
]

unassigned_indices <- unassigned_indices[
  !unassigned_indices %in% ambiguous_indices
]

print(paste0("\nUnassigned reads: ", length(unassigned_indices), "\n"))

if( args$stat != FALSE ){
  write.table(
    data.frame(
      sampleName = paste0(
        c(samples_df$sampleName, "ambiguous_reads", "unassigned_reads"), 
        ".demulti"
      ),
      metric = "reads",
      count = c(
        samples_df$read_counts, 
        length(ambiguous_indices), 
        length(unassigned_indices)
      )
    ),
    file = file.path(args$outfolder, args$stat),
    sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE
  )
}

# Create multiplex dataframe for subseting sequencing files ----
multiplexed_data <- data.frame(
  "sampleName" = S4Vectors::Rle(
    values = samples_df$sampleName, 
    length = sapply(demultiplexed_indices, length)
  ),
  "index" = unlist(demultiplexed_indices),
  row.names = NULL
)

ambiguous_data <- data.frame(
  "sampleName" = rep("ambiguous", length(ambiguous_indices)),
  "index" = ambiguous_indices,
  row.names = NULL
)

unassignedData <- data.frame(
  "sampleName" = rep("unassigned", length(unassigned_indices)),
  "index" = unassigned_indices,
  row.names = NULL
)

multiplexed_data <- rbind(multiplexed_data, ambiguous_data, unassignedData)
multiplexed_data$sampleName <- factor(
  multiplexed_data$sampleName,
  levels = c(samples_df$sampleName, "ambiguous", "unassigned")
)

stopifnot(nrow(multiplexed_data) == length(all_indices))

if( args$poolreps ){
  multiplexed_data$sampleName <- gsub("-\\d+$", "", multiplexed_data$sampleName)
}

cat(paste0("\nReads to be written to files: ", nrow(multiplexed_data), "\n"))

# Write files to read files to outfolder directory
if( args$cores > 1 ){
  
  cluster <- parallel::makeCluster(min(c(parallel::detectCores(), args$cores)))
  
  read_list <- demulti$readType[demulti$path != "NA"]
  read_paths <- demulti$path[match(read_list, demulti$readType)]

  demultiplex <- parallel::clusterMap(
    cluster,
    writeDemultiplexedSequences,
    read.file.path = read_paths,
    type = read_list,
    MoreArgs = list(
      multiplexed.data = multiplexed_data,
      read.name.pattern = args$readNamePattern,
      out.folder = args$outfolder,
      compress = args$compress
    )
  )
  
  parallel::stopCluster(cluster)
  
}else{
  
  read_list <- demulti$readType[demulti$path != "NA"]
  read_paths <- demulti$path[match(read_list, demulti$readType)]
  
  demultiplex <- mapply(
    writeDemultiplexedSequences,
    read.file.path = read_paths,
    type = read_list,
    MoreArgs = list(
      multiplexed.data = multiplexed_data,
      read.name.pattern = args$readNamePattern,
      out.folder = args$outfolder,
      compress = args$compress
    )
  )

}

cat("\nDemultiplexing complete.\n")
q()
