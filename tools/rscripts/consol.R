#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, scipen = 99, width = 999)

code_dir <- dirname(sub(
  pattern = "--file=", 
  replacement = "", 
  x = grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
))

desc <- desc <- yaml::yaml.load_file(
  file.path(code_dir, "descriptions/consol.desc.yml")
)

#' Set up and gather command line arguments
parser <- argparse::ArgumentParser(
  description = desc$program_short_description,
  usage = "nuc consol <seqFile> [-h/--help, -v/--version] [optional args]"
)

parser$add_argument(
  "seqFile", nargs = 1, type = "character", default = "NA",
  help = desc$seqFile
)

parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", default = "NA",
  help = desc$output
)

parser$add_argument(
  "-k", "--keyFile", nargs = 1, type = "character", default = "NA",
  help = desc$keyFile
)

parser$add_argument(
  "-l", "--seqName", nargs = 1, type = "character", default = "NA",
  help = desc$seqName
)

parser$add_argument(
  "--stat", nargs = 1, type = "character", default = FALSE, help = desc$stat
)

parser$add_argument(
  "--compress", action = "store_true", help = desc$compress
) 

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

if( args$seqFile == "NA" ){
  stop("\n  No sequence file specified. Please provide.\n")
}

# Check I/O file types
seq_type <- unlist(strsplit(args$seqFile, "/"))
seq_type <- seq_type[length(seq_type)]
seq_type <- stringr::str_extract(seq_type, "fa[\\w]*")

if( !seq_type %in% c("fa", "fasta", "fastq") ){
  stop(desc$unrecognized_file_type, " ", desc$compression_note)
}

seq_type <- ifelse(seq_type %in% c("fa", "fasta"), "fasta", "fastq")

if( args$output != "NA" ){
  
  outType <- unlist(strsplit(args$output, "/"))
  outType <- outType[length(outType)]
  outType <- stringr::str_extract(outType, "fa[\\w]*")
  args$output <- unlist(strsplit(args$output, outType))[1]
  
  if( !outType %in% c("fa", "fasta", "fastq") ){
    stop(desc$unrecognized_file_type)
  }
  
  outType <- ifelse(outType %in% c("fa", "fasta"), "fasta", "fastq")
  
  if( outType == "fastq" ){
    message(desc$fastq_input)
    outType <- "fasta"
  }
  
  args$output <- paste0(args$output, outType)
  
}else{
  
  stop("\n  No output file name given. Please provide.\n")
  
}

if( args$keyFile != "NA" ){
  
  key_type <- stringr::str_extract(args$keyFile, "[\\w]+$")
  
  if( !key_type %in% c("csv", "tsv", "rds", "RData") ){
    stop(desc$output_keyfile_type)
  }
  
}else{
  
  stop("\n  No key file name given. Please provide.\n")
  
}
  
# Check sequence name lead
if( args$seqName == "NA" ){
  
  parsedName <- unlist(strsplit(args$seqFile, "/"))[
    length(unlist(strsplit(args$seqFile, "/")))
  ]
  
  args$seqName <- unlist(strsplit(parsedName, "fa[\\w]*"))[1]
  
}

# Print inputs to table
input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(
    seq_along(args), function(i){
      paste(args[[i]], collapse = ", ")
    }
  )
)

input_table <- input_table[
  match(
    c("seqFile :", "output :", "keyFile :", "seqName :", "stat :"),
    input_table$Variables)
,]

cat("\nConsolidate Inputs:\n")
print(
  data.frame(input_table, row.names = NULL), 
  right = FALSE, row.names = FALSE
)


# Read sequence file
if( seq_type == "fasta" ){
  seq_pointer <- ShortRead::readFasta(args$seqFile)
}else{
  seq_pointer <- ShortRead::readFastq(args$seqFile)
}

seqs <- ShortRead::sread(seq_pointer)
names(seqs) <- ShortRead::id(seq_pointer)

# Generate blank files if inputs are empty
if( length(seqs) == 0 ){
  
  Biostrings::writeXStringSet(
    x = Biostrings::DNAStringSet(),
    filepath = args$output,
    format = "fasta",
    compress = args$compress
  )
  
  if( !is.null(args$keyFile) ){
    
    key <- data.frame("readNames" = c(), "seqID" = c())
    
    if( key_type == "csv" ){
      write.csv(key, file = args$keyFile, row.names = FALSE, quote = FALSE)
    }else if( key_type == "tsv" ){
      write.table(
        key, file = args$keyFile, sep = "\t", row.names = FALSE, quote = FALSE
      )
    }else if(key_type == "rds"){
      saveRDS(key, file = args$keyFile)
    }else if(key_type == "RData"){
      save(key, file = args$keyFile)
    }
    
  }
  
  
  if( args$stat != FALSE ){
    
    sampleName <- unlist(strsplit(args$output, "/"))
    
    sampleName <- unlist(
      strsplit(sampleName[length(sampleName)], ".fa", fixed = TRUE)
    )[1]
    
    write.table(
      x = data.frame(
        sampleName = sampleName,
        metric = "reads",
        count = length(seqs)
      ),
      file = args$stat,
      sep = ",", 
      row.names = FALSE, 
      col.names = FALSE, 
      quote = FALSE
    )
    
  }
  
  q()
  
}
  

factor_seqs <- factor(as.character(seqs))

key <- data.frame(
  "readNames" = names(factor_seqs),
  "seqID" = paste0(args$seqName, as.integer(factor_seqs))
)

consolidated_seqs <- Biostrings::DNAStringSet(levels(factor_seqs))
names(consolidated_seqs) <- paste0(args$seqName, seq_along(levels(factor_seqs)))

if( !is.null(args$output) ){
  
  if( args$compress & !grepl(".gz", args$output) ){
    args$output <- paste0(args$output, ".gz")
  }
  
}

# Stats if requested
if( args$stat != FALSE ){
  
  sampleName <- unlist(strsplit(args$output, "/"))
  sampleName <- unlist(
    strsplit(sampleName[length(sampleName)], ".fa", fixed = TRUE)
  )[1]
  
  write.table(
    x = data.frame(
      sampleName = sampleName,
      metric = "reads",
      count = length(consolidated_seqs)),
    file = args$stat,
    sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE
  )
  
}

# Write output and key files
# Output
if( is.null(args$output) ){
  
  print(
    data.frame(
      "seqID" = names(consolidated_seqs),
      "sequence" = as.character(consolidated_seqs)
    ),
    row.names = FALSE
  )
  
}else{
  
  unlink(args$output)
  
  ShortRead::writeFasta(
    consolidated_seqs, 
    file = args$output,
    width = max(Biostrings::width(consolidated_seqs)),
    compress = args$compress
  )
  
}

# Key file
if( !is.null(args$keyFile) ){
  
  if( key_type == "csv" ){
    write.csv(key, file = args$keyFile, row.names = FALSE, quote = FALSE)
  }else if( key_type == "tsv" ){
    write.table(
      key, file = args$keyFile, sep = "\t", row.names = FALSE, quote = FALSE
    )
  }else if( key_type == "rds" ){
    saveRDS(key, file = args$keyFile)
  }else if( key_type == "RData" ){
    save(key, file = args$keyFile)
  }
  
}

q()
