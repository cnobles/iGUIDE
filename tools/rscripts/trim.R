#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, scipen = 99, width = 999)

code_dir <- dirname(sub(
  pattern = "--file=", 
  replacement = "", 
  x = grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
))

desc <- yaml::yaml.load_file(
  file.path(code_dir, "descriptions/trim.desc.yml")
)

#' Set up and gather command line arguments
parser <- argparse::ArgumentParser(
  description = desc$program_short_description,
  usage = "nuc trim <seqFile> [-h/--help, -v/--version] [optional args]"
)

parser$add_argument(
  "seqFile", nargs = 1, type = "character", help = desc$seqFile
)

parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "-l", "--leadTrimSeq", nargs = 1, type = "character", default = ".",
  help = desc$leadTrimSeq
)

parser$add_argument(
  "-r", "--overTrimSeq", nargs = 1, type = "character", default = ".",
  help = desc$overTrimSeq
)

parser$add_argument(
  "--phasing", nargs = 1, type = "integer", default = 0, 
  help = desc$phasing
)

parser$add_argument(
  "--maxMismatch", nargs = 1, type = "integer", help = desc$maxMismatch
)

parser$add_argument(
  "--leadMismatch", nargs = "+", type = "integer", default = 0,
  help = desc$leadMismatch
)

parser$add_argument(
  "--overMismatch", nargs = 1, type = "integer", default = 0,
  help = desc$overMismatch
)

parser$add_argument(
  "--overMaxLength", nargs = 1, type = "integer", default = 20,
  help = desc$overMaxLength
)

parser$add_argument(
  "--overMinLength", nargs = 1, type = "integer", default = 3,
  help = desc$overMinLength
)

parser$add_argument(
  "--minSeqLength", nargs = 1, type = "integer", default = 30,
  help = desc$minSeqLength
)

parser$add_argument(
  "--collectRandomIDs", nargs = "+", type = "character", default = FALSE,
  help = desc$collectRandomIDs
)

parser$add_argument(
  "--badQualBases", nargs = 1, type = "integer", default = 5,
  help = desc$basQualBases
)

parser$add_argument(
  "--qualSlidingWindow", nargs = 1, type = "integer", default = 10,
  help = desc$qualSlidingWindow
)

parser$add_argument(
  "--qualThreshold", nargs = 1, type = "character", default = '?',
  help = desc$qualThreshold
)

parser$add_argument(
  "--stat", nargs = 1, type = "character", default = FALSE, help = desc$stat
)

parser$add_argument(
  "-c", "--cores", nargs = 1, default = 1, type = "integer", help = desc$cores
)

parser$add_argument(
  "--compress", action = "store_true", help = desc$compress
)

parser$add_argument(
  "--noFiltering", action = "store_true",
  help = desc$noFiltering
)

parser$add_argument(
  "--noQualTrimming", action = "store_true",
  help = desc$noQualTrimming
)



args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

if( is.null(args$seqFile) ){
  stop("\n  Please choose a sequence file (fasta or fastq).\n")
}

if( !is.null(args$maxMismatch) ){
  args$leadMismatch <- args$maxMismatch
  args$overMismatch <- args$maxMismatch
}

if( args$overMaxLength == 0 ){
  args$overMaxLength <- nchar(args$overTrimSeq)
}

if( all(args$collectRandomIDs != FALSE) ){
  if( !grepl("N", args$leadTrimSeq) ){
    cat(
      "\n  No random nucleotides (Ns) found in leadTrimSeq.",
      "Turning off collection of randomIDs.\n"
    )
    args$collectRandomIDs <- FALSE
  }
}

if( args$leadTrimSeq == "." ){
  args$leadTrimSeq <- ""
}

if( args$overTrimSeq == "." ){
  args$overTrimSeq <- ""
}

if( args$cores <= 0 ){
  args$cores <- 1
}

input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(
    seq_along(args), 
    function(i) paste(args[[i]], collapse = ", ")
  )
)

input_table <- input_table[
  match(
    c("seqFile :", "output :", "leadTrimSeq :", "overTrimSeq :", 
      "phasing :", "maxMismatch :", "leadMismatch :", "overMismatch :", 
      "overMaxLength :", "overMinLength :", "minSeqLength :", 
      "collectRandomIDs :", "noFiltering :", "noQualTrimming :", 
      "badQualBases :", "qualSlidingWindow :", "qualThreshold :", 
      "stat :", "compress :", "cores :"),
    input_table$Variables)
  ,]

cat("\nTrim Inputs:\n")
print(
  data.frame(input_table, row.names = NULL), 
  right = FALSE, 
  row.names = FALSE
)

# Reduce number of requested cores if needed.
if( args$cores > 1 ){
  if( args$cores > parallel::detectCores() ){
    cat(
      "\n  Requested cores is greater than availible for system.",
      "Changing cores to max allowed."
    )
    args$cores <- detectCores()
  }
}
  
# Load supporting scripts
source(file.path(code_dir, "supporting_scripts", "trimLeading.R"))

source(file.path(code_dir, "supporting_scripts", "trimOverreading.R"))

source(file.path(code_dir, "supporting_scripts", "writeSeqFiles.R"))

source(file.path(code_dir, "supporting_scripts", "utility_funcs.R"))

if( !all(
  c("trimLeading", "trimOverreading", "writeSeqFiles", "logSeqData", 
    "serialAppendS4") %in% 
  ls()
)){
  stop(
    "\n  Cannot load supporting scripts. ",
    "You may need to clone from github again.\n"
  )
}

# Determine sequence file types
seq_type <- seqFileType(args$seqFile)
out_type <- seqFileType(args$output)

# Determine random output file type
if( all(args$collectRandomIDs != FALSE) ){
  random_type <- seqFileType(args$collectRandomIDs)
}

# Read sequence file
if( seq_type == "fasta" ){
  seqs <- ShortRead::readFasta(args$seqFile)
}else{
  seqs <- ShortRead::readFastq(args$seqFile)
}

# Log info
input_tbl <- logSeqData(seqs)
cat("\nInput sequence information:\n")
print(input_tbl, row.names = FALSE)

# If no reads remaining, terminate and write output
if( length(seqs) == 0 ){
  
  cat(
    "\n  No reads remaining to trim. Terminating script after writing output.\n"
  )
  
  writeNullFile(
    file = args$output, 
    writeRandom = args$collectRandomIDs, 
    stat = args$stat,
    compress = args$compress
  )
  
  q()
  
}

# Quality trimming, trim from left to remove consecutive bad quality bases.
## Below block sets the OpenMP threads to the cores specified in args.
if( !args$noQualTrimming & seq_type == "fastq" ){
  
  nthreads <- .Call(ShortRead:::.set_omp_threads, as.integer(args$cores))
  on.exit(.Call(ShortRead:::.set_omp_threads, nthreads))

  seqs <- ShortRead::trimTailw(
    object = seqs, 
    k = args$badQualBases, 
    a = args$qualThreshold, 
    halfwidth = round(args$qualSlidingWindow/2)
  )
  
  # Log info
  qual_trimmed_tbl <- logSeqData(seqs)
  cat("\nSequence information remaining after quality trimming:\n")
  print(qual_trimmed_tbl, row.names = FALSE)
  
}

# If no reads remaining, terminate and write output
if( length(seqs) == 0 ){

  cat(
    "\n  No reads remaining to trim. Terminating script after writing output.\n"
  )
  
  writeNullFile(
    file = args$output, 
    writeRandom = args$collectRandomIDs, 
    stat = args$stat,
    compress = args$compress
  )
  
  q()
  
}

# Remove sequences that do not contain enough sequence information
seqs <- seqs[
  Biostrings::width(seqs) >= (
    args$minSeqLength + nchar(args$leadTrimSeq) + args$phasing
  )
]

len_trimmed_tbl <- logSeqData(seqs)
cat("\nSequence information remaining after minimum length trimming:\n")
print(len_trimmed_tbl, row.names = FALSE)

# Trim sequences, either on a single core or multiple cores
if( args$cores <= 1 ){
  
  # Trim 5' end or leading end. Conditionals present for added features.
  if( nchar(args$leadTrimSeq) > 0 ){

    trimmed_seqs <- trimLeading(
      seqs,
      trim.sequence = args$leadTrimSeq,
      phasing = args$phasing,
      max.mismatch = args$leadMismatch,
      collect.random = all(args$collectRandomIDs != FALSE),
      filter = !args$noFiltering
    )
    
  }else{
    
    trimmed_seqs <- seqs
    
  }
  
  
  # Collect random sequences if desired.
  if( all(args$collectRandomIDs != FALSE) ){
    random_seqs <- trimmed_seqs$randomSequences
    trimmed_seqs <- trimmed_seqs$trimmedSequences
  }
  
  # Log info
  lead_trimmed_tbl <- logSeqData(trimmed_seqs)
  cat("\nSequence information remaining after lead trimming:\n")
  print(lead_trimmed_tbl, row.names = FALSE)
  
  # Overread trimming
  if( nchar(args$overTrimSeq) > 0 ){
    
    # Determine percent identity from allowable mismatch.
    percent_id <- (nchar(args$overTrimSeq) - args$overMismatch) / 
      nchar(args$overTrimSeq)
  
    # Trim 3' end or overreading protion of sequences.
    trimmed_seqs <- trimOverreading(
      seqs = trimmed_seqs, 
      trim.sequence = args$overTrimSeq, 
      percent.id = percent_id, 
      max.seq.length = args$overMaxLength,
      min.seq.length = args$overMinLength
    )
    
    # Log info
    over_trimmed_tbl <- logSeqData(trimmed_seqs)
    cat("\nSequence information remaining after overreading trimming:\n")
    print(over_trimmed_tbl, row.names = FALSE)
    
  }
  
}else{
  
  # Split sequences up evenly across cores for trimming
  split.seqs <- split(
    seqs, ceiling(seq_along(seqs)/(length(seqs)/args$cores))
  )
  
  # Set up buster the cluster
  buster <- parallel::makeCluster(args$cores)
  
  # Trim 5' end or leading section of sequence while capturing random sequences,
  # if desired. Added features required workflow changes.
  if( nchar(args$leadTrimSeq) > 0 ){
    
    trimmed_seqs <- parallel::parLapply(
      buster,
      split.seqs,
      trimLeading,
      trim.sequence = args$leadTrimSeq,
      phasing = args$phasing,
      max.mismatch = args$leadMismatch,
      collect.random = all(args$collectRandomIDs != FALSE),
      filter = !args$noFiltering
    )
  
    if( all(args$collectRandomIDs != FALSE) ){
      
      random_seqs <- lapply(trimmed_seqs, "[[", "randomSequences")
      
      random_seqs <- lapply(seq_along(random_seqs[[1]]), function(i){
        serialAppendS4(
          lapply(seq_along(random_seqs), function(j){
            random_seqs[[j]][[i]]
          })
        )
      })

      trimmed_seqs <- lapply(trimmed_seqs, "[[", "trimmedSequences")
      
    }
    
  }else{
    
    trimmed_seqs <- split.seqs
    
  }

  trimmed_seqs <- serialAppendS4(trimmed_seqs)
  
  # Log info
  lead_trimmed_tbl <- logSeqData(trimmed_seqs)
  cat("\nSequence information remaining after lead trimming:\n")
  print(lead_trimmed_tbl, row.names = FALSE)
  
  # The method for overread trimming sequentially aligns shorter fragments of 
  # the overTrimSeq, and solely requiring mismatches could lead to some issues.
  # Therefore the same percent identity is requried across all alignments, 
  # however long.
  if( nchar(args$overTrimSeq) > 0 ){
    
    trimmed_seqs <- split(
      x = trimmed_seqs, 
      f = ceiling(seq_along(trimmed_seqs)/(length(trimmed_seqs)/args$cores))
    )
    
    percent_id <- (nchar(args$overTrimSeq) - args$overMismatch) / 
      nchar(args$overTrimSeq)
  
    # Trim 3' end or overreading protion of the sequence.
    trimmed_seqs <- parallel::parLapply(
      buster,
      trimmed_seqs,
      trimOverreading,
      trim.sequence = args$overTrimSeq, 
      percent.id = percent_id, 
      max.seq.length = args$overMaxLength,
      min.seq.length = args$overMinLength
    )

    trimmed_seqs <- serialAppendS4(trimmed_seqs)
    
    # Log info
    over_trimmed_tbl <- logSeqData(trimmed_seqs)
    cat("\nSequence information remaining after overreading trimming:\n")
    print(over_trimmed_tbl, row.names = FALSE)
    
  }
  
  # Stop buster before he gets out of control.
  parallel::stopCluster(buster)
  
}


# If no reads remaining, terminate and write output
if( length(seqs) == 0 ){
  
  cat(
    "\n  No reads remaining to trim. Terminating script after writing output.\n"
  )
  
  writeNullFile(
    file = args$output, 
    writeRandom = args$collectRandomIDs, 
    stat = args$stat,
    compress = args$compress
  )
  
  q()
  
}


# Second check for sequences below minimum length
trimmed_seqs <- trimmed_seqs[
  Biostrings::width(trimmed_seqs) >= args$minSeqLength
]

# Recover filtered reads if requested
if( args$noFiltering ){
  
  if( seq_type == "fasta" ){
    inputSeqs <- ShortRead::readFasta(args$seqFile)
  }else{
    inputSeqs <- ShortRead::readFastq(args$seqFile)
  }
  
  matched_idx <- which(id(inputSeqs) %in% id(trimmed_seqs))
  unmatched_idx <- which(!id(inputSeqs) %in% id(trimmed_seqs))
  untrimmed_seqs <- inputSeqs[unmatched_idx]
  output_seqs <- Biostrings::append(trimmed_seqs, untrimmed_seqs)
  output_seqs <- output_seqs[order(c(matched_idx, unmatched_idx))]
  
}else{
  
  output_seqs <- trimmed_seqs
  
}


# Log info
final_trimmed_tbl <- logSeqData(output_seqs)
cat("\nSequence information remaining:\n")
print(final_trimmed_tbl, row.names = FALSE)


# Write stats if requested
if( args$stat != FALSE ){
  
  sample_name <- unlist(strsplit(args$output, "/"))
  
  sample_name <- unlist(
    strsplit(sample_name[length(sample_name)], ".fa", fixed = TRUE)
  )[1]
  
  write.table(
    data.frame(
      sampleName = sample_name,
      metric = "reads",
      count = length(output_seqs)
    ),
    file = args$stat,
    sep = ",", 
    row.names = FALSE, 
    col.names = FALSE, 
    quote = FALSE
  )
  
}

# Collect RandomIDs if requested
if( all(args$collectRandomIDs != FALSE) ){
  
  random_seqs <- lapply(
    seq_along(random_seqs), 
    function(i, ids){
      
      random_seqs[[i]][
        which(as.character(ShortRead::id(random_seqs[[i]])) %in% ids)
      ]
      
    }, 
    ids = as.character(ShortRead::id(trimmed_seqs))
  )
  
}

# Sequences have been trimmed and random sequnces collected (if desired). 
# Next step is to write to output file(s).
# For fasta format, this is as simple as writing out the sequences currently in
# the environment. For fastq format, the quality scores for the trimmed bases
# must be loaded and trimmed as well.

# Write sequence file.
writeSeqFiles(
  seqs = output_seqs, 
  file = args$output,
  compress = args$compress
)

# Write randomID file.
if( all(args$collectRandomIDs != FALSE) ){
    
  if( length(args$collectRandomIDs) != length(random_seqs) ){
    
    new_file_name <- unlist(strsplit(args$collectRandomIDs[[1]], ".fa"))
    
    new_names <- paste0(
      new_file_name[[1]], ".", 1:length(random_seqs), ".", random_type
    )
    
    args$collectRandomIDs <- new_names
  }

  null <- mapply(
    writeSeqFiles,
    seqs = random_seqs,
    file = args$collectRandomIDs,
    MoreArgs = list(compress = args$compress)
  )
  
}

cat("\nScript completed.\n")

q()
