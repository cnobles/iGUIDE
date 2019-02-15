#' For anyone reviewing the code below, the following is a small style guide 
#' outlining the various formats for the code. 
#' 
#' Names with "_": objects, inlucding data.frames, GRanges, vectors, ...
#' Names in caMel format: functions or components of objects (i.e. columns 
#' within a data.frame).
#' Names with ".": arguments / options for functions

# Set Global options and load intiial packages ---------------------------------
options(stringsAsFactors = FALSE, scipen = 99, width = 999)

code_dir <- dirname(sub(
  pattern = "--file=", 
  replacement = "", 
  x = grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
))

desc <- yaml::yaml.load_file(
  file.path(code_dir, "descriptions/filt.desc.yml")
)


# Set up arguments and workflow of script --------------------------------------
## Argument parser =============================================================
parser <- argparse::ArgumentParser(
  description = desc$program_short_description,
  usage = "nuc filt <seqFile(s)> [-h/--help, -v/--version] [optional args]"
)

parser$add_argument(
  "seqFile", nargs = "+", type = "character", help = desc$seqFile
)

parser$add_argument(
  "-o", "--output", nargs = "+", type = "character", help = desc$output
)

parser$add_argument(
  "-i", "--index", nargs = "*", type = "character", help = desc$index
)

parser$add_argument(
  "-s", "--seq", nargs = "*", type = "character", help = desc$seq
)

parser$add_argument(
  "-m", "--mismatch", nargs = "+", type = "integer", default = 0, 
  help = desc$mismatch
)

parser$add_argument(
  "-r", "--refseqs", nargs = "+", type = "character", help = desc$refseqs
)

parser$add_argument(
  "--aligntype", nargs = 1, type = "character", default = "ov",
  help = desc$aligntype
)

parser$add_argument(
  "--pctID", nargs = 1, type = "integer", default = 95,
  help = desc$pctID
)

parser$add_argument(
  "--pctIDtype", nargs = 1, type = "character", default = "global", 
  help = desc$pctIDtype
)

parser$add_argument(
  "--subMatAdj", nargs = "+", type = "character", default = FALSE,
  help = desc$subMatAdj
)

parser$add_argument(
  "--gapOpen", nargs = 1, type = "integer", default = 10,
  help = desc$gapOpen
)

parser$add_argument(
  "--gapExt", nargs = 1, type = "integer", default = 4,
  help = desc$gapExt
)

parser$add_argument(
  "--minAlignLength", nargs = 1, type = "integer", default = 20,
  help = desc$minAlignLength
)

parser$add_argument(
  "--readNamePattern", nargs = 1, type = "character", default = "[\\w:-]+",
  help = desc$readNamePattern
)

parser$add_argument(
  "-c", "--cores", nargs = 1, default = 1, type = "integer", help = desc$cores
)

parser$add_argument(
  "--stat", nargs = 1, default = FALSE, type = "character", help = desc$stat
)

parser$add_argument(
  "--header", action = "store_true", help = desc$header
)

parser$add_argument(
  "-n", "--negSelect", action = "store_true", help = desc$negSelect
)

parser$add_argument(
  "--any", action = "store_true", help = desc$any
)

parser$add_argument(
  "--compress", action = "store_true", help = desc$compress
)

parser$add_argument(
  "-q", "--quiet", action = "store_true", help = desc$quiet
)



## Parse cmd line args =========================================================
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

## Checks and balance ==========================================================
if( args$cores > 1 ){
  
  # Stop code since parallel operation has not been constructed yet
  stop("\n  Parallel options have not yet been implemented.\n")
  
  if( args$cores > parallel::detectCores() ){
  
    cat(
      "\n  Requested cores is greater than availible for system.",
      "Changing cores to max allowed.\n"
    )
    args$cores <- detectCores()
    
  }
  
}else if( args$cores < 1 ){
  
  args$cores <- 1
  
}

if( length(args$seqFile) != length(args$output) ){
  stop(
    "\n  The same number of input and output file names need to be provided.\n")
}

if( length(args$index) > 1 ){
  stop(
    "\n  Only one index file can be used at a time. ",
    "Please consolidate indices.\n"
  )
}

if( length(args$mismatch) != length(args$seq) ){
  args$mismatch <- rep(args$mismatch[1], length(args$seq))
}

if( length(args$seq) > 0 ){
  
  args$seq <- toupper(gsub("U", "T", args$seq))
  
  if( 
    any(!unlist(strsplit(paste(args$seq, collapse = ""), "")) %in% 
      names(Biostrings::IUPAC_CODE_MAP)) 
  ){
    stop("\n  Unknown nucleotides detected in input filtering sequence(s).\n")
  }
  
}

if( !args$pctIDtype %in% c("global", "local") ){
  stop("\n  Input '--pctIDtype' must be either 'local' or 'global' [default].")
}

# Determine input sequence file type(s)
seq_type <- unlist(strsplit(args$seqFile, "/"))
seq_type <- seq_type[length(seq_type)]
seq_type <- stringr::str_extract(seq_type, ".fa[\\w]*")

if( any(!seq_type %in% c(".fa", ".fq", ".fasta", ".fastq")) ){
  
  stop(
    "\n  Unrecognized sequence file type, please convert to '*.fasta' or ", 
    "'*.fastq'. Gzip compression is acceptable as well.\n"
  )
  
}

seq_type <- ifelse(seq_type %in% c(".fa", ".fasta"), "fasta", "fastq")

# Determine sequence output file type(s)
if( length(args$output) > 0 ){
  
  out_type <- unlist(strsplit(args$output, "/"))
  out_type <- out_type[length(out_type)]
  out_type <- stringr::str_extract(out_type, ".fa[\\w]*")
  
  if( any(!out_type %in% c(".fa", ".fq", ".fasta", ".fastq")) ){
    
    stop(
      "\n  Unrecognized output sequence file type, please change to ", 
      "'*.fasta' or '*.fastq'.\n"
    )
    
  }
  
  out_type <- ifelse(out_type %in% c(".fa", ".fasta"), "fasta", "fastq")
  
}

# Identify filtering type
select_methods <- c()
if( length(args$index) == 1 ) select_methods <- c(select_methods, 1)
if( length(args$seqFile) > 1 ) select_methods <- c(select_methods, 2)
if( length(args$seq) > 0 ) select_methods <- c(select_methods, 3)
if( length(args$refseqs) > 0 ) select_methods <- c(select_methods, 4)

methods <- c(
  "input indices", "multiple file input indices", 
  "sequence content", "sequence matching reference(s)"
)[select_methods]

filt_type <- paste0(
  ifelse(args$negSelect, "negative", "positive"), 
  " selection using ", 
  paste(methods, collapse = ifelse(args$any, " or ", " and ")), "."
)


## Input arguments table =======================================================
input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(seq_along(args), function(i){
    paste(args[[i]], collapse = ", ")
  })
)

input_table <- input_table[
  match(
    c("seqFile :", "output :", "index :", "header :", "negSelect :", "seq :", 
      "mismatch :", "refseqs :", "aligntype :", "pctID :", "pctIDtype :", 
      "subMatAdj :", "gapOpen :", "gapExt :", "minAlignLength :", 
      "readNamePattern :", "compress :", "cores :"), 
    input_table$Variables)
  ,]

if( !args$quiet ){
  
  cat("\nFilter Inputs:\n")
  print(
    data.frame(input_table, row.names = NULL), 
    right = FALSE, 
    row.names = FALSE
  )
  cat("\n  Filtering methods include", filt_type, "\n")
  
}


# Additional supporting functions ----------------------------------------------
source(file.path(code_dir, "supporting_scripts", "writeSeqFiles.R"))
source(file.path(code_dir, "supporting_scripts", "nucleotideScoringMatrices.R"))
source(file.path(code_dir, "supporting_scripts", "substituteAdjustments.R"))
source(file.path(code_dir, "supporting_scripts", "utility_funcs.R"))

#' Filter sequences based on input arguments
#' This function is the basis for the script.
filterSeqFile <- function(input.seqs, args){

  ## Identify sequence names matching across multiple sequence files
  if( length(input.seqs) > 1 ){
    
    multi_input_ids <- lapply(input_seqs, function(seq){
      stringr::str_extract(
        string = as.character(unique(ShortRead::id(seq))), 
        pattern = args$readNamePattern
      )
    })
    
    multi_input_tbl <- table(unlist(multi_input_ids))
    
    if( args$negSelect ){
      multi_input_names <- names(multi_input_tbl)[which(multi_input_tbl == 1)]
    }else if( args$any ){
      multi_input_names <- names(multi_input_tbl)[which(multi_input_tbl > 1)]
    }else{
      multi_input_names <- names(multi_input_tbl)[
        which(multi_input_tbl == length(input_seqs))
      ]
    }
    
    multi_filter_idx <- lapply(input_seqs, function(seqs, idx){
        
        ids <- stringr::str_extract(
          string = as.character(ShortRead::id(seqs)), 
          pattern = args$readNamePattern
        )
        
        which(ids %in% idx)
        
      }, 
      idx = multi_input_names
    )
    
  }
  
  
  ## Identify sequence names by matching to index file
  if( length(args$index) == 1 ){
    
    input_ids <- lapply(input_seqs, function(seq){
      stringr::str_extract(
        string = as.character(ShortRead::id(seq)), 
        pattern = args$readNamePattern
      )
    })
    
    index_df <- read.delim(args$index, header = args$header)
    
    index <- stringr::str_extract(
      string = as.character(index_df[,1]), 
      pattern = args$readNamePattern
    )
    
    index_filter_idx <- lapply(input_ids, function(ids, idx){
        which(ids %in% idx)
      }, 
      idx = index
    )
    
  }
  
  
  ## Identify sequences by matching input nucleotide sequence
  if( length(args$seq) > 0 ){
    
    seq_filter_idx <- lapply(
      input_seqs, 
      function(seqs, pattern, mismatch, neg, any){
        
        vcp <- mapply(function(pat, mis, seqs, neg){
          
            v <- Biostrings::vcountPattern(
              pat, ShortRead::sread(seqs), max.mismatch = mis, fixed = FALSE)
          
            if( neg ){
              return(which(v == 0))
            }else{
              return(which(v > 0))
            }
          
          }, 
          pat = pattern, mis = mismatch, 
          MoreArgs = list(seqs = seqs, neg = neg),
          SIMPLIFY = FALSE
        )
        
        vcp_tbl <- table(unlist(vcp))
        
        if( any ){
          return(as.numeric(names(vcp_tbl[which(vcp_tbl >= 1)])))
        }else{
          return(as.numeric(names(vcp_tbl[which(vcp_tbl == length(pattern))])))
        }
        
      }, 
      pattern = args$seq, 
      mismatch = args$mismatch, 
      neg = args$negSelect, 
      any = args$any
    )
    
  }
  
  
  ## Identify sequence that match to reference sequence(s)
  if( length(args$refseqs) > 0 ){
    
    ref_filter_idx <- lapply(
      input.seqs,
      function(seqs, refs, alntype, pctID, idtype, subadj, gapOpen, gapExt, 
               minAlignLength, neg){
        
        # Load reference sequences
        ref_types <- unlist(strsplit(refs, "/"))
        ref_types <- ref_types[length(ref_types)]
        ref_types <- stringr::str_extract(ref_types, ".fa[\\w]*")
        
        if( any(!ref_types %in% c(".fa", ".fq", ".fasta", ".fastq")) ){
          
          stop(
            "\n  Unrecognized sequence file type, please convert to '*.fasta' or ", 
            "'*.fastq'. Gzip compression is acceptable as well.\n"
          )
          
        }
        
        ref_types <- ifelse(ref_types %in% c(".fa", ".fasta"), "fasta", "fastq")
        
        refs <- mapply(
          function(file, file_type){
            
            if( file_type == "fasta" ){
              return(ShortRead::readFasta(file))
            }else{
              return(ShortRead::readFastq(file))
            }
            
          }, 
          file = refs, 
          file_type = ref_types, 
          SIMPLIFY = FALSE
        )
        
        if( length(refs) > 1 ){
          refs <- serialAppendS4(refs)
        }else{
          refs <- refs[[1]]
        }

        # Alignment type
        align_types <- structure(
          c("global", "local", "overlap", "global-local", "local-global"),
          names = c("gg", "ll", "ov", "gl", "lg")
        )
        
        alntype <- align_types[alntype]
        
        # Score only?
        SO <- idtype == 'global'
        
        # Interpret adjustment if any
        input_adjs <- stringr::str_extract(
          subadj[grep("^i", subadj)], "[\\w]{2}$"
        )
        
        refer_adjs <- stringr::str_extract(
          subadj[grep("^r", subadj)], "[\\w]{2}$"
        )
        
        # Apply adjustments and convert to character vectors
        seqs <- substituteAdjustments(ShortRead::sread(seqs), input_adjs)
        refs <- substituteAdjustments(ShortRead::sread(refs), refer_adjs)
        
        alignments <- lapply(
          refs, 
          function(ref, seqs, alntype, gapOpen, gapExt, SO){
            
            Biostrings::pairwiseAlignment(
              pattern = seqs, 
              subject = ref, 
              type = alntype, 
              gapOpening = gapOpen,
              gapExtension = gapExt,
              substitutionMatrix = usanmat(),
              scoreOnly = SO
            )
            
          },
          seqs = Biostrings::DNAStringSet(seqs),
          alntype = alntype,
          gapOpen = gapOpen,
          gapExt = gapExt,
          SO = SO
        )
        
        if( idtype == "global" ){
          
          max_score <- 100 * apply(
            matrix(unlist(alignments), ncol = length(refs)), 1, max
          ) / nchar(seqs)
          
        }else if( idtype == "local" ){
          
          local_score <- matrix(
            unlist(lapply(alignments, function(x) Biostrings::score(x))), 
            ncol = length(refs)
          )
          
          local_size <- matrix(
            unlist(lapply(alignments, function(x) x@pattern@range@width)),
            ncol = length(refs)
          )
          
          top_score_idx <- apply(100 * local_score / local_size, 1, function(x){ 
            which(x == max(x)) 
          })
          
          top_score_len <- unlist(lapply(
            seq_along(top_score_idx), 
            function(i){ 
              unique(local_size[i, top_score_idx[[i]], drop = TRUE])
            }
          ))
          
          top_score <- unlist(lapply(
            seq_along(top_score_idx), 
            function(i){ 
              unique(local_score[i, top_score_idx[[i]], drop = TRUE])
            }
          ))
          
          max_score <- 100 * top_score / top_score_len
          
        }else{
          
          stop("\n  Input error, pctIDtype must be either 'local' or 'global'.")
          
        }
        
        if( neg ){
          return( which(max_score < pctID | top_score_idx < minAlignLength) )
        }else{
          return( which(max_score >= pctID & top_score_len >= minAlignLength) )
        }
        
      },
      refs = args$refseqs,
      alntype = args$aligntype,
      pctID = args$pctID,
      idtype = args$pctIDtype,
      subadj = args$subMatAdj,
      gapOpen = args$gapOpen,
      gapExt = args$gapExt,
      minAlignLength = args$minAlignLength,
      neg = args$negSelect
    )
    
  }
  
  # Consolidate indices from each method employed 
  lapply(seq_along(input_seqs), function(i){
    
    idx <- NULL
    cnt <- 0
    
    if( exists("multi_filter_idx") ){
      cnt <- cnt + 1
      idx <- c(idx, multi_filter_idx[[i]]) 
    }
    
    if( exists("index_filter_idx") ){ 
      cnt <- cnt + 1
      idx <- c(idx, index_filter_idx[[i]]) 
    }
    
    if( exists("seq_filter_idx") ){ 
      cnt <- cnt + 1
      idx <- c(idx, seq_filter_idx[[i]]) 
    }
    
    if( exists("ref_filter_idx") ){
      cnt <- cnt + 1
      idx <- c(idx, ref_filter_idx[[i]]) 
    }
    
    if( args$any ){
      return(unique(idx))
    }else{
      idx_tbl <- table(idx)
      return(as.numeric(names(idx_tbl)[idx_tbl == cnt]))
    }
    
  })
  
}


# Identify indices of input file(s) for filtering ------------------------------
input_seqs <- mapply(
  function(file, file_type){
    
    if( file_type == "fasta" ){
      return(ShortRead::readFasta(file))
    }else{
      return(ShortRead::readFastq(file))
    }
    
  }, 
  file = args$seqFile, 
  file_type = seq_type, 
  SIMPLIFY = FALSE
)

output_indices <- filterSeqFile(input_seqs, args)

output_seqs <- mapply(
  function(seqs, idx){ seqs[idx] }, 
  seqs = input_seqs, 
  idx = output_indices, 
  SIMPLIFY = FALSE
)


# Write output files -----------------------------------------------------------
if( args$stat != FALSE ){
  
  sample_name <- strsplit(args$output, "/", fixed = TRUE)
  sample_name <- mapply("[[", sample_name, lengths(sample_name))
  sample_name <- strsplit(sample_name, ".fa", fixed = TRUE)
  sample_name <- mapply("[[", sample_name, 1)
  
  write.table(
    data.frame(
      sampleName = sample_name,
      metric = "reads",
      count = lengths(output_seqs)
    ),
    file = args$stat,
    sep = ",", 
    row.names = FALSE, 
    col.names = FALSE, 
    quote = FALSE
  )
  
}
 
 
null <- sapply(args$output, unlink)

null <- mapply(
  writeSeqFiles, 
  seqs = output_seqs, 
  file = args$output, 
  MoreArgs = list(compress = args$compress)
)

q()

# Notes:
## Work flow: R1xR2 --> index --> sequence
## 
## For output from each step, pass on a list of indices for the ShortRead object
## This will allow for filtering and then just returning an object that contains
## the correct indices. 
## 
## For R1xR2 type filtering: intersect / c(id()) | duplicated / table / which(bool)
## For index type filtering: which(id() %in% index) or which(!id() %in% index)
## For seq type filtering  : v[match/count]pattern(fixed = FALSE) | which(lengths() >= 1 or lengths() == 0)
