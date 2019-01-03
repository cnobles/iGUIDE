#' Trim beginning or leading ends of nucleotide sequences
#'
#' @param seqs ShortReadQ object of reads or unique sequences
#' @param trim.sequence character string of lenth 1, such as "GAAAATC". This
#' string will be used to match to the beginning of sequences, upon which
#' non-matching sequences will be discarded and the matching portion will be
#' trimmed from the leading side of the sequence. Ambiguous nucleotides within
#' the sequence will be used to determine random sequences and can be used for
#' alignment or collecting random nucleotide sequences embedded within the
#' trim.sequence structure.
#' @param phasing integer/numeric value, denoting the number of nucleotides used
#' for phasing while sequencing. This number of nucleotides will be removed from
#' the beginning of the sequence before any alignment.
#' @param max.mismatch integer/numeric value or vector. Values indicate the
#' number of mismatches allowed within the DNA segment. Integer / numeric 
#' vectors can be used to indicate the number of allowable mismatches within 
#' segments of non-ambiguous nucleotices. The length of input vecors must be 
#' equal to the number of non-ambiguous segments within the trim.sequence. If 
#' no ambiguous segments are present in trim.sequence vectors for max.mismatch 
#' will be summed to determine the maximum allowable mismatch over the entire 
#' string.
#' @param collect.random logical should random / ambiguous protions of the
#' trim.sequence be collected? They will be returned as a listed DNAStringSet
#' under 'randomSequences' in order from left to right of appearance within
#' trim.sequence.
#' @param filter logical (default: TRUE). If TRUE, sequences not matching the 
#' trim.sequence will be dropped or filtered from the output. If FALSE, all
#' input sequences are returned, but matching sequences are trimmed.
#' 
#' @return DNAStringSet of sequences with trim.sequence removed or a listed
#' object with the first postion being the DNAStringSet of trimmed sequences and
#' the second being the random sequences collected during trimming.
#' 
#' @author Christopher Nobles, Ph.D.
#'  

trimLeading <- function(seqs, trim.sequence, phasing = 0L, max.mismatch = 1L,
                         collect.random = FALSE, filter = TRUE){
  
  # Checks and requirements
  stopifnot(!is.null(ShortRead::id(seqs)))
  
  # Change scientific notation switch to inhibit indexing errors
  ori_scipen <- getOption("scipen")
  options(scipen = 99)
  
  # Phasing will ignore the first number of nucleotides of the sequence
  seqs <- IRanges::narrow(seqs, start = 1L + phasing)
  
  # Determine the structure of ambiguous sequences within the trim.sequence
  ambi_present <- stringr::str_detect(trim.sequence, pattern = "[^A^T^G^C]")
  
  if( ambi_present & collect.random ){
    
    trim_segments <- unlist(strsplit(trim.sequence, "[^A^T^G^C]+"))
    
    trim_seg_ir <- sapply(
      trim_segments, Biostrings::vmatchPattern, subject = trim.sequence
    )
    
    trim_seg_ir <- IRanges::IRanges(
      start = sapply(trim_seg_ir, function(x) as.integer(IRanges::start(x))),
      end = sapply(trim_seg_ir, function(x) as.integer(IRanges::end(x))),
      names = names(trim_seg_ir)
    )
    
    rand_ir <- IRanges::setdiff(
      IRanges::IRanges(start = 1, end = nchar(trim.sequence)),
      trim_seg_ir)
    
  }else{
    
    trim_segments <- trim.sequence
    trim_seg_ir <- unlist(
      Biostrings::vmatchPattern(trim_segments, trim.sequence))
    names(trim_seg_ir) <- trim_segments
    
  }
  
  # Set allowable mismatches for each range
  if( length(max.mismatch) == 1 ){
  
    trim_seg_ir@metadata$misMatch <- rep(max.mismatch, length(trim_seg_ir))
    
  }else if(
    length(max.mismatch) == length(trim_seg_ir) & is.numeric(max.mismatch)
  ){
    
      trim_seg_ir@metadata$misMatch <- max.mismatch
      
  }else{
    
    stop(
      "\nThe variable max.mismatch needs to be either a single integer or 
      integer vector of length equal to fixed fragments within the 
      trim.sequence."
    )
    
  }
  
  # Remove seqs that do not have enough sequence for analysis
  # Cutoff = length(trim.sequence)
  seqs <- seqs[Biostrings::width(seqs) > IRanges::width(range(trim_seg_ir))]
  
  lead_seqs <- IRanges::narrow(
    seqs, 
    start = 1, 
    end = ifelse(
      Biostrings::width(seqs) >= nchar(trim.sequence) + 1, 
      rep(nchar(trim.sequence) + 1, length(seqs)), 
      Biostrings::width(seqs)
    )
  )
  
  # Align whole sequence to 5' end of sequence
  aln <- Biostrings::vmatchPattern(
    pattern = trim.sequence, 
    subject = ShortRead::sread(lead_seqs), 
    max.mismatch = sum(max.mismatch), 
    fixed = FALSE
  )
  
  matched_idx <- which(S4Vectors::lengths(aln) == 1)
  
  # Serially align segment(s) from trim.sequence to seqs
  if( length(trim_seg_ir) > 1 ){
  
    aln_segs <- lapply(
      seq_along(trim_seg_ir), 
      function(i, trim_seg_ir, seqs){
        
        t_seq <- names(trim_seg_ir[i])
        mismatch_i <- trim_seg_ir@metadata$misMatch[i]
        
        aln_seqs_i <- IRanges::narrow(
          seqs,
          start = ifelse(
            IRanges::start(trim_seg_ir[i]) == 1L, 
            1, 
            IRanges::start(trim_seg_ir[i]) - 1
          ),
          end = IRanges::end(trim_seg_ir[i]) + 1
        )
        
        aln <- Biostrings::vmatchPattern(
          pattern = t_seq, 
          subject = ShortRead::sread(aln_seqs_i), 
          max.mismatch = mismatch_i, 
          fixed = FALSE
        )
      
        if( any(S4Vectors::lengths(aln) > 1) ){
          stop(
            "\nAlignment too permissive. Ambiguous mapping of sequences.
             Please adjust max.mismatch criteria."
          )
        }
        
        idx <- S4Vectors::lengths(aln) == 1 
        
        return(list("match" = aln, "idx" = idx))
        
      },
      trim_seg_ir = trim_seg_ir,
      seqs = lead_seqs
    )
    
    # Identify and trim only sequences that have all required alignments
    seg_idx <- table(unlist(
      lapply(lapply(aln_segs, "[[", "idx"), which)
    ))
    
    seg_idx <- as.numeric(names(seg_idx)[seg_idx == length(trim_seg_ir)])
    
    matched_idx <- matched_idx[matched_idx %in% seg_idx]
    
  }
  
  # Isolate sequences matching the input criteria
  matched_seqs <- seqs[matched_idx]
  trim_shift <- ifelse(length(trim_seg_ir) > 1, 2, 1)
  matched_starts <- unlist(aln@ends[matched_idx]) + 1
  trimmed_seqs <- IRanges::narrow(matched_seqs, start = matched_starts)

  if( !filter ){
    
    unmatched_idx <- which(!seq_along(seqs) %in% matched_idx)
    untrimmed_seqs <- seqs[unmatched_idx]
    trimmed_seqs <- Biostrings::append(trimmed_seqs, untrimmed_seqs)
    trimmed_seqs <- trimmed_seqs[order(c(matched_idx, unmatched_idx))]
    
  }
  
  if( !collect.random | !ambi_present ){
    
    # Return scipen option to original value
    options(scipen = ori_scipen)
    return(trimmed_seqs)
    
  }else{
    
    random_seqs <- lapply(
      seq_along(rand_ir), 
      function(k, lead_seqs){
        
        region <- rand_ir[k]
        
        reg_seq <- Biostrings::subseq(
          trim.sequence, 
          start = IRanges::start(region), 
          end = IRanges::end(region)
        )
        
        rand_seqs <- IRanges::narrow(
          lead_seqs[matched_idx], 
          start = IRanges::start(region), 
          end = IRanges::end(region)
        )
        
        matched_random_idx <- Biostrings::vmatchPattern(
          pattern = reg_seq, 
          subject = ShortRead::sread(rand_seqs), 
          fixed = FALSE
        )
        
        matched_random_idx <- which(S4Vectors::lengths(matched_random_idx) == 1)
        
        rand_seqs[matched_random_idx]
        
      }, 
      lead_seqs = lead_seqs
    )

    # Return scipen option to original value
    options(scipen = ori_scipen)
    
    return(
      list(
        "trimmedSequences" = trimmed_seqs,
        "randomSequences" = random_seqs
      )
    )
    
  }
}
