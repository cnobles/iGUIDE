#' Trim 3' ends of sequences when the sequence has overread into
#' synthetic sequence
#' \code{trimOverreading} removes 3' ends of nucleotide sequences completely or
#' partially matching the trim.sequence.
#' @param seqs ShortReadQ object of reads or unique sequences
#' @param trim.sequence character string of lenth 1, such as "GAAAATC". This
#' string will be used to match the end of sequences, upon which the matching
#' portion will be trimmed from the start of the match to the end of the
#' sequence. Ambiguous nucleotides within the sequence can be used for 
#' alignment, but matching sequences are not recored.
#' @param percent.id numeric between 0 and 1.0 denoting the minimum percent
#' identity acceptable for an matching alignment.
#' @param max.seq.length integer the maximum length to consider of the
#' trim.sequence to use for alignments. Using the full length of sequence
#' avaliable for matching can many times be computationally intensive and
#' unnessesarily time consuming. Further, identical results can be obtained
#' using only a portion of the sequence. For example, setting the max.seq.length
#' to 15L will only use the first 15 nucleotides of the trim.sequence.
#' @param min.seq.length integer the minimum length to consider of the 
#' trim.sequence to use for alignments. Default 3.
#' @author Christopher Nobles, Ph.D.

trimOverreading <- function(seqs, trim.sequence, percent.id, 
                             max.seq.length = NULL, min.seq.length = 3){
  
  stopifnot(class(seqs) %in% c("ShortReadQ", "ShortRead"))
  stopifnot(!is.null(ShortRead::id(seqs)))

  # Trim down trim.sequence if max.seq.length provided
  if( !is.null(max.seq.length) ){
    
    trim_sequence <- as.character(Biostrings::DNAStringSet(
      x = trim.sequence, 
      start = 1L, 
      end = min(nchar(trim.sequence), max.seq.length)
    ))
    
  }
  
  trim_seqs <- sapply(
    0:( nchar(trim_sequence) - min.seq.length ), 
    function(i){
      substr(trim_sequence, 1, nchar(trim_sequence) - i)
    }
  )

  alignments <- do.call(
    c, 
    lapply(trim_seqs, function(trim_seq, seqs, percent.id){
    
        mismatch <- round( nchar(trim_seq) - percent.id*nchar(trim_seq) )
        
        vmp <- Biostrings::vmatchPattern(
          pattern = trim_seq, 
          subject = ShortRead::sread(seqs), 
          max.mismatch = mismatch, 
          fixed = FALSE
        )
        
        idx <- which( BiocGenerics::lengths(vmp) >= 1 )
        len <- BiocGenerics::lengths(vmp[idx])
        len <- len[ which(len >= 1) ]
        idx <- S4Vectors::Rle(values = idx, lengths = len)
        ir <- unlist(vmp)
        names(ir) <- idx
        ir
        
      }, 
      seqs = seqs, 
      percent.id = percent.id
    )
  )
  
  seq_widths <- Biostrings::width(seqs[as.numeric(names(alignments))])
  
  alignments <- alignments[
    IRanges::width(alignments) == nchar(trim_sequence) | 
    IRanges::end(alignments) == seq_widths
  ]
  
  alignments <- IRanges::IRangesList(
    split(x = alignments, f = names(alignments))
  )
  
  alignments <- unlist(IRanges::reduce(alignments, min.gapwidth = 1000000))
  alignments <- alignments[order(as.numeric(names(alignments)))]
  idx <- as.numeric(names(alignments))
  
  if( any(table(idx) > 1) ){
    stop("Issues with overread trimming, please adjust input parameters.")
  }
  
  # Trim sequences
  seqs[idx] <- IRanges::narrow(
    seqs[idx], 
    end = IRanges::start(alignments) - 1
  )
  
  return(seqs)
}
