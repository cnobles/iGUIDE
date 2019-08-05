#' Convert matches to GRanges objects
#'
#' @usage mindexToGranges(mindex, strand, ref = NULL)
#'
#' @param mindex an MIndex object, such as produced by
#' Biostrings::vmatchPattern.
#'
#' @param strand character to specify the strand of the alignment, either "+" or
#' "-" or "*".
#'
#' @param ref reference genome object, sequence information will be included
#' with output GRanges object if this parameter is specified.
#'
#' @description Converts MIndex objects into GRange objects. Helpful conversion
#' if using the Biostrings::vmatchPattern function and genome objects.
#'
#' @author Christopher Nobles, Ph.D.

mindexToGranges <- function(mindex, strand, ref = NULL){

  if( is.null(mindex@NAMES) ) stop("NAMES column not found for seqnames.")
  ir <- unlist(mindex)

  if( is.null(ref) ){

    return(unname(GenomicRanges::GRanges(
      seqnames = names(ir),
      ranges = ir,
      strand = rep(strand, length(ir))
    )))

  }else{

    return(unname(GenomicRanges::GRanges(
      seqnames = names(ir),
      ranges = ir,
      strand = rep(strand, length(ir)),
      seqinfo = GenomeInfoDb::seqinfo(ref)
    )))

  }

}
