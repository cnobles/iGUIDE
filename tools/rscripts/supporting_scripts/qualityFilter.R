#' Quality filter read alignments
#' @param alignments GRanges object of alignment ranges with alignment metrics
#' of "matches", "repMatches", "qSize", "qStart", and "tBaseInsert".
#' @param q.start.max integer Maximum allowable nucleotide position to start read
#' alignment.
#' @param global.identity.min numeric Between 0 and 100 denoting the percent of 
#' global identity for the alignment needed to pass the filter.
#' @param base.insert.max integer Maximum number of allowable inserted nucleotides
#' within the alignment.
#' @return Subset of input GRanges Object passing filter criteria
#' @author Christopher Nobles, Ph.D. 
qualityFilter <- function(alignments, q.start.max = NULL, 
                          global.identity.min = NULL, base.insert.max = 5){
  
  stopifnot( class(alignments) == "data.frame" )
  
  statsNeeded <- c("matches", "repMatches", "qSize", "qStart", "tBaseInsert")
  
  stopifnot( all(statsNeeded %in% names(alignments)) )
  
  alignments$percIdent <- 
    100 * ( alignments$matches + alignments$repMatches ) / alignments$qSize
  
  if( length(q.start.max) > 0 ){
    alignments <- subset(
      alignments, alignments$percIdent >= global.identity.min
    )
  }
  
  if( length(global.identity.min) > 0 ){
    alignments <- subset(alignments, alignments$qStart <= q.start.max)
  }
  
  if( length(base.insert.max) > 0 ){
    alignments <- subset(alignments, alignments$tBaseInsert <= base.insert.max)
  }
  
  alignments
  
}
