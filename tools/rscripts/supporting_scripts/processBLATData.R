#' Converts alignment data into a GRanges object without information loss. BLAT 
#' aligns on a 0-base system, while all genomic data is read on a 1-base system.
#' 
#' @param algns data.frame of psl table, containing labeled columns. 
#' @param from character, which read is the algns object from? ("anchor" or 
#' "adrift")
#' @param ref.genome BSgenome object with seqInfo of reference genome.
#' @author Christopher Nobles, Ph.D.
#' 

processBLATData <- function(algns, from, ref.genome){
  
  stopifnot(from == "anchor" | from == "adrift")
  
  algns$from <- from
  
  algns$qtStart <- ifelse(
    algns$strand == "+",
    ( algns$tStart - (algns$qStart) ),
    ( algns$tStart - (algns$qSize - algns$qEnd - 1) )
  )
  
  algns$qtEnd <- ifelse(
    algns$strand == "+",
    ( algns$tEnd + (algns$qSize - algns$qEnd - 1) ),
    ( algns$tEnd + (algns$qStart) )
  )    
  
  algns_gr <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(algns$tName),
    ranges = IRanges::IRanges(
      start = ( algns$qtStart + 1 ), 
      end = ( algns$qtEnd ) #Convert to 1-base
    ), 
    strand = S4Vectors::Rle(algns$strand),
    seqinfo = GenomeInfoDb::seqinfo(ref.genome))
  
  GenomicRanges::mcols(algns_gr) <- algns[,c(
    "from", "qName", "matches", "repMatches", "misMatches", 
    "qStart", "qEnd", "qSize", "tBaseInsert"
  )]
  
  algns_gr
  
}
