#' Deterimine all the sites present in a GRanges Object where each row is one
#' read. Read counts ('counts') for each site will be returned in the metadata
#' column of the resulting GRanges object.
#' @param sites GRanges object
#' @param keep.cols character vector with names of metadata columns to retain 
#' after condensing the ranges. Unique values will be split upon and will 
#' delimit data found with different metadata.
#' @param list.bps logical Whether breakpoint lengths from the start site should
#' be retained in an IntegerList object in the metadata after condensing the 
#' ranges. Default is FALSE.
#' @param list.bp.counts logical Whether the read counts for each breakpoint
#' length should be included as an IntegerList along with list.bps. Default is
#' FALSE, but TRUE will force list.bps to be TRUE as well.
#' @return GRanges object with one or more metadata columns.
#' @author Christopher Nobles, Ph.D.
#' 

condenseSites <- function(sites, keep.cols = NULL, list.bps = FALSE, 
                          list.bp.counts = FALSE){
  
  stopifnot(class(sites) == "GRanges")
  
  if( list.bp.counts ) list.bps <- TRUE
  
  if( !is.null(keep.cols) ){
  
    keep_cols <- as.data.frame(
      GenomicRanges::mcols(sites)[
        match(keep.cols, names(GenomicRanges::mcols(sites)))
      ]
    )

    split_vec <- split( seq_along(sites), keep_cols )
    
    sites_list <- lapply( split_vec, function(x){
      GenomicRanges::granges(sites[x]) 
    })
    
    col_list <- lapply( split_vec, function(x){
      GenomicRanges::as.data.frame(keep_cols[x,])
    } )
    
    col_list <- lapply( col_list, function(x){
      names(x) <- keep.cols
      return(x)
    } )

  }else{
    
    sites_list <- GenomicRanges::GRangesList(sites)
    col_list <- NULL
    
  }
  
  unlist(GenomicRanges::GRangesList(lapply(
    seq_along(sites_list), 
    function(i, bps, bp.counts){
      
      gr <- sites_list[[i]]
      cols <- col_list[[i]]
      
      gr_red <- GenomicRanges::reduce(
        flank(gr, -1, start = TRUE), 
        min.gapwidth = 0L, 
        with.revmap = TRUE
      )
    
      revmap <- gr_red$revmap
      gr_red$revmap <- NULL
      
      bp_df <- data.frame(
        redid = S4Vectors::Rle(
          values = seq_along(revmap), 
          lengths = width(revmap@partitioning)
        ),
        oriid = unlist(revmap),
        width = width(gr[unlist(revmap)])
      )
      
      bp_df <- bp_df[order(bp_df$width),]
      bp_rle <- split(S4Vectors::Rle(bp_df$width), bp_df$redid)
      
      if( !is.null(cols) ) GenomicRanges::mcols(gr_red) <- unique(cols)      
      
      gr_red$counts <- S4Vectors::width(revmap@partitioning)
      gr_red$fragLengths <- S4Vectors::width(
        S4Vectors::runValue(bp_rle)@partitioning
      )
      
      if( bps ){
        gr_red$bpWidths <- S4Vectors::runValue(bp_rle)
        if(bp.counts) gr_red$bpCounts <- S4Vectors::runLength(bp_rle)
      }
      
      return(gr_red)
      
    },
    bps = list.bps, 
    bp.counts = list.bp.counts
  )))

}