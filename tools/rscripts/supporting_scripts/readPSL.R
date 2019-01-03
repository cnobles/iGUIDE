#' read psl.gz files, assuming psl gz files don't have column header
#' @param psl.file character vector of file name(s).
#' @param to.null  character vector of column names to remove.
#' @return data.frame, data.table of the psl table
#' @author Christopher Nobles, Ph.D. and Yinghua Wu, Ph.D.

readPSL <- function(psl.file, to.null = NULL) {
  
  cols <- c(
    "matches", "misMatches", "repMatches", "nCount", "qNumInsert", 
    "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName", "qSize", 
    "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd", "blockCount", 
    "blockSizes", "qStarts", "tStarts")
  
  cols_class <- c(
    rep("numeric",8), rep("character",2), rep("numeric",3), "character", 
    rep("numeric",4), rep("character",3))
  
  psl <- suppressMessages(
    lapply(psl.file, function(f){
      message("Reading ",f)
      try(data.table::fread(paste("zcat", f), sep = "\t"), silent = TRUE)
    })
  )
  
  psl <- lapply(psl, function(p){
    
    if( any(class(p) == "try-error") ){

      if(grepl("File is empty:", p[1])){
        p <- read.table(text = "", col.names = cols, colClasses = cols_class)
        return(data.table::as.data.table(p))
      }else{
        stop("Error in loading psl files. Check output from alignment (blat).")
      }
      
    }else{
      return(p)
    }
  
  })
  
  psl <- data.table::rbindlist(psl)
  colnames(psl) <- cols
  
  if(length(to.null)>0) psl[, to.null] <- NULL
  
  return(as.data.frame(psl))
}
