#' Read key files from multiple formats.
#' @param key.file character File name (path/to/name.format).
#' @param format character Acceptable formats include: "csv", "tsv", "rds", and 
#' "RData".
#' @author Christopher Nobles, Ph.D.
#'

readKeyFile <- function(key.file, format){
  
  cols <- c("readNames", "seqID")
  cols.class <- c("character", "character")
  
  if( format == "csv" ){
    
    file <- try(
      data.table::fread(key.file, sep = ",", data.table = FALSE), silent = TRUE
    )
    
    if( any(class(file) == "try-error") ){

      if(grepl("Input is either empty", file[1])){
        file <- read.table(text = "", col.names = cols, colClasses = cols.class)
      }else{
        stop("Error in loading key files. Check input files.")
      }
      
    }
    
    return( file )
    
  }else if( format == "tsv" ){
    
    file <- try(
      data.table::fread(key.file, sep = "\t", data.table = FALSE), silent = TRUE
    )
    
    if(any(class(file) == "try-error")){
      
      if(grepl("Input is either empty", file[1])){
        file <- read.table(text = "", col.names = cols, colClasses = cols.class)
      }else{
        stop("Error in loading key files. Check input files.")
      }
      
    }
    
    return( file )
    
  }else if( format == "rds" ){
    
    file <- readRDS(key.file)
    return(as.data.frame(file))
    
  }else if( format == "RData" ){
    
    env <- new.env()
    load(key.file, envir = env)
    if( length(env) > 1 ){
      stop("More than one object present in key.file (RData).")
    }
    
    return( env[[ls(env)]] )
    
  }else{
    
    stop(
      "Key file not in a supported file type,",
      " convert to csv, tsv, rds, or RData."
    )
    
  }
  
}
