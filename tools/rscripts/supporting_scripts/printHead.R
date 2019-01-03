#' Pander out head portions of objects for checking.
#' @param object object which will be coerced into a data.frame and will be 
#' passed to head().
#' @param title character Title to be included above table
#' @param caption character Caption to be included beneath table.
#' @param row.names logical TRUE includes while FALSE (default) drops row names.
#' @author Christopher Nobles, Ph.D.
#' 

printHead <- function(object, title = NULL, caption = NULL, row.names = FALSE){
  
  stopifnot( !class(object) == "list" )
  
  if( !is.null(title) ) cat(paste0("\n", title, ":"))
  
  if(!row.names){
    df <- as.data.frame( head(object), row.names = NULL )
  }else{
    df <- as.data.frame( head(object) )
  }
  
  acceptable_classes <- c("factor", "character", "numeric", "integer")
  
  df <- df[, which(
    sapply( seq_len(ncol(df)), function(i) class(df[,i]) ) %in% 
      acceptable_classes
  )]
  
  print(df, row.names = FALSE)
  cat(paste0("Table Caption: ", caption))
    
}
