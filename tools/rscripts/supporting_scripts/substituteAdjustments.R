#' @param seqs sequences coercible to character vectors. 
#' @param adjs character vector where each value is two characters long, the 
#' first is the character to match, the second the character to substitute. For
#' example, the input of 'AR' refers to the adjustment of substituting 'R' for 
#' all 'A' characters.

substituteAdjustments <- function(seqs, adjs, ...){
  
  # Stop if input format is wrong
  if( any(nchar(adjs) != 2) ){
    stop(
      "\n  Incorrect input format for substituteAdjustments. Please read ",
      "documentation."
    )
  }
  
  # Create find and substitute vectors
  adj_f <- stringr::str_extract(adjs, "[\\w]")
  adj_s <- stringr::str_extract(adjs, "[\\w]$")
  
  # Create new environment for modifying input sequences
  mod_env <- new.env()
  mod_env$seqs <- as.character(seqs)
  
  # Iteratively change the sequence given adjustments
  null <- mapply(
    function(f, s){
      mod_env$seqs <- gsub(f, s, mod_env$seqs, fixed = TRUE)
    },
    f = adj_f,
    s = adj_s,
    SIMPLIFY = FALSE
  )
  
  # Output the modified sequences
  mod_env$seqs
  
}