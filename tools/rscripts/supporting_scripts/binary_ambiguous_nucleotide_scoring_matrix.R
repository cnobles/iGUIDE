#' A Binary Ambiguous Nucleotide scoring Matrix (BAN Mat)
#' 
#' Constructed based on NUC4.4.
#' 
#' Meant for comparing ambiguous sequences against "A", "T", "G", "C", and "N"
#' containing sequences. Currently matches between ambiuous nucleotides are 
#' considered mismatch.
#' 
#' @author Christopher Nobles, Ph.D.

banmat <- function(){
  matrix(c(
    1,0,0,0,0,1,1,0,0,1,0,1,1,1,1,
    0,1,0,0,0,1,0,1,1,0,1,0,1,1,1,
    0,0,1,0,1,0,1,0,1,0,1,1,0,1,1,
    0,0,0,1,1,0,0,1,0,1,1,1,1,0,1,
    0,0,1,1,1,0,0,0,0,0,0,0,0,0,1,
    1,1,0,0,0,1,0,0,0,0,0,0,0,0,1,
    1,0,1,0,0,0,1,0,0,0,0,0,0,0,1,
    0,1,0,1,0,0,0,1,0,0,0,0,0,0,1,
    0,1,1,0,0,0,0,0,1,0,0,0,0,0,1,
    1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,
    0,1,1,1,0,0,0,0,0,0,1,0,0,0,1,
    1,0,1,1,0,0,0,0,0,0,0,1,0,0,1,
    1,1,0,1,0,0,0,0,0,0,0,0,1,0,1,
    1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
    ncol = 15,
    nrow = 15,
    byrow = TRUE,
    dimnames = list(
      c("A", "T", "G", "C", "S", "W", "R", "Y", 
        "K", "M", "B", "V", "H", "D", "N"),
      c("A", "T", "G", "C", "S", "W", "R", "Y", 
        "K", "M", "B", "V", "H", "D", "N")))
}
