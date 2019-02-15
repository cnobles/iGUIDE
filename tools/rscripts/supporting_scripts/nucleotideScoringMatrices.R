#' A Binary Ambiguous Nucleotide scoring Matrix (BAN Mat)
#' 
#' @usage banmat()
#' 
#' @description Meant for comparing ambiguous sequences against "A", "T", "G", 
#' "C", and "N" containing sequences. Currently matches between ambiuous 
#' nucleotides are considered mismatch. Constructed based on NUC4.4. Function
#' outputs a single square matrix.
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

#' A Symmetric Scaled Ambiguous Nucleotide scoring Matrix (SSAN Mat)
#' 
#' @usage ssanmat()
#' 
#' @description Meant for comparing ambiguous sequences against other sequences
#' possibly containing ambiguous nucleotides. The scores are scaled between 0
#' and 1 for ambiguous matches based on overlaping proportions and the matrix is
#' symmetric. Constructed originally from NUC4.4. Function outputs a single 
#' square matrix.
#' 
#' @author Christopher Nobles, Ph.D.

ssanmat <- function(){
  matrix(c(
    1,0,0,0, 0   , 1   , 1   , 0   , 0   , 1   , 0   , 1   , 1   , 1   , 1,
    0,1,0,0, 0   , 1   , 0   , 1   , 1   , 0   , 1   , 0   , 1   , 1   , 1,
    0,0,1,0, 1   , 0   , 1   , 0   , 1   , 0   , 1   , 1   , 0   , 1   , 1,
    0,0,0,1, 1   , 0   , 0   , 1   , 0   , 1   , 1   , 1   , 1   , 0   , 1,
    0,0,1,1, 1   , 0   , 0.5 , 0.5 , 0.5 , 0.5 , 0.67, 0.67, 0.33, 0.33, 1,
    1,1,0,0, 0   , 1   , 0.5 , 0.5 , 0.5 , 0.5 , 0.33, 0.33, 0.67, 0.67, 1,
    1,0,1,0, 0.5 , 0.5 , 1   , 0   , 0.5 , 0.5 , 0.33, 0.67, 0.33, 0.67, 1,
    0,1,0,1, 0.5 , 0.5 , 0   , 1   , 0.5 , 0.5 , 0.67, 0.33, 0.67, 0.33, 1,
    0,1,1,0, 0.5 , 0.5 , 0.5 , 0.5 , 1   , 0   , 0.67, 0.33, 0.33, 0.67, 1,
    1,0,0,1, 0.5 , 0.5 , 0.5 , 0.5 , 0   , 1   , 0.33, 0.67, 0.67, 0.33, 1,
    0,1,1,1, 0.67, 0.33, 0.33, 0.67, 0.67, 0.33, 1   , 0.67, 0.67, 0.67, 1,
    1,0,1,1, 0.66, 0.33, 0.67, 0.33, 0.33, 0.67, 0.67, 1   , 0.67, 0.67, 1,
    1,1,0,1, 0.33, 0.67, 0.33, 0.67, 0.33, 0.67, 0.67, 0.67, 1   , 0.67, 1,
    1,1,1,0, 0.33, 0.67, 0.67, 0.33, 0.67, 0.33, 0.67, 0.67, 0.67, 1   , 1,
    1,1,1,1, 1   , 1   , 1   , 1   , 1   , 1   , 1   , 1   , 1   , 1   , 1),
    ncol = 15,
    nrow = 15,
    byrow = TRUE,
    dimnames = list(
      c("A", "T", "G", "C", "S", "W", "R", "Y", 
        "K", "M", "B", "V", "H", "D", "N"),
      c("A", "T", "G", "C", "S", "W", "R", "Y", 
        "K", "M", "B", "V", "H", "D", "N")))
}

#' A Unsymmetric Scaled Ambiguous Nucleotide scoring Matrix (USAN Mat)
#' 
#' @usage usanmat()
#' 
#' @description Meant for comparing ambiguous sequences against other sequences
#' possibly containing ambiguous nucleotides. The scores are scaled between 0
#' and 1 for ambiguous matches based on overlaping proportions and the matrix is
#' unsymmetric. This is meant for proper scoreing between different ambiguous
#' bases. For example, comparing a "S" (meaning a "C" or "G") to a "V" (meaning
#' a "C", "G", or "A") would be if either nucleotide was chosen (score = 1), but
#' comparing a "V" to an "S" would only be correct 2 out of three times 
#' (score = 0.67).Constructed originally from NUC4.4. Function outputs a single 
#' square matrix.
#' 
#' @author Christopher Nobles, Ph.D.

usanmat <- function(){
  matrix(c(
    # A T G C  S     W     R     Y     K     M     B     V     H     D     N
    1,0,0,0, 0   , 1   , 1   , 0   , 0   , 1   , 0   , 1   , 1   , 1   , 1, #A
    0,1,0,0, 0   , 1   , 0   , 1   , 1   , 0   , 1   , 0   , 1   , 1   , 1, #T
    0,0,1,0, 1   , 0   , 1   , 0   , 1   , 0   , 1   , 1   , 0   , 1   , 1, #G
    0,0,0,1, 1   , 0   , 0   , 1   , 0   , 1   , 1   , 1   , 1   , 0   , 1, #C
    0,0,1,1, 1   , 0   , 0.5 , 0.5 , 0.5 , 0.5 , 1   , 1   , 0.5 , 0.5 , 1, #S
    1,1,0,0, 0   , 1   , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 , 1   , 1   , 1, #W
    1,0,1,0, 0.5 , 0.5 , 1   , 0   , 0.5 , 0.5 , 0.5 , 1   , 0.5 , 1   , 1, #R
    0,1,0,1, 0.5 , 0.5 , 0   , 1   , 0.5 , 0.5 , 1   , 0.5 , 1   , 0.5 , 1, #Y
    0,1,1,0, 0.5 , 0.5 , 0.5 , 0.5 , 1   , 0   , 1   , 0.5 , 0.5 , 1   , 1, #K
    1,0,0,1, 0.5 , 0.5 , 0.5 , 0.5 , 0   , 1   , 0.5 , 1   , 1   , 0.5 , 1, #M
    0,1,1,1, 0.67, 0.33, 0.33, 0.67, 0.67, 0.33, 1   , 0.67, 0.67, 0.67, 1, #B
    1,0,1,1, 0.66, 0.33, 0.67, 0.33, 0.33, 0.67, 0.67, 1   , 0.67, 0.67, 1, #V
    1,1,0,1, 0.33, 0.67, 0.33, 0.67, 0.33, 0.67, 0.67, 0.67, 1   , 0.67, 1, #H
    1,1,1,0, 0.33, 0.67, 0.67, 0.33, 0.67, 0.33, 0.67, 0.67, 0.67, 1   , 1, #D
    1,1,1,1, 1   , 1   , 1   , 1   , 1   , 1   , 1   , 1   , 1   , 1   , 1),#N
    ncol = 15,
    nrow = 15,
    byrow = TRUE,
    dimnames = list(
      c("A", "T", "G", "C", "S", "W", "R", "Y", 
        "K", "M", "B", "V", "H", "D", "N"),
      c("A", "T", "G", "C", "S", "W", "R", "Y", 
        "K", "M", "B", "V", "H", "D", "N")))
}