#!/usr/bin/env Rscript

# Check for installed packages
options(stringsAsFactors = FALSE, scipen = 99)

# Capture commandline files
parser <- argparse::ArgumentParser(
  description = "Script to check for an installed package.",
  usage = "Rscript tools/rscripts/check_pkgs.R <pkgs>"
)

parser$add_argument(
  "pkg", nargs = "+", type = "character", default = "NA",
  help = "Package(s) name."
)

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

pkgs <- args$pkg

pkgs_present <- pkgs %in% rownames(installed.packages())

if( all(pkgs_present) ){
  
  q(save = "no", status = 0)
  
}else{
  
  cat(
    " Packages not installed:\n  ", 
    paste(pkgs[!pkgs_present], collapse = "\n   "), 
    "\n"
  )
  
  q(save = "no", status = 1)
  
}
