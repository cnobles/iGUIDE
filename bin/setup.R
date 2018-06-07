# Options
options(stringsAsFactors = FALSE)
.libPaths(new = grep("miniconda3/env/", .libPaths(), value = TRUE))

r_libs <- c(
  "argparse", "data.table", "devtools", "igraph", "knitr", "magrittr", "Matrix",
  "pander", "RColorBrewer", "reshape2", "rmarkdown", "scales", "tidyverse", 
  "yaml")
bioc_libs <- c(
  "BiocGenerics", "Biostrings", "BSgenome", "BSgenome.Hsapiens.UCSC.hg38",
  "GenomicRanges", "hiAnnotator", "IRanges", "Rsamtools", "ShortRead")

# Install developmental based packages from github
options(unzip = "internal")
devtools::install_github(
  "cnobles/gintools",
  repos = "https://mran.microsoft.com/snapshot/2018-05-01/",
  upgrade_dependencies = FALSE,
  dependencies = FALSE)

all_required <- c(r_libs, bioc_libs, "gintools")
is_installed <- suppressMessages(
  sapply(all_required, require, character.only = TRUE))

print(data.frame(
  "Loaded" = is_installed,
  "R-Package" = names(is_installed)),
  right = FALSE, row.names = FALSE)

q()

if(!all(is_installed)){
  # Install required R-packages, downloaded from CRAN Mirror snapshot
  get_r_libs <- r_libs[!r_libs %in% row.names(installed.packages())]
  if(length(get_r_libs) > 0){
    install.packages(
      get_r_libs,
      repos = "https://mran.microsoft.com/snapshot/2018-05-01/",
      dependencies = c("Depends", "Imports")) }
  
  # Install BioConductoR-based packages
  suppressMessages(source("https://bioconductor.org/biocLite.R"))
  get_bioc_libs <- bioc_libs[!bioc_libs %in% row.names(installed.packages())]
  if(length(get_bioc_libs) > 0){
    biocLite(
      get_bioc_libs,
      suppressUpdates = TRUE, ask = FALSE,
      siteRepos = "https://mran.microsoft.com/snapshot/2018-05-01/") }
  
  # Install developmental based packages from github
  options(unzip = "internal")
  devtools::install_github(
    "cnobles/gintools",
    repos = "https://mran.microsoft.com/snapshot/2018-05-01/",
    upgrade_dependencies = FALSE,
    dependencies = FALSE)
  
  # Check for installed packages
  all_required <- c(r_libs, bioc_libs, "gintools")
  is_installed <- suppressMessages(
    sapply(all_required, require, character.only = TRUE))
  
  print(data.frame(
    "Loaded" = is_installed,
    "R-Package" = names(is_installed)),
    right = FALSE, row.names = FALSE)
  
  is_installed <- suppressMessages(
    sapply(all_required, require, character.only = TRUE))
  
  if(!all(is_installed)){
    stop("Not all required R-packages have been installed. Check dependencies.")
  }else{
    message("\nAll required packages installed from extra sources.")
  }
}else{
  message("\nAll required packages installed from conda.")
}

q()