# This R script is for checking for missing dependencies and for installing
# those missing through CRAN or a specified mirror. Options should be passed
# as follows:
#   Usage: Rscript path/to/check_for_required_packages [options]
#   Options:
#     --cran               Install missing dependencies from CRAN.
#     --cran_mirror [url]  Install missing dependencies from CRAN mirror.
#     --conda              Install within an active conda environment.
#     -q, --quiet          Execute script in a non-verbose mode.

# Read command line inputs ----
options(stringsAsFactors = FALSE, unzip = "internal")

cmd_args <- unlist(strsplit(c("", commandArgs(trailingOnly = TRUE)), " "))


# Set arguments ----
cran_install <- any(grepl("--cran$", cmd_args, perl = TRUE))

mirror_install <- any(grepl("--cran_mirror$", cmd_args, perl = TRUE))

mirror_url <- cmd_args[
  which(grepl("--cran_mirror$", cmd_args, perl = TRUE)) + 1
]

within_conda <- any(grepl("--conda$", cmd_args, perl = TRUE)) 

quiet <- any(grepl("-q", cmd_args))


# Check installed packages for dependencies ----
r_packs <- c(
  "argparse", "data.table", "devtools", "digest", "igraph", "ggforce", 
  "knitr", "magrittr", "Matrix", "pander", "RColorBrewer", "rmarkdown", 
  "scales", "tidyverse", "yaml")

bioc_packs <- c(
  "BiocGenerics", "Biostrings", "BSgenome", "BSgenome.Hsapiens.UCSC.hg38",
  "GenomicRanges", "hiAnnotator", "IRanges", "Rsamtools", "ShortRead"
)

packs <- c(r_packs, bioc_packs)

present <- packs %in% row.names(installed.packages())

if( !quiet ){
  print(data.frame(row.names = packs, "Installed" = present))
}

if( !cran_install | !mirror_install ){
  stopifnot(all(present))
  q()
}


# Install from CRAN or from CRAN mirror ----

if( within_conda ){

  .libPaths(
    new = grep(
      pattern = "conda./envs/", x = .libPaths(), 
      perl = TRUE, value = TRUE
    )
  )

}

if( mirror_install ){
  repo <- mirror_url
}else{
  repo <- getOption("repos")
}

r_packs_to_get <- r_packs[!r_packs %in% row.names(installed.packages())]

if( length(r_packs_to_get) > 0 ){

  install.packages(
    r_packs_to_get, 
    repos = repo, 
    dependencies = c("Depends", "Imports"),
    quiet = TRUE
  )

}


# Install from BioConductor ----
bioc_packs_to_get <- bioc_packs[
  !bioc_packs %in% row.names(installed.packages())
]

if( length(bioc_packs_to_get) > 0 ){

  suppressMessages(source("https://bioconductor.org/biocLite.R"))

  biocLite(
    bioc_packs_to_get,
    suppressUpdates = TRUE, 
    ask = FALSE,
    siteRepos = repo
  )

}


# Check for installed packages again and close out

if( !quiet ){
  print(data.frame(row.names = packs, "Installed" = present))
}

stopifnot(all(present))

q()

