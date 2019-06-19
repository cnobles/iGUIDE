# This script writes / generates a short report with the stat.csv file generated
# from processing a run.
# 
# usage: 
#   Rscript write_stat_report <stat.csv> -o <output> -c <config> -f <format(html/pdf)> -t <template.path> -i <iguide_dir>
# 
# Input arguments for this script include: 
#   input stat.csv file 
#   output file name / path
#   format, either 'html' (default) or 'pdf'
#   iguide installation directory, will look for sys argument 'IGUIDE_DIR'
# 

options(stringsAsFactors = FALSE, scipen = 99, width = 120)

args <- commandArgs(trailingOnly = TRUE)

code_dir <- dirname(sub(
  pattern = "--file=", 
  replacement = "", 
  x = grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
))

# Check input file ----
core_file <- args[1]
eval_file <- args[2]

if( !file.exists(core_file) | !file.exists(eval_file) ){
  stop("\n  Cannot find input stat files. Check inputs.")
}

# Check output file ----
if( length(grep("-o", args)) > 0){
  output_file <- args[grep("-o", args) + 1]
}else{
  stop("\n  Output file is required, please specify with '-o'.\n")
}

if( file.exists(output_file) ){
  cat("Removing existing output file of the same name: ", args[2], "\n")
  unlink(output_file)
}


# Check config input ----
if( length(grep("-c", args)) > 0 ){
  config_path <- args[grep("-c", args) + 1]
}else{
  stop("\n  Config input is required, please specify with '-c'.\n")
}

if( !file.exists(config_path) ){
  stop("\n  Cannot find config file: ", config_path, ".\n")
}

config <- yaml::yaml.load_file(config_path)

# Check for format ----
if( length(grep("-f", args)) > 0 ){
  output_format <- args[grep("-f", args) + 1]
}else{
  output_format <- "html"
}

report_formats <- c("html" = "html_document", "pdf" = "pdf_document")

if( !output_format %in% names(report_formats) ){
  stop(
    "\n  Please input either 'html' or 'pdf' for format.\n",
    "  Other formats not supported."
  )
}

output_format <- report_formats[output_format]

# Check for install dir ----
if( length(grep("-i", args)) > 0 ){
  iguide_dir <- args[grep("-i", args) + 1]
}else{
  iguide_dir <- Sys.getenv("IGUIDE_DIR")
}

if( !dir.exists(iguide_dir) ){
  stop(paste0("\n  Cannot find install path to iGUIDE: ", iguide_dir, ".\n"))
}
  
# Check for template path ----
if( length(grep("-t", args)) > 0 ){
  template_path <- args[grep("-t", args) + 1]
}else{
  template_path <- file.path(
    iguide_dir, "tools/rscripts/report_templates/iGUIDE_stat_template.Rmd"
  )
}

## Resolve template file path.
if( file.exists(file.path(iguide_dir, template_path)) ){
  
  template_path <- normalizePath(file.path(iguide_dir, template_path))
  
}else if( file.exists(file.path(template_path)) ){
  
  template_path <- normalizePath(file.path(template_path))
  
}else{
  
  stop("\n  Cannot find template file: ", template_path, ".\n")
  
}


# Load required r-packages ----
packs_loaded <- sapply(
  c("magrittr", "knitr"), require, character.only = TRUE
)

if( !any(packs_loaded) ){
  stop(
    "\n  Could not find required r-package: ", 
    paste(names(packs_loaded)[packs_loaded], collapse = ", "), 
    ".\n"
  )
}


# Get versioning ----
soft_version <- as.character(read.delim(
  file = file.path(iguide_dir, ".version"), header = FALSE
))

build_version <- list.files(file.path(iguide_dir, "etc")) %>%
  grep(pattern = "build.b[0-9\\.]+.*", x = ., value = TRUE) %>%
  stringr::str_extract(pattern = "b[0-9]+\\.[0-9]+\\.[0-9]+")

signature <- config[["signature"]]

# Load input data ----
core_stat_df <- read.csv(core_file)
eval_stat_df <- read.csv(eval_file)

stat_df <- dplyr::full_join(core_stat_df, eval_stat_df, by = "sampleName") %>%
  dplyr::mutate_all(function(x) ifelse(is.na(x), rep(0, length(x)), x))

sampleName_levels <- unique(stat_df$sampleName)

sampleNames <- sampleName_levels[
  -match(
    c("ambiguous_reads", "degenerate_reads", "unassigned_reads"), 
    sampleName_levels
  )
]

sampleName_levels <- c(
  sampleNames, "ambiguous_reads", "degenerate_reads", "unassigned_reads"
)


# Read attrition table ----
read_tbl <- dplyr::select(
    stat_df, sampleName, demulti.reads, R1.trim.reads, R2.primer.trim.reads, 
    R2.trim.reads, umitags.reads, filt.reads, R1.consol.reads, R2.consol.reads, 
    align.unique.reads, align.chimera.reads, align.multihit.reads
  ) %>%
  dplyr::mutate(sampleName = factor(sampleName, levels = sampleName_levels)) %>%
  dplyr::arrange(sampleName)

names(read_tbl) <- stringr::str_replace(names(read_tbl), ".reads$", "")


# Alignment outcome table ----
algn_tbl <- dplyr::select(
    stat_df, sampleName, align.unique.reads, align.unique.algns, 
    align.unique.loci, align.multihit.reads, align.multihit.lengths, 
    align.multihit.clusters, align.chimera.reads
  ) %>%
  dplyr::filter(sampleName %in% sampleNames) %>%
  dplyr::mutate(sampleName = factor(sampleName, levels = sampleNames)) %>%
  dplyr::arrange(sampleName)

names(algn_tbl) <- stringr::str_replace(names(algn_tbl), "align.", "")


# Incorporation breakdown table ----
incorp_tbl <- dplyr::select(
    stat_df, sampleName, eval.total.algns, eval.combined.algns, 
    eval.pileup.algns, eval.paired.algns, eval.matched.algns, 
    eval.ontarget.algns, eval.offtarget.algns
  ) %>%
  dplyr::filter(sampleName %in% sampleNames) %>%
  dplyr::mutate(sampleName = factor(sampleName, levels = sampleNames)) %>%
  dplyr::arrange(sampleName)

names(incorp_tbl) <- stringr::str_replace(names(incorp_tbl), "eval.", "")
names(incorp_tbl) <- stringr::str_replace(names(incorp_tbl), ".algns$", "")


# Initiate report generation ----
# Normalize file output path
write(c(), file = output_file)
output_file <- normalizePath(output_file)
unlink(output_file)

output_path <- unlist(strsplit(output_file, "/"))
output_dir <- paste(output_path[seq_len(length(output_path)-1)], collapse = "/")
output_file <- output_path[length(output_path)]

if( output_format == "html_document" & !stringr::str_detect(output_file, ".html$") ){
  output_file <- paste0(output_file, ".html")
}

if( output_format == "pdf_document" & !stringr::str_detect(output_file, ".pdf$") ){
  output_file <- paste0(output_file, ".pdf")
}

if( output_format == "html_document" ){
  
  css_path <- normalizePath(file.path(code_dir, "report_templates/iguide.css"))
  
  rmarkdown::render(
    input = template_path,
    output_format = output_format, 
    output_file = output_file,
    output_dir = output_dir,
    output_options = list("css" = css_path)
  )
  
}else{
  
  rmarkdown::render(
    input = template_path,
    output_format = output_format, 
    output_file = output_file,
    output_dir = output_dir
  )
  
}

q()
