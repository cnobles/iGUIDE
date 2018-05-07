#' For anyone reviewing the code below, the following is a small style guide 
#' outlining the various formats for the code. 
#' 
#' Names with "_": objects, inlucding data.frames, GRanges, vectors, ...
#' Names in caMel format: functions or components of objects (i.e. columns 
#' within a data.frame).
#' Names with ".": arguments / options for functions

options(stringsAsFactors = FALSE, scipen = 99)
suppressMessages(library("argparse"))
suppressMessages(library("pander"))
suppressMessages(library("yaml"))
panderOptions("table.style", "simple")
panderOptions("table.split.table", Inf)


# Set up and gather command line arguments -------------------------------------
parser <- ArgumentParser(
  description = "Generate an iGUIDE report for input run(s).")
parser$add_argument(
  "input", nargs = "+", type = "character",
  help = "Post-processing output .rds file(s).")
parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", 
  help = "Output report file, pdf format.")
parser$add_argument(
  "-c", "--config", nargs = "+", type = "character",
  help = "Run specific config file(s) in yaml format.")
parser$add_argument(
  "-s", "--support", nargs = "+", type = "character",
  help = "Supplementary data input, csv or tsv format.")
parser$add_argument(
  "-f", "--figures", action = "store_true",
  help = "Figures will be written to a directory with output.")
parser$add_argument(
  "-d", "--data", action = "store_true",
  help = "Data to generate the report will be saved as an R image with output.")


args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

if(length(args$input) != length(args$config)){
  stop("Must supply one config file for each input.")
}

input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(1:length(args), function(i){
    paste(args[[i]], collapse = ", ")}))
input_table <- input_table[
  match(
    c("input :", "output :", "config :", "support :", "figures :", "data :"),
    input_table$Variables),]
pandoc.title("iGUIDE Report Inputs:")
pandoc.table(data.frame(input_table, row.names = NULL), 
             justify = c("left", "left"), 
             split.tables = Inf)


# Load dependancies ------------------------------------------------------------
add_packs <- c(
  "stringr", "magrittr", "dplyr", "data.table", "Biostrings", "GenomicRanges",
  "knitr", "markdown")

add_packs_loaded <- suppressMessages(
  sapply(add_packs, require, character.only = TRUE))
if(!all(add_packs_loaded)){
  pandoc.table(data.frame(
    "R-Packages" = names(add_packs_loaded), 
    "Loaded" = add_packs_loaded, 
    row.names = NULL))
  stop("Check dependancies.")
}

code_dir <- dirname(sub(
  "--file=", "", 
  grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))


# Import metadata and consolidate into report objects --------------------------
## Load config files
configs <- lapply(args$config, yaml.load_file)
names(configs) <- sapply(configs, "[[", "Run_Name")

umitag_option <- all(unlist(lapply(configs, "[[", "UMItags")))

## Combine sampleInfo files
sample_info <- bind_rows(lapply(
    paste0(
      sapply(configs, "[[", "Install_Directory"), "/", 
      sapply(configs, "[[", "Sample_Info")), 
    function(x){
      data.table::fread(x, data.table = FALSE)
    }), .id = "run_set")
sample_name_col <- unique(sapply(configs, "[[", "Sample_Name_Column"))
if(length(sample_name_col) != 1) stop("SampleInfo files not in same format.")
sample_info$specimen <- stringr::str_extract(
  sample_info[,sample_name_col], "[\\w]+")

## Identify all gRNAs used from config files
gRNAs <- lapply(
  do.call(c, lapply(configs, "[[", "Guide_RNA_Sequences")), toupper)
pams <- lapply(
  do.call(c, lapply(configs, "[[", "PAM_Sequence")), toupper)

gRNAs <- split(gRNAs, stringr::str_extract(names(gRNAs), "[\\w\\-\\_]+"))
pams <- split(pams, stringr::str_extract(names(pams), "[\\w\\-\\_]+"))

gRNAs <- bind_rows(mapply(function(seqs, pam){
    pam_mat <- sapply(pam, function(pat){
      stringr::str_detect(seqs, paste0(pat,"$")) })
    if(any(rowSums(pam_mat) > 1)){ 
      stop("Multiple PAM sequences detected on a single guide RNA.") }
    rownames(pam_mat) <- str_extract(names(seqs), "[\\w\\-\\_]+$")
    colnames(pam_mat) <- unlist(pam)
    data.frame(
      row.names = rownames(pam_mat),
      "Guide" = rownames(pam_mat),
      "gRNA" = sapply(1:length(seqs), function(i){
        stringr::str_replace(
          seqs[[i]], paste0(colnames(pam_mat)[pam_mat[i,]], "$"), "")}),
      "PAM" = colnames(pam_mat)[
        sapply(1:nrow(pam_mat), function(i) which(pam_mat[i,]))])
  }, seqs = gRNAs, pam = pams, SIMPLIFY = FALSE), .id = "run.set") %>%
  distinct(Guide, gRNA, PAM)


## Identify on-target edit sites from config files
on_targets <- unlist(lapply(configs, "[[", "On_Target_Sites"))
names(on_targets) <- stringr::str_extract(names(on_targets), "[\\w\\_\\-]+$")
on_targets <- on_targets[match(unique(names(on_targets)), names(on_targets))]

## Treatment across runs
treatments <- lapply(configs, "[[", "Treatment")

if(any(grepl("sampleInfo:", treatments[[1]]))){
  info_col <- match(
    stringr::str_extract(treatments[[1]], "[\\w]+$"), names(sample_info))
  if(length(info_col) != 1){
    stop("Cannot parse treatment data. Check config yaml and sampleInfo.")
  }
  treatment_df <- data.frame(
    sampleName = sample_info[,sample_name_col], 
    treatment = sample_info[,info_col])
  treatment_df$specimen <- str_extract(treatment_df$sampleName, "[\\w]+")
  treatment_df <- unique(treatment_df[,c("specimen", "treatment")]) %>%
    arrange(specimen)
  treatment <- strsplit(treatment_df$treatment, ";")
  names(treatment) <- treatment_df$specimen
}else if(any(grepl("all", names(treatments[[1]])))){
  treatment_df <- data.frame(
    sampleName = sample_info[,sample_name_col], 
    treatment = unique(unlist(treatments)))
  treatment_df$specimen <- str_extract(treatment_df$sampleName, "[\\w]+")
  treatment_df <- unique(treatment_df[,c("specimen", "treatment")]) %>%
    arrange(specimen)
  treatment <- strsplit(treatment_df$treatment, ";")
  names(treatment) <- treatment_df$specimen
}

## Load in supporting information ==============================================
if(length(args$support) > 0){
  supp_data <- data.table::fread(args$support, data.table = FALSE)
  supp_data <- filter(supp_data, specimen %in% sample_info$specimen)
}else{
  supp_data <- data.frame()
}

# Read in experimental data and contatenate different sets ---------------------
input_data <- lapply(args$input, readRDS)
names(input_data) <- names(configs)
data_type_names <- names(input_data[[1]])
input_data <- lapply(1:length(input_data[[1]]), function(i){
  bind_rows(lapply(input_data, "[[", i), .id = "run.set")
})
names(input_data) <- data_type_names

# Data passed to Rmd for report generation -------------------------------------
set_names <- ifelse(
  length(configs) ==  1, 
  names(configs), 
  paste0(paste(names(configs)[1:(length(configs)-1)], collapse = ", "), 
         ", and ", names(configs)[length(configs)]))

# Normalize file output path
write(c(), file = args$output)
args$output <- normalizePath(args$output)
unlink(args$output)

output_path <- unlist(strsplit(args$output, "/"))
output_dir <- paste(output_path[1:(length(output_path)-1)], collapse = "/")
output_file <- output_path[length(output_path)]

if(!str_detect(output_file, ".pdf$")){
  output_file <- paste0(output_file, ".pdf")
}

if(args$data){
  save.image(file = file.path(
    output_dir, stringr::str_replace(output_file, ".pdf$", ".RData")))
}

rmarkdown::render(
  input = file.path(code_dir, "iGUIDE_report_template.Rmd"),
  output_format = "all", 
  output_file = output_file,
  output_dir = output_dir)

q()
