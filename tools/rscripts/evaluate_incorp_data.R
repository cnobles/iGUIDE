#' For anyone reviewing the code below, the following is a small style guide 
#' outlining the various formats for the code. 
#' 
#' Names with "_": objects, inlucding data.frames, GRanges, vectors, ...
#' Names in caMel format: functions or components of objects (i.e. columns 
#' within a data.frame).
#' Names with ".": arguments / options for functions

options(stringsAsFactors = FALSE, scipen = 99, width = 180)


# Set up and gather command line arguments ----
parser <- argparse::ArgumentParser(
  description = "Evaluation of iGUIDE data from input run(s).",
  usage = paste(
    "iguide eval <config(s)> -o <output> [-h/--help, -v/--version]",
    "[optional args]"
  )
)

parser$add_argument(
  "config", nargs = "+", type = "character",
  help = paste(
    "Run specific config file(s) in yaml format. Can specify more than",
    "one to combine several runs together for evaluation."
  )
)

parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", required = TRUE,
  help = "Output eval file, .rds format. i.e. output.rds or output"
)

parser$add_argument(
  "-s", "--support", nargs = 1, type = "character",
  help = paste(
    "Supplementary data input, csv or tsv format. Only one file. Must have",
    "'specimen' column and only specimens matching data in this column will",
    "be considered for evaluation."
  )
)

parser$add_argument(
  "-q", "--quiet", action = "store_true", 
  help = "Hide standard output messages."
)

parser$add_argument(
  "--iguide_dir", nargs = 1, type = "character", default = "IGUIDE_DIR",
  help = "iGUIDE install directory path, do not change for normal applications."
)


args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

if( !dir.exists(args$iguide_dir) ){
  root_dir <- Sys.getenv(args$iguide_dir)
}else{
  root_dir <- args$iguide_dir
}

if( !dir.exists(root_dir) ){
  stop(paste0("\n  Cannot find install path to iGUIDE: ", root_dir, ".\n"))
}else{
  args$iguide_dir <- root_dir
}

## Determine output file name and path
if( !stringr::str_detect(args$output, ".rds$") ){
  args$output <- paste0(args$output, ".rds")
}

write(c(), file = args$output)
args$output <- normalizePath(args$output)
unlink(args$output)

## Construct input table and print to terminal
input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(seq_along(args), function(i){
    paste(args[[i]], collapse = ", ")
  })
)

input_table <- input_table[
  match(
    c("config :", "output :", "support :", "iguide_dir :"),
    input_table$Variables),
]

if( !args$quiet ){
  
  cat("\niGUIDE Evaluation Inputs:\n")
  
  print(
    data.frame(input_table),
    right = FALSE, 
    row.names = FALSE
  )

}

# Load dependancies ----
if( !args$quiet ) cat("\nLoading dependencies.\n")

add_packs <- c("magrittr", "knitr")

add_packs_loaded <- suppressMessages(
  sapply(add_packs, require, character.only = TRUE)
)

if( !all(add_packs_loaded) ){
  
  print(
    data.frame(
      "R-Packages" = names(add_packs_loaded), 
      "Loaded" = add_packs_loaded
    ), 
    right = FALSE,
    row.names = FALSE
  )
  
  stop("Check dependancies.\n")
  
}

code_dir <- dirname(sub(
  pattern = "--file=", 
  replacement = "", 
  x = grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
))

source(file.path(code_dir, "supporting_scripts/iguide_support.R"))


# Import metadata and consolidate objects ----
if( !args$quiet ) cat("Importing experimental data and configurations.\n\n")

## Load config files
configs <- lapply(args$config, function(x){
  if( file.exists(file.path(root_dir, x)) ){
    return(yaml::yaml.load_file(file.path(root_dir, x)))
  }else if( file.exists(x) ){
    return(yaml::yaml.load_file(x))
  }else{
    stop("\n  Cannot find config file: ", x, ".\n")
  }
})

names(configs) <- sapply(configs, "[[", "Run_Name")

## Load reference genome 
if( grepl(".fa", unique(sapply(configs, "[[", "Ref_Genome"))) ){
  
  if( !(
    file.exists(
      file.path(root_dir, unique(sapply(configs, "[[", "Ref_Genome")))
    ) | file.exists(unique(sapply(configs, "[[", "Ref_Genome")))
  ) ){
    stop("\n  Specified reference genome file not found.\n")
  }
  
  ref_file_type <- ifelse(
    grepl(".fastq", unique(sapply(configs, "[[", "Ref_Genome"))), 
    "fastq", 
    "fasta"
  )
  
  if( file.exists(
    file.path(root_dir, unique(sapply(configs, "[[", "Ref_Genome"))) 
    ) ){

    ref_genome <- Biostrings::readDNAStringSet(
      filepath = file.path(
        root_dir, unique(sapply(configs, "[[", "Ref_Genome"))
      ),
      format = ref_file_type
    )
    
  }else{
    
    ref_genome <- Biostrings::readDNAStringSet(
      filepath = unique(sapply(configs, "[[", "Ref_Genome")), 
      format = ref_file_type
    )
  }
  
}else{
  
  ref_genome <- unique(sapply(configs, "[[", "Ref_Genome"))
  
  genome <- grep(
    pattern = ref_genome, 
    x = unique(BSgenome::installed.genomes()), 
    value = TRUE
  )
  
  if( length(genome) == 0 ){
    
    cat("\nInstalled genomes include:")
    print(unique(BSgenome::installed.genomes()))
    cat("\n  Selected reference genome not in list.\n")
    stop("\n  Genome not available.\n")
    
  }else if( length(genome) > 1 ){
    
    cat("\nInstalled genomes include:\n")
    print(unique(BSgenome::installed.genomes()))
    cat(
      "\n  Please be more specific about reference genome.",
      "Multiple matches to input."
    )
    stop("\n  Multiple genomes requested.\n")
    
  }
  
  suppressMessages(library(genome, character.only = TRUE))
  
  ref_genome <- get(genome)
  
}

## Get versioning

soft_version <- as.character(read.delim(
  file = file.path(root_dir, ".version"), header = FALSE))

build_version <- list.files(file.path(root_dir, "etc")) %>%
  grep(pattern = "build.b[0-9\\.]+.*", x = ., value = TRUE) %>%
  stringr::str_extract(pattern = "b[0-9]+\\.[0-9]+\\.[0-9]+")


## Load reference files
ref_genes <- suppressMessages(loadRefFiles(
  configs[[1]]$refGenes, 
  type = "GRanges", 
  freeze = configs[[1]]$Ref_Genome,
  root = root_dir
))

onco_genes <- suppressMessages(loadRefFiles(
  configs[[1]]$oncoGeneList, 
  type = "gene.list", 
  freeze = configs[[1]]$Ref_Genome,
  root = root_dir
))

special_genes <- suppressMessages(loadRefFiles(
  configs[[1]]$specialGeneList, 
  type = "gene.list", 
  freeze = config[[1]]$Ref_Genome,
  root = root_dir
))

umitag_option <- all(unlist(lapply(configs, "[[", "UMItags")))

upstream_dist <- unique(sapply(configs, function(x) x$upstreamDist))
downstream_dist <- unique(sapply(configs, function(x) x$downstreamDist))

if( length(upstream_dist) > 1 | length(downstream_dist) > 1 ){
  stop(
    "\n  Inconsistant upstream or downstream distances between config files.\n",
    "  Comparisons between groups with different run specific criteria\n", 
    "  is not recommended.\n")
}

## Combine sampleInfo files

sample_info <- dplyr::bind_rows(lapply(
    sapply(configs, "[[", "Sample_Info"), 
    function(x){
      
      if( file.exists(file.path(root_dir, x)) ){
        return(data.table::fread(file.path(root_dir, x), data.table = FALSE))
      }else if( file.exists(x) ){
        return(data.table::fread(x, data.table = FALSE))
      }else{
        stop("\n  Cannot find Sample_Info: ", x, ".\n")
      }

    }
  ), 
  .id = "run_set"
)

sample_name_col <- unique(sapply(configs, "[[", "Sample_Name_Column"))

if( length(sample_name_col) != 1 ){
  stop("\n  Sample_Info files not in same format.\n")
}

sample_info$specimen <- stringr::str_extract(
  string = sample_info[,sample_name_col], 
  pattern = "[\\w]+"
)

specimen_levels <- unique(sample_info$specimen)

sample_info$specimen <- factor(sample_info$specimen, levels = specimen_levels)


## Load in supporting information ----
if( length(args$support) > 0 ){
  
  if( file.exists(file.path(root_dir, args$support)) ){
    support_path <- file.path(root_dir, args$support)
  }else if( file.exists(args$support) ){
    support_path <- args$support
  }else{
    stop("\n  Cannot find supporting data file: ", args$support, ".\n")
  }
  
  supp_data <- data.table::fread(support_path, data.table = FALSE) %>%
    dplyr::mutate(run_set = "supp_data")
  
  specimen_levels <- supp_data$specimen[supp_data$specimen %in% specimen_levels]
  
  supp_data <- dplyr::filter(supp_data, specimen %in% specimen_levels) %>%
    dplyr::mutate(specimen = factor(specimen, levels = specimen_levels))
  
  sample_info <- dplyr::filter(sample_info, specimen %in% specimen_levels) %>%
    dplyr::mutate(
      specimen = factor(as.character(specimen), levels = specimen_levels)
    ) %>%
    dplyr::arrange(specimen)
  
}else{
  
  supp_data <- data.frame()
  
}


## Identify on-target edit sites from config files
on_targets <- unlist(lapply(configs, "[[", "On_Target_Sites"))
names(on_targets) <- stringr::str_replace(
  string = names(on_targets), pattern = stringr::fixed("."), replacement = ":"
)

names(on_targets) <- stringr::str_extract(
    string = names(on_targets), 
    pattern = "[\\w\\_\\-\\'\\.]+$"
  ) %>%
  stringr::str_extract(pattern = "[\\w\\_\\-\\.]+")

on_targets <- structure(
  unique(on_targets), 
  names = names(on_targets)[match(unique(on_targets), on_targets)]
)


## Treatment across runs
treatments <- lapply(configs, "[[", "Treatment")

if( any(grepl("sampleInfo:", treatments[[1]])) ){
  
  info_col <- match(
    stringr::str_extract(string = treatments[[1]], pattern = "[\\w]+$"), 
    names(sample_info)
  )
  
  if( length(info_col) != 1 ){
    stop("\n  Cannot parse treatment data. Check config yaml and sampleInfo.\n")
  }
  
  treatment_df <- data.frame(
      run_set = sample_info$run_set,
      sampleName = sample_info[,sample_name_col], 
      treatment = sample_info[,info_col]
    ) %>%
    dplyr::mutate(
      specimen = stringr::str_extract(string = sampleName, pattern = "[\\w]+")
    ) %>%
    dplyr::filter(specimen %in% specimen_levels) %>%
    dplyr::mutate(specimen = factor(specimen, levels = specimen_levels)) %>%
    dplyr::distinct(run_set, specimen, treatment) %>%
    dplyr::arrange(specimen) %>%
    dplyr::mutate(treatment = ifelse(is.na(treatment), "Mock", treatment)) %>%
    dplyr::select(run_set, specimen, treatment)
  
  treatment <- strsplit(treatment_df$treatment, ";")
  names(treatment) <- as.character(treatment_df$specimen)
  
}else if( any(grepl("all", names(treatments[[1]]))) ){
  
  treatment_df <- data.frame(
      run_set = sample_info$run_set,
      sampleName = sample_info[,sample_name_col], 
      treatment = unique(unlist(treatments))
    ) %>%
    dplyr::mutate(
      specimen = stringr::str_extract(string = sampleName, pattern = "[\\w]+")
    ) %>%
    dplyr::filter(specimen %in% specimen_levels) %>%
    dplyr::mutate(specimen = factor(specimen, levels = specimen_levels)) %>%
    dplyr::distinct(run_set, specimen, treatment) %>%
    dplyr::arrange(specimen) %>%
    dplyr::mutate(treatment = ifelse(is.na(treatment), "Mock", treatment)) %>%
    dplyr::select(run_set, specimen, treatment)
  
  treatment <- strsplit(treatment_df$treatment, ";")
  names(treatment) <- treatment_df$specimen
  
}else{
  
  stop(
    "\n  Treatment information not accurately parsed from config(s).\n", 
    "  Check config(s) formating.\n"
  )
  
}


## Nucleases used across runs
nuc_profiles <- unlist(
  unname(lapply(configs, "[[", "Nuclease_Profiles")), 
  recursive = FALSE
)

nuc_profiles <- nuc_profiles[
  match(unique(names(nuc_profiles)), names(nuc_profiles))
]

nucleases <- lapply(configs, "[[", "Nuclease")

if( any(grepl("sampleInfo:", nucleases[[1]])) ){
  
  info_col <- match(
    stringr::str_extract(string = nucleases[[1]], pattern = "[\\w]+$"), 
    names(sample_info)
  )
  
  if( length(info_col) != 1 ){
    stop("Cannot parse nuclease data. Check config yaml and sampleInfo.")
  }
  
  nuclease_df <- data.frame(
    run_set = sample_info$run_set,
    sampleName = sample_info[,sample_name_col], 
    nuclease = sample_info[,info_col]
    ) %>%
    dplyr::mutate(
      specimen = stringr::str_extract(string = sampleName, pattern = "[\\w]+")
    ) %>%
    dplyr::filter(specimen %in% specimen_levels) %>%
    dplyr::mutate(specimen = factor(specimen, levels = specimen_levels)) %>%
    dplyr::distinct(run_set, specimen, nuclease) %>%
    dplyr::arrange(specimen) %>%
    dplyr::mutate(nuclease = ifelse(is.na(nuclease), "Mock", nuclease)) %>%
    dplyr::select(run_set, specimen, nuclease)
  
  nuclease <- strsplit(nuclease_df$nuclease, ";")
  names(nuclease) <- nuclease_df$specimen
  
}else if( any(grepl("all", names(nucleases[[1]]))) ){
  
  nuclease_df <- data.frame(
    run_set = sample_info$run_set,
    sampleName = sample_info[,sample_name_col], 
    nuclease = unique(unlist(nucleases))
    ) %>%
    dplyr::mutate(
      specimen = stringr::str_extract(string = sampleName, pattern = "[\\w]+")
    ) %>%
    dplyr::filter(specimen %in% specimen_levels) %>%
    dplyr::mutate(specimen = factor(specimen, levels = specimen_levels)) %>%
    dplyr::distinct(run_set, specimen, nuclease) %>%
    dplyr::arrange(specimen) %>%
    dplyr::mutate(nuclease = ifelse(is.na(nuclease), "Mock", nuclease)) %>%
    dplyr::select(run_set, specimen, nuclease)
  
  nuclease <- strsplit(nuclease_df$nuclease, ";")
  names(nuclease) <- nuclease_df$specimen
  
}else{
  
  stop(
    "\n  Nuclease information not accurately parsed from config(s).\n", 
    "  Check config(s) formating.\n"
  )
  
}

nuclease_treaments <- dplyr::left_join(
  nuclease_df, treatment_df, by = c("run_set", "specimen")
)

target_combn <- structure(
  strsplit(nuclease_treaments$treatment, ";"), 
  names = as.character(nuclease_treaments$specimen)
)

combn_tbl <- data.frame(
  run_set = nuclease_treaments$run_set[
    as.vector(match(
      S4Vectors::Rle(names(target_combn), lengths(target_combn)), 
      nuclease_treaments$specimen
    ))
  ],
  nuclease = nuclease_treaments$nuclease[
    as.vector(match(
      S4Vectors::Rle(names(target_combn), lengths(target_combn)), 
      nuclease_treaments$specimen
    ))
  ],
  target = unlist(target_combn),
  row.names = NULL
  ) %>%
  dplyr::filter(target != "Mock") %>%
  dplyr::distinct()


## Identify all target sequences used from config files
target_seqs <- lapply(
  do.call(c, lapply(configs, "[[", "Target_Sequences")), 
  toupper
)

target_grps <- stringr::str_extract(
  string = names(target_seqs), 
  pattern = "[\\w\\-\\_]+"
)

names(target_seqs) <- sub("[\\w\\-\\_]+.", "", names(target_seqs), perl = TRUE)
target_seqs <- split(target_seqs, target_grps)

target_seqs_df <- data.frame(
  run_set = as.character(
    S4Vectors::Rle(names(target_seqs), lengths(target_seqs))
  ),
  target = as.character(unlist(lapply(target_seqs, names))),
  sequence = as.character(unlist(target_seqs))
)


## Identify PAM sequences associated with nucleases
pam_seqs <- do.call(c, lapply(configs, function(x){
  toupper(unlist(lapply(x$Nuclease_Profiles, "[[", "PAM")))
}))

pam_grps <- stringr::str_extract(
  string = names(pam_seqs), 
  pattern = "[\\w\\-\\_]+"
)

names(pam_seqs) <- sub("[\\w\\-\\_]+.", "", names(pam_seqs), perl = TRUE)
pam_seqs <- split(pam_seqs, pam_grps)

pam_seqs_df <- data.frame(
  run_set = as.character(S4Vectors::Rle(names(pam_seqs), lengths(pam_seqs))),
  nuclease = as.character(unlist(lapply(pam_seqs, names))),
  PAM = as.character(unlist(pam_seqs))
)


## Combine into a single table for output
considered_target_seqs <- unique(unlist(treatment))
considered_nucleases <- unique(unlist(nuclease))

target_tbl <- combn_tbl %>%
  dplyr::left_join(target_seqs_df, by = c("run_set", "target")) %>%
  dplyr::left_join(pam_seqs_df, by = c("run_set", "nuclease")) %>%
  dplyr::filter(
    target %in% considered_target_seqs & nuclease %in% considered_nucleases
    )

### Log combination treatment table
if( !args$quiet ){
  cat("\nTarget Sequence Table:\n")
  print(target_tbl, right = FALSE, row.names = FALSE)
}


## Consolidate supplementary data ----
if( is.null(args$support) ){
  spec_overview <- treatment_df
}else{
  spec_overview <- supp_data
}

cond_overview <- spec_overview %>%
  dplyr::mutate(
    condition = vcollapse(
      d = dplyr::select(spec_overview, -run_set, -specimen), 
      sep = " - ", 
      fill = "NA"
    ),
    condition = factor(condition, levels = c(unique(c(condition, "Mock"))))
  ) %>%
  dplyr::select(specimen, condition)


## Read in experimental data and contatenate different sets ----
input_data_paths <- lapply(configs, function(x){
  name <- x$Run_Name
  file.path(
    "analysis", name, paste0("output/incorp_sites.", name ,".rds")
  )
})

input_data <- lapply(input_data_paths, function(x){
  if( file.exists(file.path(root_dir, x)) ){
    return(readRDS(file.path(root_dir, x)))
  }else if( file.exists(x) ){
    return(readRDS(x))
  }else{
    stop("\n  Cannot find edited_sites file: ", x, ".\n")
  }
})

names(input_data) <- names(configs)
data_type_names <- names(input_data[[1]])

input_data <- lapply(
  seq_along(input_data[[1]]), 
  function(i) dplyr::bind_rows(lapply(input_data, "[[", i), .id = "run.set")
)

names(input_data) <- data_type_names

input_data <- lapply(
  input_data, 
  function(x) x[x$specimen %in% spec_overview$specimen,] 
)

## Updating on-target data if needed
on_targets <- on_targets[names(on_targets) %in% considered_target_seqs]


# Beginnin analysis ----
if( !args$quiet ) cat("\nStarting analysis...\n")

## Specimen summary ----
# Summarize components and append to specimen table
tbl_algn_counts <- input_data$algnmts %>% 
  dplyr::mutate(specimen = factor(specimen, levels = specimen_levels)) %>%
  dplyr::group_by(specimen)


if( umitag_option ){
  
  tbl_algn_counts <- dplyr::summarise(
    tbl_algn_counts, 
    Reads = sum(count), UMItags = sum(umitag), Alignments = sum(contrib)
  )
  
}else{
  
  tbl_algn_counts <- dplyr::summarise(
    tbl_algn_counts, Reads = sum(count), Alignments = sum(contrib)
  )
  
}


spec_overview <- dplyr::left_join(
  spec_overview, tbl_algn_counts, by = "specimen"
) 


## Annotate incorporation data ----
input_data$matched_summary <- suppressMessages(dplyr::mutate(
  input_data$matched_summary,
  gene_id = assignGeneID( 
    seqnames = stringr::str_extract(edit.site, "[\\w]+"), 
    positions = as.numeric(stringr::str_extract(edit.site, "[\\w]+$")), 
    reference = ref_genome, 
    ref.genes = ref_genes, 
    onco.genes = onco_genes, 
    special.genes = special_genes
  )
))

input_data$paired_regions <- suppressMessages(dplyr::mutate(
  input_data$paired_regions,     
  gene_id = assignGeneID(
    seqnames = seqnames, 
    positions = mid, 
    reference = ref_genome, 
    ref.genes = ref_genes, 
    onco.genes = onco_genes, 
    special.genes = special_genes
  )
))

input_data$pile_up_summary <- suppressMessages(dplyr::mutate(
  input_data$pile_up_summary,
  gene_id = assignGeneID( 
    seqnames = stringr::str_extract(clus.ori, "[\\w]+"), 
    positions = as.numeric(stringr::str_extract(clus.ori, "[\\w]+$")), 
    reference = ref_genome, 
    ref.genes = ref_genes, 
    onco.genes = onco_genes, 
    special.genes = special_genes
  )
))


## On-target summary ----
# Algnmts
tbl_ot_algn <- input_data$algnmts %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(
    ot_algns = pNums(
      sum(contrib * as.integer(edit.site %in% expandPosStr(on_targets)))
    ),
    ot_algns_pct = 100 * sum(
        contrib * as.integer(edit.site %in% expandPosStr(on_targets))
      ) /
      sum(contrib)
  ) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# Probable edited sites
tbl_ot_prob <- input_data$probable_algns %>% 
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(
    ot_prob = pNums(
      sum(contrib * as.integer(edit.site %in% expandPosStr(on_targets)))
    ),
    ot_prob_pct = 100 * sum(
        contrib * as.integer(edit.site %in% expandPosStr(on_targets))
      ) /
      sum(contrib)
  ) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# Pile ups of read alignments
tbl_ot_pile <- input_data$pile_up_algns %>% 
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(
    ot_pile = pNums(
      sum(contrib * as.integer(edit.site %in% expandPosStr(on_targets)))
    ),
    ot_pile_pct = 100 * sum(
        contrib * as.integer(edit.site %in% expandPosStr(on_targets))
      ) /
      sum(contrib)
  ) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# Paired or flanking algnments
tbl_ot_pair <- input_data$paired_regions %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen))),
    on.off.target = factor(
      on.off.target, levels = c("On-target", "Off-target")
    )
  ) %>%
  dplyr::group_by(specimen, on.off.target) %>%
  dplyr::summarise(cnt = sum(algns)) %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(
    ot_pair = pNums(sum(ifelse(on.off.target == "On-target", cnt, 0))),
    ot_pair_pct = 100 * sum(ifelse(on.off.target == "On-target", cnt, 0)) /
      sum(cnt)
  ) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# Guide RNA matched within 6 mismatches
tbl_ot_match <- input_data$matched_summary %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::group_by(specimen, on.off.target) %>%
  dplyr::summarise(cnt = sum(algns)) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(specimen) %>%
  dplyr::summarise(
    ot_match = pNums(sum(ifelse(on.off.target == "On-target", cnt, 0))),
    ot_match_pct = 100 * sum(ifelse(on.off.target == "On-target", cnt, 0)) /
      sum(cnt)
  ) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# Summary table
ot_tbl_summary <- dplyr::left_join(
    treatment_df, cond_overview, by = "specimen"
  ) %>%
  dplyr::mutate(
    condition = factor(
      ifelse(is.na(condition), "Mock", paste(condition)), 
      levels = levels(condition)
    )
  )

ot_tbl_summary <- Reduce(
    function(x,y) dplyr::left_join(x, y, by = "specimen"),
    list(
      tbl_ot_algn[,c(1,3)], tbl_ot_pile[,c(1,3)],
      tbl_ot_pair[,c(1,3)], tbl_ot_match[,c(1,3)]
    ),
    init = ot_tbl_summary
  ) %>%
  dplyr::mutate(specimen = factor(specimen, levels = specimen_levels)) %>%
  dplyr::arrange(specimen) %>%
  dplyr::select(-treatment)



## On-target incorporation distribution ----
on_tar_dists <- input_data$matched_algns %>%
  dplyr::filter(on.off.target == "On-target") %>%
  dplyr::mutate(
    target = stringr::str_extract(string = target.match, pattern = "[\\w]+"),
    pos = as.numeric(
      stringr::str_extract(string = edit.site, pattern = "[0-9]+$")
    ),
    edit.site.dist = ifelse(strand == "+", start - pos, end - pos),
    specimen = factor(specimen, levels = levels(cond_overview$specimen))
  ) %>%
  dplyr::left_join(cond_overview, by = "specimen") %>%
  dplyr::select(
    run.set, specimen, target, condition, 
    edit.site, edit.site.dist, strand, contrib)

on_tar_dens <- lapply(
  split(on_tar_dists, on_tar_dists$condition), 
  function(x){
    
    if( nrow(x) >= 10 ){
      return(
        density(abs(x$edit.site.dist), from = 0, to = upstream_dist, bw = 1)
      )
    }else{
      return(NA)
    }
    
  }
)

on_tar_dists <- dplyr::group_by(
    on_tar_dists, 
    condition, target, edit.site.dist, strand
  ) %>%
  dplyr::summarise(cnt = sum(contrib)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    strand.cnt = ifelse(
      strand == "+", log(cnt, base = 10), -log(cnt, base = 10))
  )

if( length(unique(cond_overview$condition)) == 1 ){
  on_tar_dists$condition <- " "
}

if( is.null(args$support) ){
  sites_included <- on_tar_dists %>% dplyr::group_by(target)
}else{
  sites_included <- on_tar_dists %>% dplyr::group_by(condition, target)
}

sites_included <- dplyr::summarise(
    sites_included,
    prop = 100 * sum(cnt[
      abs(edit.site.dist) <= upstream_dist & 
        abs(edit.site.dist) >= -downstream_dist
    ]) / sum(cnt),
    x_pos = upstream_dist,
    y_pos = 0.8 * min(strand.cnt[
      abs(edit.site.dist) <= upstream_dist & 
        abs(edit.site.dist) >= -downstream_dist
    ])
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(prop = paste0(pNums(prop, digits = 4), "%"))


## Off-target summary ----
# All alignments
tbl_ft_algn <- input_data$algnmts %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::filter(!edit.site %in% expandPosStr(on_targets)) %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(ft_algns = dplyr::n_distinct(clus.ori)) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# Probable edit sites
tbl_ft_prob <- input_data$probable_algns %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::filter(on.off.target == "Off-target") %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(ft_prob = dplyr::n_distinct(clus.ori)) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# Pile ups
tbl_ft_pile <- input_data$pile_up_algns %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::filter(on.off.target == "Off-target") %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(ft_pile = dplyr::n_distinct(clus.ori)) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# Paired or flanked loci
tbl_ft_pair <- input_data$paired_regions %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::filter(on.off.target == "Off-target") %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(ft_pair = n()) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# target sequence matched
tbl_ft_match <- input_data$matched_summary %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::filter(on.off.target == "Off-target") %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(ft_match = n()) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# Summary table
ft_tbl_summary <- dplyr::left_join(
    treatment_df, cond_overview, by = "specimen"
  ) %>%
  dplyr::mutate(
    condition = factor(
      ifelse(is.na(condition), "Mock", paste(condition)), 
      levels = levels(condition)
    )
  )

ft_tbl_summary <- Reduce(
    function(x,y) dplyr::left_join(x, y, by = "specimen"),
    list(tbl_ft_algn, tbl_ft_pile, tbl_ft_pair, tbl_ft_match),
    init = ft_tbl_summary
  ) %>%
  dplyr::mutate(specimen = factor(specimen, levels = specimen_levels)) %>%
  dplyr::arrange(specimen) %>%
  dplyr::select(-treatment)


## Onco-gene enrichment analysis ----
rand_sites <- selectRandomSites(
  num = nrow(input_data$paired_regions) + nrow(input_data$matched_summary), 
  ref.genome = ref_genome, 
  drop.extra.seqs = TRUE, 
  rnd.seed = 1
)

rand_sites$gene_id <- suppressMessages(assignGeneID(
  seqnames = seqnames(rand_sites), 
  positions = start(rand_sites), 
  reference = ref_genome, 
  ref.genes = ref_genes, 
  onco.genes = onco_genes, 
  special.genes = special_genes
))

rand_df <- data.frame(
  condition = "Random", 
  "total" = length(rand_sites), 
  "onco" = sum(stringr::str_detect(rand_sites$gene_id, "~")), 
  "special" = sum(stringr::str_detect(rand_sites$gene_id, "!"))
)

paired_list <- split(
  x = input_data$paired_regions, 
  f = cond_overview$condition[
    match(input_data$paired_regions$specimen, cond_overview$specimen)
  ]
)

paired_df <- dplyr::bind_rows(
  lapply(
    paired_list, 
    function(df){
    
    data.frame(
      "total" = nrow(df), 
      "onco" = sum(stringr::str_detect(df$gene_id, "~")), 
      "special" = sum(stringr::str_detect(df$gene_id, "!"))
    )
    
    }
  ), 
  .id = "condition"
)

matched_list <- split(
  x = input_data$matched_summary, 
  f = cond_overview$condition[
    match(input_data$matched_summary$specimen, cond_overview$specimen)
])

matched_df <- dplyr::bind_rows(
  lapply(
    matched_list, 
    function(df){
      
      data.frame(
        "total" = nrow(df), 
        "onco" = sum(stringr::str_detect(df$gene_id, "~")), 
        "special" = sum(stringr::str_detect(df$gene_id, "!"))
      )
      
    }
  ), 
  .id = "condition"
)

enrich_df <- dplyr::bind_rows(
    list(
      "Reference" = rand_df, 
      "Flanking Pairs" = paired_df, 
      "Target Matched" = matched_df), 
    .id = "origin"
  ) %>%
  dplyr::filter(total > 0)

enrich_df$onco.p.value <- p.adjust(
  sapply(
    seq_len(nrow(enrich_df)), 
    function(i){
      
      ref <- enrich_df[1, c("total", "onco")]
      query <- enrich_df[i, c("total", "onco")]
      ref$diff <- abs(diff(as.numeric(ref)))
      query$diff <- abs(diff(as.numeric(query)))
      
      fisher.test(as.matrix(rbind(
        ref[,c("diff", "onco")], query[,c("diff", "onco")]
      )))$p.value
      
    }
  ), 
  method = "BH"
)

enrich_df$special.p.value <- p.adjust(
  sapply(
    seq_len(nrow(enrich_df)), 
    function(i){
      
      ref <- enrich_df[1, c("total", "special")]
      query <- enrich_df[i, c("total", "special")]
      ref$diff <- abs(diff(as.numeric(ref)))
      query$diff <- abs(diff(as.numeric(query)))
      
      fisher.test(as.matrix(rbind(
        ref[,c("diff", "special")], query[,c("diff", "special")]
      )))$p.value
      
    }
  ), 
  method = "BH"
)


## Off-target sequence analysis ----
ft_MESL <- input_data$matched_algns %>%
  dplyr::mutate(
    edit.site.dist = abs(ifelse(
      strand == "+", 
      start - as.numeric(stringr::str_extract(edit.site, "[0-9]+$")), 
      as.numeric(stringr::str_extract(edit.site, "[0-9]+$")) - end
    )),
    specimen = factor(specimen, levels = levels(cond_overview$specimen))
  ) %>%
  dplyr::left_join(cond_overview, by = "specimen") %>%
  dplyr::mutate(order = seq_len(n())) 

if( nrow(ft_MESL) > 0 ){
  
  ft_MESL <- ft_MESL %>%
    dplyr::group_by(order) %>%
    dplyr::mutate(
      ESL = predictESProb(
        x = edit.site.dist, density = on_tar_dens[[condition]]
      ),
      gene_id = input_data$matched_summary$gene_id[
        match(edit.site, input_data$matched_summary$edit.site)
      ]
    ) %>%
    dplyr::group_by(condition, edit.site, gene_id) %>%
    dplyr::summarise(MESL = 100 * max(c(0,ESL), na.rm = TRUE)) %>%
    dplyr::ungroup()
  
}else{
  
  ft_MESL <- ft_MESL %>%
    dplyr::mutate(
      MESL = vector(mode = "numeric"), gene_id = vector(mode = "character")) %>%
    dplyr::select(condition, edit.site, gene_id, MESL)
  
}

ft_seqs <- input_data$matched_summary %>%
  dplyr::select(
    specimen, aligned.sequence, target.match, edit.site,
    target.mismatch, on.off.target, algns, gene_id
  ) %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = levels(cond_overview$specimen))
  ) %>%
  dplyr::left_join(cond_overview, by = "specimen") %>%
  dplyr::left_join(ft_MESL, by = c("condition", "edit.site", "gene_id"))

if( is.null(args$support) ){
  
  ft_seqs <- dplyr::group_by(
      ft_seqs, 
      target.match, edit.site, aligned.sequence, 
      target.mismatch, on.off.target, gene_id
    ) %>%
    dplyr::summarise(algns = sum(algns), MESL = max(MESL, na.rm = TRUE))
  
}else{
  
  ft_seqs <- dplyr::group_by(
      ft_seqs,
      condition, target.match, edit.site, aligned.sequence, 
      target.mismatch, on.off.target, gene_id
    ) %>%
    dplyr::summarise(algns = sum(algns), MESL = max(MESL, na.rm = TRUE))
  
}

ft_seqs <- dplyr::arrange(
    ft_seqs, desc(algns), desc(MESL), target.mismatch
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    on.off.target = stringr::str_extract(on.off.target, "[\\w]+")
  ) %>%
  dplyr::rename(
    "target" = on.off.target, 
    "mismatch" = target.mismatch, 
    "target.seq" = target.match,
    "aligns" = algns
  )

if( is.null(args$support) ){
  
  ft_seqs_list <- split(ft_seqs, ft_seqs$target.seq)
  
}else{
  
  ft_seqs_conds <- dplyr::arrange(ft_seqs, condition, target.seq) %$% 
    unique(paste0(condition, " - ", target.seq))
  
  ft_seqs_list <- split(
    x = ft_seqs, 
    f = factor(
      paste0(ft_seqs$condition, " - ", ft_seqs$target.seq), 
      levels = ft_seqs_conds
    )
  )
  
}

if( !args$quiet ) cat("Analysis complete.\nStarting report generation...\n")

# Data consolidated for output object ----
set_names <- ifelse(
  length(configs) ==  1, 
  names(configs), 
  paste0(
    paste(names(configs)[seq_len(length(configs)-1)], collapse = ", "),
    ", and ", 
    names(configs)[length(configs)]
  )
)

## Write output file
saveRDS(
  object = list(
    "params" = list(
      "set_names" = set_names, 
      "configs" = configs, 
      "soft_version" = soft_version, 
      "build_version" = build_version,
      "specimen_levels" = specimen_levels
    ),
    "spec_info" = list(
      "sample_info" = sample_info, 
      "target_seqs" = target_seqs,
      "target_tbl" = target_tbl,
      "on_targets" = on_targets, 
      "treatment" = treatment, 
      "treatment_df" = treatment_df,
      "nuclease" = nuclease,
      "nuclease_df" = nuclease_df,
      "nuclease_profiles" = nuc_profiles,
      "supp_data" = supp_data, 
      "spec_overview" = spec_overview, 
      "cond_overview" = cond_overview
    ),
    "incorp_data" = input_data, 
    "summary_tbls" = list(
      "ot_tbl_summary" = ot_tbl_summary, 
      "ft_tbl_summary" = ft_tbl_summary
    ), 
    "edit_models" = list(
      "on_tar_dists" = on_tar_dists, 
      "on_tar_dens" = on_tar_dens, 
      "sites_included" = sites_included
    ),
    "enrich_data" = list(
      "rand_sites" = rand_sites, 
      "rand_df" = rand_df, 
      "enrich_df" = enrich_df
    ),
    "ft_data" = ft_seqs_list
  ),
  file = args$output
)
  
if( !file.exists(args$output) ){
  stop("\n  Cannot verify existence of output file: ", args$output, "\n")
}else{
  if( !args$quiet ){
    cat("Evaluation complete, output writen to:", args$output, "\n")
  }
  q(status = 0)
}

