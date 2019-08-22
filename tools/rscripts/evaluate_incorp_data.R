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
  "--stat", nargs = 1, type = "character", default = FALSE, 
  help = paste(
    "File name to be written in output directory of read couts for each",
    "sample. CSV file format. ie. test.stat.csv."
  )
)

parser$add_argument(
  "--override", action = "store_true", 
  help = "Override software and build version control checks."
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

add_packs <- c("magrittr", "knitr", "iguideSupport")

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

submat <- banmat()

## Determine processing parameters
## Some parameters will need to be an "all or nothing" approach, including:
##   - UMItags
##   - recoverMultihits
## Depending on these parameters others (upstream/downstream_dist, ...) may need
## to be consistent between runs otherwise, the primary config file (first one),
## will be used for parameterization.

umitag_option <- all(unlist(lapply(configs, "[[", "UMItags")))
multihit_option <- all(unlist(lapply(configs, "[[", "recoverMultihits")))

if( multihit_option ){
  
  upstream_dist <- unique(sapply(configs, function(x) x$upstreamDist))
  downstream_dist <- unique(sapply(configs, function(x) x$downstreamDist))
  pile_up_min <- unique(sapply(configs, function(x) x$pileUpMin))

  if( 
    length(upstream_dist) > 1 | 
    length(downstream_dist) > 1 | 
    length(pile_up_min) > 1 
  ){
    
    stop(
      "\n  Inconsistant upstream or downstream distances between config files.",
      "\n  Comparisons between groups with different run specific criteria", 
      "\n  is not recommended when considering the recover multihit option.\n"
    )
    
  }
  
}else{
  
  upstream_dist <- configs[[1]]$upstreamDist
  downstream_dist <- configs[[1]]$downstreamDist
  pile_up_min <- configs[[1]]$pileUpMin
  
}

max_target_mismatch <- configs[[1]]$maxTargetMismatch


## Combine sampleInfo files

sample_info <- lapply(configs, loadSampleInfo, root_dir) %>%
  dplyr::bind_rows(.id = "run_set")

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
  
  supp_data <- data.table::fread(support_path, data.table = FALSE)
  
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
on_targets <- unlist(lapply(configs, "[[", "On_Target_Sites")) %>%
  data.frame(id = names(.), target = ., row.names = NULL) %>%
  dplyr::mutate(
    id = stringr::str_replace(
      string = id, pattern = stringr::fixed("."), replacement = ":"
    ),
    id = stringr::str_extract(string = id, pattern = "[\\w\\_\\-\\'\\.]+$"),
    id = stringr::str_extract(string = id, pattern = "[\\w\\_\\-\\.]+")
  ) %>%
  dplyr::distinct() %$%
  structure(target, names = id)


## Treatment across runs
treatment_df <- lapply(configs, getTreatmentInfo, root_dir) %>%
  dplyr::bind_rows(.id = "run_set") %>%
  dplyr::filter(specimen %in% specimen_levels) %>%
  dplyr::mutate(specimen = factor(specimen, levels = specimen_levels)) %>%
  dplyr::arrange(specimen)

### If only one condition present besides 'Mock', then switch 'Mock' treated
### specimens to condition for analysis background
### This block of code may not be appropriate for all circumstances.
if( any(tolower(treatment_df$treatment) == "mock") ){
  
  all_treatments <- unique(unlist(strsplit(treatment_df$treatment, ";")))
  
  all_treatments <- all_treatments[
    !grepl("mock", tolower(all_treatments))
  ]
  
  treatment_df$treatment <- ifelse(
    tolower(treatment_df$treatment) == "mock", 
    paste(all_treatments, collapse = ";"), 
    treatment_df$treatment
  )
  
}


treatment <- structure(
  strsplit(treatment_df$treatment, ";"), 
  names = as.character(treatment_df$specimen)
)
  

## Nucleases used across runs
nuc_profiles <- unlist(
  unname(lapply(configs, "[[", "Nuclease_Profiles")), 
  recursive = FALSE
)

nuc_profiles <- nuc_profiles[
  match(unique(names(nuc_profiles)), names(nuc_profiles))
]

nuclease_df <- lapply(configs, getNucleaseInfo, root_dir) %>%
  dplyr::bind_rows(.id = "run_set") %>%
  dplyr::filter(specimen %in% specimen_levels) %>%
  dplyr::mutate(specimen = factor(specimen, levels = specimen_levels)) %>%
  dplyr::arrange(specimen)

### If only one condition present besides 'Mock', then switch 'Mock' treated
### specimens to condition for analysis background
### This block of code may not be appropriate for all circumstances.
if( any(tolower(nuclease_df$nuclease) == "mock") ){
  
  all_nucleases <- unique(unlist(strsplit(nuclease_df$nuclease, ";")))
  
  all_nucleases <- all_nucleases[
    !grepl("mock", tolower(all_nuclease))
  ]
  
  nuclease_df$nuclease <- ifelse(
    tolower(nuclease_df$nuclease) == "mock", 
    paste(all_nucleases, collapse = ";"), 
    nuclease_df$nuclease
  )
  
}

nuclease <- structure(
  strsplit(nuclease_df$nuclease, ";"), 
  names = as.character(nuclease_df$specimen)
)


## Combine treatment and nuclease profiles
nuclease_treaments <- dplyr::left_join(
  nuclease_df, treatment_df, by = c("run_set", "specimen")
)

target_combn <- structure(
  strsplit(nuclease_treaments$treatment, ";"), 
  names = as.character(nuclease_treaments$specimen)
)

combn_tbl <- data.frame(
    run_set = nuclease_treaments$run_set[
      as.vector(S4Vectors::match(
        S4Vectors::Rle(names(target_combn), lengths(target_combn)), 
        nuclease_treaments$specimen
      ))
    ],
    nuclease = nuclease_treaments$nuclease[
      as.vector(S4Vectors::match(
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

on_targets <- on_targets[names(on_targets) %in% considered_target_seqs]

target_tbl <- combn_tbl %>%
  dplyr::left_join(target_seqs_df, by = c("run_set", "target")) %>%
  dplyr::left_join(pam_seqs_df, by = c("run_set", "nuclease")) %>%
  dplyr::filter(
    target %in% considered_target_seqs & nuclease %in% considered_nucleases
  )

uniq_target_df <- target_tbl %>%
  dplyr::distinct(target, sequence, PAM)

uniq_target_seqs <- Biostrings::DNAStringSet(
  structure(uniq_target_df$sequence, names = uniq_target_df$target),
  use.names = TRUE
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
  spec_overview <- supp_data %>%
    dplyr::mutate(run_set = "supp_data")
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


# Beginning analysis ----
if( !args$quiet ) cat("\nStarting analysis...\n")

## Read in experimental data and contatenate different sets
input_data <- lapply(configs, function(x){
    name <- x$Run_Name
    
    path <- file.path(
      "analysis", name, paste0("output/incorp_sites.", name ,".rds")
    )
    
    if( file.exists(file.path(root_dir, path)) ){
      y <- readRDS(file.path(root_dir, path))
    }else if( file.exists(path) ){
      y <- readRDS(path)
    }else{
      stop("\n  Cannot find incorp_sites file: ", x, ".\n")
    }

    y$reads %>%
      dplyr::mutate(
        soft.version = y$soft_version,
        build.version = y$build_version
      )
    
  }) %>%
  dplyr::bind_rows(.id = "run.set") %>%
  dplyr::mutate(
    specimen = stringr::str_extract(sampleName, pattern = "[\\w]+")
  ) %>%
  dplyr::filter(specimen %in% spec_overview$specimen)

if( !multihit_option ){
  input_data <- dplyr::filter(input_data, type == "uniq")
}

## Check versioning for imported data ----
vc_check <- input_data %>%
  dplyr::distinct(run.set, soft.version, build.version)

input_data <- dplyr::select(input_data, -soft.version, -build.version)

cat("\nVersioning:\n")
print(vc_check, right = FALSE, row.names = FALSE)

if( dplyr::n_distinct(vc_check$soft.version) > 1 | 
      dplyr::n_distinct(vc_check$build.version) > 1 ){

  if( args$override ){
    warning("Data processed under different software versions.")
  }else{
    stop("\n  Data processed with inconsistent software versions.")
  }

}

## Format input alignments ----
## Determine abundance metrics, with or without UMItags
algnmts_summaries <- list(
  count = dplyr::quo(sum(contrib)),
  umitag = if( umitag_option ){
    dplyr::quo(sum(as.integer(!duplicated(umitag[!is.na(umitag)])) * contrib))
  },
  contrib = dplyr::quo(max(contrib))
)

algnmts_summaries <- algnmts_summaries[!sapply(algnmts_summaries, is.null)]

algnmts <- input_data %>%
  dplyr::arrange(desc(contrib)) %>%
  dplyr::group_by(seqnames, start, end, strand, specimen, sampleName) %>%
  dplyr::summarise(!!! algnmts_summaries) %>%
  dplyr::ungroup() %>%
  as.data.frame()

## Generate a sample table of the data for log purposes
sample_index <- ifelse(nrow(algnmts) > 10, 10, nrow(algnmts))
sample_index <- sample(seq_len(nrow(algnmts)), sample_index, replace = FALSE)

cat("\nSample of aligned templates:\n")

print(
  data.frame(algnmts[sample_index,]),
  right = FALSE,
  row.names = FALSE
)

cat(paste0("\nNumber of alignments: ", nrow(algnmts), "\n"))

rm(sample_index)

## Transform the data into a GRanges object
algnmts_gr <- GenomicRanges::GRanges(
  seqnames = algnmts$seqnames,
  ranges = IRanges::IRanges(start = algnmts$start, end = algnmts$end),
  strand = algnmts$strand,
  seqinfo = GenomeInfoDb::seqinfo(ref_genome)
)

GenomicRanges::mcols(algnmts_gr) <- dplyr::select(
  algnmts, c(specimen, sampleName, count, if( umitag_option ) umitag, contrib)
)

# Analyze alignments ----
## Identify groups of alignments or pileups of aligned fragments
## These pileups give strong experimental evidence of directed incorporation of
## the dsODN into a region. Initially, pileups are identified and then checked 
## for pairing, or if there is another pileup on the opposite strand in close 
## proximity.
algnmts_gr$clus.ori <- pileupCluster(
  gr = algnmts_gr, 
  grouping = "specimen", 
  maxgap = 0L, 
  return = "simple"
)

algnmts_gr$paired.algn <- identifyPairedAlgnmts(
  gr = algnmts_gr, 
  grouping = "specimen", 
  maxgap = upstream_dist * 2
)

algnmts_grl <- split(algnmts_gr, unlist(nuclease)[algnmts_gr$specimen])

annot_clust_info <- dplyr::bind_rows(lapply(
  seq_along(algnmts_grl), 
  function(i, grl){

    gr <- grl[[i]]
    nuc <- names(grl)[i]
    
    if( !nuc %in% names(nuc_profiles) ){
      nuc_profile <- NULL
    }else{
      nuc_profile <- nuc_profiles[[nuc]]
    }
    
    ## Create a GRange with only the unique cluster origins
    split_clus_id <- stringr::str_split(
      string = unique(paste0(gr$specimen, ":", gr$clus.ori)), 
      pattern = ":", 
      simplify = TRUE
    )
    
    algn_clusters <- GenomicRanges::GRanges(
      seqnames = split_clus_id[,2],
      ranges = IRanges::IRanges(
        start = as.numeric(split_clus_id[,4]), width = 1
      ),
      strand = split_clus_id[,3],
      seqinfo = GenomeInfoDb::seqinfo(ref_genome)
    )
    
    algn_clusters$specimen <- split_clus_id[,1]
    algn_clusters$clus.ori <- vcollapse(split_clus_id[, 2:4], sep = ":")
    
    algn_clusters$clus.seq <- getSiteSeqs(
      gr = algn_clusters, 
      upstream.flank = upstream_dist, 
      downstream.flank = downstream_dist, 
      ref.genome = ref_genome
    )
    
    ## Identify which target sequences binding near clusters
    if( !is.null(nuc_profile) ){
      
      algn_clusters <- compareTargetSeqs(
        gr.with.sequences = algn_clusters, 
        seq.col = "clus.seq", 
        target.seqs = uniq_target_seqs,
        tolerance = max_target_mismatch,
        nuc.profile = nuc_profile,
        submat = submat, 
        upstream.flank = upstream_dist, 
        downstream.flank = downstream_dist
      )
      
    }else{
      
      algn_clusters$target.match <- "No_valid_match"
      algn_clusters$target.mismatch <- NA
      algn_clusters$target.score <- NA
      algn_clusters$aligned.sequence <- NA
      algn_clusters$edit.site <- NA
      
    }
    
    as.data.frame(GenomicRanges::mcols(algn_clusters))
    
  },
  grl = algnmts_grl
))


## Merge the target sequence alignment information from the clusters back to all
## unique alignments
algnmts <- as.data.frame(merge(
  x = as.data.frame(algnmts_gr), 
  y = dplyr::select(annot_clust_info, -clus.seq),
  by = c("specimen", "clus.ori")
))

## Change guideRNA.match to No_Valid_Match if an inappropriate gRNA is annotated
algnmts$target.match <- filterInappropriateComparisons(
  guideRNA.match = algnmts$target.match, 
  specimen = algnmts$specimen, 
  treatment = treatment
)

## Fragment pileups, paired clustering, and guideRNA alignments have been used 
## to characterize the incorporation sites analyzed here. Each metric will be 
## used to create a list of incorporation sites that may be nuclease cut sites. 
## The following identifies which alignments are associated with each of these 
## criteria.
tbl_clus_ori <- algnmts %>% 
  dplyr::group_by(specimen, clus.ori) %>%
  dplyr::filter(n() >= pile_up_min) %>%
  dplyr::ungroup() %$%
  table(clus.ori)

idx_clus_ori <- which(algnmts$clus.ori %in% names(tbl_clus_ori))

tbl_paired_algn <- algnmts %>%
  dplyr::filter(!is.na(paired.algn)) %$%
  table(paired.algn)

idx_paired_algn <- which(algnmts$paired.algn %in% names(tbl_paired_algn))

idx_matched <- which(algnmts$target.match != "No_valid_match")

idx_combined <- sort(unique(c(idx_clus_ori, idx_paired_algn, idx_matched)))

idx_df <- data.frame(
  "Type" = c("PileUp", "Paired", "gRNA_Matched", "Combined"),
  "Counts" = sapply(
    list(idx_clus_ori, idx_paired_algn, idx_matched, idx_combined), 
    length
  )
)

cat("\nTable of uniquely aligned template counts:\n")
print(idx_df, right = FALSE, row.names = FALSE) 
cat(paste0("\nTotal number of alignments: ", nrow(algnmts), "\n"))

probable_algns <- algnmts[idx_combined,]

probable_algns$on.off.target <- ifelse(
  probable_algns$edit.site %in% expandPosStr(on_targets), 
  "On-target", 
  "Off-target"
)

cat("\nOn / Off target alignment counts:\n")
print(table(probable_algns$on.off.target))


## Create summary and output formated object related to each of the criteria for
## edited site detection.

## Matched alignments
matched_algns <- probable_algns[
  probable_algns$target.match != "No_valid_match",
  ]

matched_summaries <- list(
  on.off.target = dplyr::quo(
    paste(sort(unique(on.off.target)), collapse = ";")
  ),
  paired.algn = dplyr::quo(paste(sort(unique(paired.algn)), collapse = ";")),
  count = dplyr::quo(sum(count)), 
  umitag = if( umitag_option ) dplyr::quo(sum(umitag)),
  algns = dplyr::quo(sum(contrib)),
  orient = dplyr::quo(paste(sort(unique(as.character(strand))), collapse = ";"))
)

matched_summaries <- matched_summaries[!sapply(matched_summaries, is.null)]

matched_summary <- matched_algns %>%
  dplyr::mutate(
    target.match = stringr::str_replace(
      string = target.match, 
      pattern = "\\:\\([\\w]+\\)$",
      replacement = ""
    )
  ) %>%
  dplyr::group_by(
    specimen, edit.site, aligned.sequence, target.match, target.mismatch
  ) %>%
  dplyr::summarise(!!! matched_summaries) %>%
  dplyr::ungroup() %>% 
  dplyr::arrange(specimen, target.match, desc(algns)) %>%
  as.data.frame()

## Paired alignments
paired_algns <- probable_algns[
  probable_algns$paired.algn %in% names(tbl_paired_algn),
]

paired_summaries <- list(
  seqnames = dplyr::quo(unique(seqnames)),
  start = dplyr::quo(min(pos)), 
  end = dplyr::quo(max(pos)), 
  mid = dplyr::quo(start + (end-start)/2),
  strand = dplyr::quo("*"), 
  width = dplyr::quo(end - start), 
  count = dplyr::quo(sum(count)), 
  umitag = if( umitag_option ) dplyr::quo(sum(umitag)), 
  algns = dplyr::quo(sum(contrib))
)

paired_summaries <- paired_summaries[!sapply(paired_summaries, is.null)]

paired_regions <- paired_algns %>%
  dplyr::group_by(specimen, paired.algn, strand) %>%
  dplyr::mutate(pos = ifelse(strand == "+", min(start), max(end))) %>%
  dplyr::group_by(specimen, paired.algn) %>%
  dplyr::summarise(!!! paired_summaries) %>%
  dplyr::ungroup()
  

if( nrow(paired_regions) > 0 ){
  
  paired_regions <- paired_regions %>%
    dplyr::group_by(specimen, paired.algn) %>%
    dplyr::mutate(
      on.off.target = ifelse(
        any(sapply(
          expandPosStr(unlist(on_targets[
            which(
              stringr::str_extract(
                names(on_targets), "[\\w\\-\\_\\.]+") %in% 
                treatment[[specimen]]
            )
            ])),
          function(x, seq, st, en){
            
            match_seq <- seq == stringr::str_extract(x, "[\\w]+")
            
            within_start <- st <= 
              as.numeric(stringr::str_extract(x, "[\\w]+$")) + downstream_dist
            
            within_end <- en >= 
              as.numeric(stringr::str_extract(x, "[\\w]+$")) - downstream_dist
            
            match_seq & within_start & within_end
            
          }, 
          seq = seqnames, 
          st = start, 
          en = end
        )), 
        "On-target", 
        "Off-target"
      )
    ) %>%
    dplyr::ungroup() %>% 
    as.data.frame()
  
}else{
  
  paired_regions <- dplyr::mutate(
    paired_regions,
    on.off.target = vector(mode = "character")
  )
  
}

## Pile up alignments
pile_up_algns <- probable_algns[
  probable_algns$clus.ori %in% names(tbl_clus_ori),
]

pile_up_summaries <- list(
  on.off.target = dplyr::quo(
    paste(sort(unique(on.off.target)), collapse = ";")
  ),
  paired.algn = dplyr::quo(paste(sort(unique(paired.algn)), collapse = ";")),
  count = dplyr::quo(sum(count)), 
  umitag = if( umitag_option ) dplyr::quo(sum(umitag)),
  algns = dplyr::quo(sum(contrib))
)

pile_up_summaries <- pile_up_summaries[!sapply(pile_up_summaries, is.null)]

pile_up_summary <- pile_up_algns %>%
  dplyr::mutate(
    target.match = stringr::str_replace(
      string = target.match, 
      pattern = "\\:\\([\\w]+\\)$",
      replacement = ""
    )
  ) %>%
  dplyr::group_by(specimen, clus.ori) %>%
  dplyr::summarise(!!! pile_up_summaries) %>%
  dplyr::ungroup() %>% 
  dplyr::arrange(specimen, desc(algns)) %>%
  as.data.frame()


# Generate stats if requested ----
## If requested, generate stats from the analysis for qc.

if( args$stat != FALSE ){
  
  stat_summary <- function(x, y){
    
    x %>%
      dplyr::mutate(metric = y) %>%
      dplyr::group_by(sampleName, metric) %>%
      dplyr::summarize(count = sum(contrib)) %>%
      dplyr::ungroup()
    
  }
  
  total_stat <- stat_summary(algnmts, "total.algns")
  combined_stat <- stat_summary(probable_algns, "combined.algns")
  pileup_stat <- stat_summary(pile_up_algns, "pileup.algns")
  paired_stat <- stat_summary(paired_algns, "paired.algns")
  matched_stat <- stat_summary(matched_algns, "matched.algns")
  
  on_tar_stat <- dplyr::filter(
    matched_algns, on.off.target == "On-target"
  ) %>%
    stat_summary("ontarget.algns")
  
  off_tar_stat <- dplyr::filter(
    matched_algns, on.off.target == "Off-target"
  ) %>%
    stat_summary("offtarget.algns")
  
  stat <- dplyr::bind_rows(
    total_stat, combined_stat, pileup_stat, paired_stat, 
    matched_stat, on_tar_stat, off_tar_stat)
  
  write.table(
    x = stat, file = args$stat, 
    sep = ",", row.names = FALSE, 
    col.names = FALSE, quote = FALSE
  )
  
}


## Specimen summary ----
# Summarize components and append to specimen table
tbl_algn_summaries <- list(
  Reads = dplyr::quo(sum(count)), 
  UMItags = if( umitag_option ) dplyr::quo(sum(umitag)), 
  Alignments = dplyr::quo(sum(contrib))
)

tbl_algn_summaries <- tbl_algn_summaries[!sapply(tbl_algn_summaries, is.null)]

tbl_algn_counts <- algnmts %>% 
  dplyr::mutate(specimen = factor(specimen, levels = specimen_levels)) %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(!!! tbl_algn_summaries)

spec_overview <- dplyr::left_join(
  spec_overview, tbl_algn_counts, by = "specimen"
) 


## Annotate incorporation data ----
matched_summary <- suppressMessages(dplyr::mutate(
  matched_summary,
  gene_id = assignGeneID( 
    seqnames = stringr::str_extract(edit.site, "[\\w]+"), 
    positions = as.numeric(stringr::str_extract(edit.site, "[\\w]+$")), 
    reference = ref_genome, 
    ref.genes = ref_genes, 
    onco.genes = onco_genes, 
    special.genes = special_genes
  )
))

paired_regions <- suppressMessages(dplyr::mutate(
  paired_regions,     
  gene_id = assignGeneID(
    seqnames = seqnames, 
    positions = mid, 
    reference = ref_genome, 
    ref.genes = ref_genes, 
    onco.genes = onco_genes, 
    special.genes = special_genes
  )
))

pile_up_summary <- suppressMessages(dplyr::mutate(
  pile_up_summary,
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
tbl_ot_algn <- algnmts %>%
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
tbl_ot_prob <- probable_algns %>% 
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
tbl_ot_pile <- pile_up_algns %>% 
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
tbl_ot_pair <- paired_regions %>%
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
tbl_ot_match <- matched_summary %>%
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

tbl_ot_eff <- matched_summary %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::group_by(specimen, on.off.target, target.match) %>%
  dplyr::summarise(cnt = sum(algns)) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(specimen, target.match) %>%
  dplyr::summarise(
    ot_eff_pct = 100 * sum(ifelse(on.off.target == "On-target", cnt, 0)) /
      sum(cnt)
  ) %>%
  dplyr::ungroup() %>%
  tidyr::spread(key = target.match, value = ot_eff_pct) %>%
  tidyr::complete(specimen) %>%
  as.data.frame() %>%
  dplyr::left_join(
    cond_overview, by = "specimen"
  ) %>%
  dplyr::select(specimen, condition, dplyr::everything())


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
on_tar_dists <- matched_algns %>%
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
    specimen, target, condition, 
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

sites_included <- on_tar_dists %>% 
  dplyr::group_by(condition, target) %>%
  dplyr::summarise(
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
tbl_ft_algn <- algnmts %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::filter(!edit.site %in% expandPosStr(on_targets)) %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(ft_algns = dplyr::n_distinct(clus.ori)) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# Probable edit sites
tbl_ft_prob <- probable_algns %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::filter(on.off.target == "Off-target") %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(ft_prob = dplyr::n_distinct(clus.ori)) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# Pile ups
tbl_ft_pile <- pile_up_algns %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::filter(on.off.target == "Off-target") %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(ft_pile = dplyr::n_distinct(clus.ori)) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# Paired or flanked loci
tbl_ft_pair <- paired_regions %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::filter(on.off.target == "Off-target") %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(ft_pair = n()) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# target sequence matched
tbl_ft_match <- matched_summary %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::filter(on.off.target == "Off-target") %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(ft_match = n()) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# Off-target summary table
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

# Evaluation summary ----
ot_eff_range <- tbl_ot_eff %>%
  tidyr::gather(key = "target", value = "eff", -specimen, -condition) %>%
  dplyr::group_by(specimen, condition) %>%
  dplyr::summarise(
    min = round(min(eff, na.rm = TRUE), digits = 1),
    max = round(max(eff, na.rm = TRUE), digits = 1),
    eff_rg = ifelse(
      min == max,
      sprintf("%.1f%%", max),
      sprintf("%1$.1f - %2$.1f%%", min, max)
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(eff_rg = ifelse(grepl("Inf", eff_rg), NA, eff_rg)) %>%
  dplyr::select(-min, -max)

eval_summary <- ot_eff_range %>%
  dplyr::left_join(
    ft_tbl_summary, by = c("specimen", "condition")
  ) %>%
  dplyr::left_join(
    tbl_algn_counts, by = "specimen"
  ) %>%
  dplyr::select(specimen, condition, Alignments, eff_rg, ft_match) %>%
  dplyr::rename(
    Specimen = specimen, Condition = condition, 
    "On-target\nEfficiency" = eff_rg, "Predicted\nOff-targets" = ft_match
  )

## Onco-gene enrichment analysis ----
rand_sites <- selectRandomSites(
  num = nrow(paired_regions) + nrow(matched_summary), 
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

ref_df <- data.frame(
  condition = "--",
  "total" = length(unique(ref_genes$annot_sym)),
  "onco" = sum(unique(onco_genes) %in% ref_genes$annot_sym),
  "special" = sum(unique(special_genes) %in% ref_genes$annot_sym)
)

paired_list <- split(
  x = paired_regions, 
  f = cond_overview$condition[
    match(paired_regions$specimen, cond_overview$specimen)
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
  x = matched_summary, 
  f = cond_overview$condition[
    match(matched_summary$specimen, cond_overview$specimen)
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
      "Reference" = ref_df, 
      "Flanking Pairs" = paired_df, 
      "Target Matched" = matched_df), 
    .id = "origin"
  ) %>%
  dplyr::filter(total > 0)

enrich_df$onco.p.value <- p.adjust(
  sapply(
    seq_len(nrow(enrich_df)), 
    function(i){
      
      ref <- enrich_df[1, c("total", "onco"), drop = TRUE]
      query <- enrich_df[i, c("total", "onco"), drop = TRUE]
      ref$diff <- abs(diff(as.numeric(ref)))
      query$diff <- abs(diff(as.numeric(query)))
      
      mat <- matrix(
        c(
          ref$total - ref$onco - query$total + query$onco,
          query$diff, ref$onco - query$onco, query$onco
        ),
        nrow = 2
      )
      
      fisher.test(mat)$p.value
      
    }
  ), 
  method = "BH"
)

enrich_df$special.p.value <- p.adjust(
  sapply(
    seq_len(nrow(enrich_df)), 
    function(i){
      
      ref <- enrich_df[1, c("total", "special"), drop = TRUE]
      query <- enrich_df[i, c("total", "special"), drop = TRUE]
      ref$diff <- abs(diff(as.numeric(ref)))
      query$diff <- abs(diff(as.numeric(query)))
      
      mat <- matrix(
        c(
          ref$total - ref$special - query$total + query$special,
          query$diff, ref$special - query$special, query$special
        ),
        nrow = 2
      )
      
      fisher.test(mat)$p.value
      
    }
  ), 
  method = "BH"
)

enrich_df <- enrich_df %>%
  dplyr::mutate(
    onco.power = sapply(seq_len(n()), function(i){

      statmod::power.fisher.test(
        p1 = onco[1] / total[1],
        p2 = onco[i] / total[i],
        n1 = total[1], n2 = total[i]
      )

    }),
    special.power = sapply(seq_len(n()), function(i){

      statmod::power.fisher.test(
        p1 = special[1] / total[1],
        p2 = special[i] / total[i],
        n1 = total[1], n2 = total[i]
      )

    })
  ) %>%
  dplyr::select(
    origin, condition, total, 
    onco, onco.p.value, onco.power, 
    special, special.p.value, special.power
  ) 

## Off-target sequence analysis ----
ft_MESL <- matched_algns %>%
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
        z = edit.site.dist, density = on_tar_dens[[condition]]
      ),
      gene_id = matched_summary$gene_id[
        match(edit.site, matched_summary$edit.site)
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

ft_seqs <- matched_summary %>%
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
      "input_vc" = vc_check,
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
    "incorp_data" = list(
      "algnmts" = algnmts, 
      "probable_algns" = probable_algns,
      "matched_algns" = matched_algns,
      "matched_summary" = matched_summary,
      "paired_algns" = paired_algns,
      "paired_regions" = paired_regions,
      "pile_up_algns" = pile_up_algns,
      "pile_up_summary" = pile_up_summary
    ), 
    "summary_tbls" = list(
      "ot_tbl_summary" = ot_tbl_summary,
      "ot_eff_summary" = tbl_ot_eff,
      "ft_tbl_summary" = ft_tbl_summary,
      "eval_summary" = eval_summary
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
  
  stop("\n  Cannot verify existence of output file:\n  ", args$output, "\n")
  
}else{
  
  if( !args$quiet ){
    cat("Evaluation complete, output writen to:\n  ", args$output, "\n")
  }
  
  q(status = 0)
  
}

