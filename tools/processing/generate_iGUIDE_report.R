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
  help = "Output report file, extension not required.")
parser$add_argument(
  "-c", "--config", nargs = "+", type = "character",
  help = "Run specific config file(s) in yaml format.")
parser$add_argument(
  "-s", "--support", nargs = "+", type = "character",
  help = "Supplementary data input, csv or tsv format.")
parser$add_argument(
  "-f", "--figures", action = "store_true",
  help = "Generate figures along with output report (pdf and png formats).")
parser$add_argument(
  "-d", "--data", action = "store_true",
  help = "Data to generate the report will be saved as an R image with output.")
parser$add_argument(
  "-t", "--format", nargs = 1, type = "character", default = "html",
  help = "Output format for report. Either 'pdf' or 'html' (default).")

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

if(length(args$input) != length(args$config)){
  stop("Must supply one config file for each input.")
}

report_formats <- c("html" = "html_document", "pdf" = "pdf_document")
if(!args$format %in% names(report_formats)){
  stop("Please input either 'html' or 'pdf' for format.\n",
       "Other formats not supported.")
}
output_format <- report_formats[args$format]

input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(1:length(args), function(i){
    paste(args[[i]], collapse = ", ")}))
input_table <- input_table[
  match(
    c("input :", "output :", "config :", "support :", 
      "figures :", "data :", "format :"),
    input_table$Variables),]
pandoc.title("iGUIDE Report Inputs:")
pandoc.table(data.frame(input_table, row.names = NULL), 
             justify = c("left", "left"), 
             split.tables = Inf)


# Load dependancies ------------------------------------------------------------
message("Loading dependencies.")
add_packs <- c(
  "stringr", "magrittr", "tidyverse", "data.table", "Biostrings", 
  "GenomicRanges", "knitr")

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

source(file.path(code_dir, "post_process_support.R"))

## Supporting functions
generate_genomic_regions <- function(ref, res, drop.alt.chr = TRUE){
  if(class(ref) == "BSgenome") ref <- GenomicRanges::seqinfo(ref)
  if(drop.alt.chr){
    ref <- ref[names(ref)[which(!stringr::str_detect(names(ref), "_"))]]
  }
  
  ref_len <- GenomeInfoDb::seqlengths(ref)
  
  unlist(GenomicRanges::GRangesList(lapply(seq_along(ref_len), function(i){
    seq <- names(ref_len)[i]
    starts <- seq(1, ref_len[seq], res)
    if(res > ref_len[seq]){
      ends <- ref_len[seq]
    }else{
      ends <- seq(res, ref_len[seq], res)
    }
    if(length(starts) - length(ends) == 1) ends <- c(ends, ref_len[seq])
    GenomicRanges::GRanges(
      seqnames = seq, 
      ranges = IRanges(start = starts, end = ends), 
      strand = "*",
      seqinfo = ref)
  })))
}

genomic_density <- function(gr, res, cutoff = 2, drop.alt.chr = TRUE){
  stopifnot(class(gr) == "GRanges")
  if(is.na(sum(is.numeric(
    GenomeInfoDb::seqlengths(GenomicRanges::seqinfo(gr)))))){
    stop("SeqInfo should be present for input GRanges.")
  }
  
  ref_regions <- generate_genomic_regions(
    GenomicRanges::seqinfo(gr), res, drop.alt.chr = drop.alt.chr)
  ref_regions$count <- GenomicRanges::countOverlaps(ref_regions, gr)
  ref_regions$log.count <- log(ref_regions$count, base = 10)
  ref_regions$norm.log.count <- ref_regions$log.count/max(ref_regions$log.count)
  ref_regions <- ref_regions[ref_regions$count >= cutoff]
  ref_regions
}

vcollapse <- function(d, sep, fill = "NA"){
  if(is.vector(d)){
    stop("Function vcollapse() is not used on vectors, use paste(collapse = ...).")
  }
  if(class(d) != "matrix") d <- as.matrix(d)
  mat <- d
  if(!is.null(fill)) mat <- ifelse(is.na(mat), fill, mat)
  mat <- do.call(
    cbind, 
    lapply(seq_len(ncol(mat)), function(i){
      if(i < ncol(mat)){
        cbind(mat[,i], rep(sep, nrow(mat)))
      }else{
        mat[,i]
      } 
    }))
  mat <- cbind(mat, rep(">|<", nrow(mat)))
  div_str <- stringr::str_c(t(mat), collapse = "")
  unlist(strsplit(div_str, split = ">\\|<"))
}

plot_genomic_density <- function(grl, res, grp.col = NULL, cutoff = 2, 
                                 drop.alt.chr = TRUE, clean = FALSE){
  if(class(grl) == "GRanges") grl <- GenomicRanges::GRangesList(grl)
  if(!is.null(grp.col)){
    stopifnot(all(unlist(lapply(grl, function(x){
      grp.col %in% names(GenomicRanges::mcols(x)) }))))
  }
  
  ref_len <- GenomeInfoDb::seqlengths(GenomicRanges::seqinfo(grl))
  ref_cum_len <- structure(
    c(0, cumsum(as.numeric(ref_len))), names = c("start", names(ref_len)))
  
  gen_densities <- lapply(grl, function(x){
    if(is.null(grp.col)){
      x <- list(x)
    }else{
      x <- split(x, vcollapse(
        GenomicRanges::mcols(x)[,grp.col, drop = FALSE], sep = " - "))
    }
      
    dplyr::bind_rows(lapply(x, function(y){
      y <- genomic_density(
        y, res = res, cutoff = cutoff, drop.alt.chr = drop.alt.chr)
      df <- as.data.frame(y, row.names = NULL) 
      cum_adj_pos <- ref_cum_len[match(df$seqnames, names(ref_cum_len)) - 1]
      df$adj.start <- cum_adj_pos + df$start
      df$adj.end <- cum_adj_pos + df$end
      df
    }), .id = "cond")
  })
  
  gen_den <- dplyr::bind_rows(gen_densities, .id = "grp")
  if(!is.null(names(grl))){
    gen_den$type <- factor(names(grl)[
      match(gen_den$grp, names(grl))], levels = names(grl))
  }else{
    gen_den$type <- factor(" ")
  }
  
  gen_den$score <- as.numeric(gen_den$type) + gen_den$norm.log.count
  
  # Grid layout
  x_breaks <- ref_cum_len[
    1:max(match(unique(gen_den$seqnames), names(ref_cum_len)))]
  x_lab_pos <- structure(
    sapply(2:length(x_breaks), function(i){
      mean(c(x_breaks[i-1], x_breaks[i])) }),
    names = names(x_breaks)[2:length(x_breaks)])
  
  y_breaks <- seq_along(grl)
  
  p <- ggplot(gen_den) 
  
  if(!clean){
    p <- p +
      geom_hline(yintercept = y_breaks, color = "grey90") +
      geom_vline(xintercept = x_breaks, color = "grey90") +
      scale_x_continuous(
        breaks = x_lab_pos,
        labels = gsub("chr", "", names(x_lab_pos)))
  }else{
    p <- p + scale_x_continuous(labels = NULL)
  }
  
  p <- p + 
    geom_rect(
      aes(
        xmin = adj.start, xmax = adj.end, 
        ymin = as.integer(type), ymax = score, 
        fill = type)) +
    scale_y_continuous(
      limits = c(0, length(grl)+1), breaks = seq_along(grl), labels = NULL) +
    scale_fill_brewer(type = "qual", palette = "Set1", direction = -1) +
    labs(x = "Chromosome", fill = "Levels") + 
    coord_polar() +
    theme_bw() +
    theme(
      panel.background = element_rect(color = "white"),
      panel.border = element_rect(color = "white"),
      panel.grid = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks = element_blank())
  
  if(clean){
    p + theme(legend.position = "none", axis.title = element_blank())
  }else{
    p
  }
}

div_seq <- function(seqs, ref, match.chr = "."){
  seqs <- as.character(seqs)
  ref <- as.character(ref)
  
  stopifnot(all(nchar(seqs) == nchar(ref)))
  
  seq_mat <- stringr::str_split(seqs, pattern = "", simplify = TRUE)
  ref_splt <- unlist(stringr::str_split(ref, pattern = ""))
  div_mat <- sapply(seq_along(ref_splt), function(i){
    ifelse(seq_mat[,i] == ref_splt[i], match.chr, seq_mat[,i]) })
  div_str <- stringr::str_c(t(div_mat), collapse = "")
  unlist(strsplit(
    div_str, split = paste0("(?<=.{", nchar(ref), "})"), perl = TRUE))
}

seq_diverge_plot <- function(df, ref, nuc.col = NULL, padding = 4, 
                             text.size = 2, convert.seq = TRUE, 
                             force.sq = FALSE, font.family = "Courier",
                             font.face = "bold"){
  if(is.null(nuc.col)){ nuc.col <- names(df)[1] }
  seqs <- dplyr::pull(df, var = match(nuc.col, names(df)))
  stopifnot(all(nchar(as.character(seqs)) == nchar(ref)))
  nuc_len <- nchar(ref)
  
  # Convert seqs
  if(convert.seq) seqs <- div_seq(seqs, ref)
  
  # Nucleotide color
  nucleotide_levels <- c("A", "T", "G", "C", ".", "N")
  nucleotide_colors <- RColorBrewer::brewer.pal(6, "Set1")
  nucleotide_colors <- c(nucleotide_colors[c(1,3,6,2)], "#FFFFFF", "#DCDCDC")
  names(nucleotide_colors) <- nucleotide_levels
  
  # Sequence matrix
  nuc_melt <- stringr::str_split(
    c(ref, paste(rep(" ", nuc_len), collapse = ""), seqs), 
    pattern = "", simplify = TRUE) %>%
    as.data.frame() %>%
    mutate(pos.y = -(1:n())) %>%
    melt(id.vars = "pos.y") %>%
    mutate(
      pos.x = as.numeric(stringr::str_extract(variable, "[0-9]+$")),
      color = nucleotide_colors[value],
      color = ifelse(is.na(color), "#FFFFFF", color)) %>%
    select(pos.x, pos.y, value, color)
  
  # Format remaining cols of input
  sup_df <- df[,-match(nuc.col, names(df))] %>%
    mutate_all(format, big.mark = ",", justify = "centre")
  
  sup_names <- names(sup_df)
  
  sup_df <- bind_rows(
    as.data.frame(
      t(matrix(
        c(names(sup_df), rep(" ", ncol(sup_df))), 
        ncol = 2, dimnames = list(names(sup_df))))),
    sup_df) %>%
    mutate_all(format, justify = "centre") %>%
    mutate(pos.y = -(1:n()))
  
  sup_melt <- melt(sup_df, id.vars = "pos.y") %>%
    mutate(
      pos.x = nuc_len + (match(variable, names(sup_df))) * padding,
      color = "#FFFFFF") %>%
    select(pos.x, pos.y, value, color)
  
  plot_melt <- bind_rows(nuc_melt, sup_melt)
  
  plot_colors <- structure(
    plot_melt$color[match(unique(plot_melt$value), plot_melt$value)],
    names = unique(plot_melt$value))
  
  p <- ggplot(plot_melt, aes(x = pos.x, y = pos.y)) +
    geom_tile(aes(fill = value)) +
    geom_text(
      aes(label = value), size = text.size, 
      family = font.family, fontface = font.face) +
    scale_fill_manual(values = plot_colors) +
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "none")
  
  if(force.sq){
    p <- p + theme(aspect.ratio = with(plot_melt, max(abs(pos.y))/max(pos.x)))
  }
  
  p
}


# Import metadata and consolidate into report objects --------------------------
message("Importing experimental data and configurations.")
## Load config files
configs <- lapply(args$config, yaml.load_file)
names(configs) <- sapply(configs, "[[", "Run_Name")

## Load reference genome 
if(grepl(".fa", unique(sapply(configs, "[[", "RefGenome")))){
  if(!file.exists(config$RefGenome)){
    stop("Specified reference genome file not found.")}
  ref_file_type <- ifelse(
    grepl(".fastq", unique(sapply(configs, "[[", "RefGenome"))), 
    "fastq", "fasta")
  ref_genome <- readDNAStringSet(
    unique(sapply(configs, "[[", "RefGenome")), format = ref_file_type)
}else{
  RefGenome <- unique(sapply(configs, "[[", "RefGenome"))
  genome <- grep(
    RefGenome, unique(BSgenome::installed.genomes()), value = TRUE)
  if(length(genome) == 0){
    pandoc.strong("Installed genomes include")
    pandoc.list(unique(BSgenome::installed.genomes()))
    stop("Selected reference genome not in list.")
  }else if(length(genome) > 1){
    pandoc.strong("Installed genomes include")
    pandoc.list(unique(BSgenome::installed.genomes()))
    stop(
      "Please be more specific about reference genome. Multiple matches to input.")
  }
  suppressMessages(library(genome, character.only = TRUE))
  ref_genome <- get(genome)
}

## Load reference files
ref_genes <- load_ref_files(
  configs[[1]]$refGenes, 
  type = "GRanges", freeze = configs[[1]]$RefGenome)
onco_genes <- load_ref_files(
  configs[[1]]$oncoGeneList, 
  type = "gene.list", freeze = configs[[1]]$RefGenome)
special_genes <- load_ref_files(
  configs[[1]]$specialGeneList, 
  type = "gene.list", freeze = config[[1]]$RefGenome)

umitag_option <- all(unlist(lapply(configs, "[[", "UMItags")))

upstream_dist <- unique(sapply(configs, function(x) x$upstreamDist))
downstream_dist <- unique(sapply(configs, function(x) x$downstreamDist))
if(length(upstream_dist) > 1 | length(downstream_dist) > 1){
  stop(
    "Inconsistant upstream or downstream distances between config files.\n",
    "Comparisons between groups with different run specific criteria\n", 
    "is not recommended.")
}

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
    pam_mat <- matrix(unlist(lapply(pam, function(pat){
        stringr::str_detect(seqs, paste0(pat,"$")) })), 
      ncol = length(pam))
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
names(on_targets) <- stringr::str_extract(
    names(on_targets), "[\\w\\_\\-\\']+$") %>%
  stringr::str_extract("[\\w\\_\\-]+")
on_targets <- structure(
  unique(on_targets), 
  names = names(on_targets)[match(unique(on_targets), on_targets)])

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
}else{
  stop(
    "Treatment information not accurately parsed from config(s).\n", 
    "Check config(s) formating.")
}

## Load in supporting information ==============================================
if(length(args$support) > 0){
  supp_data <- data.table::fread(args$support, data.table = FALSE)
  supp_data <- filter(supp_data, specimen %in% sample_info$specimen)
}else{
  supp_data <- data.frame()
}


## Consolidate supplementary data ==============================================
if(is.null(args$support)){
  spec_overview <- treatment_df
}else{
  spec_overview <- supp_data
}

cond_overview <- spec_overview %>%
  mutate(condition = vcollapse(
    select(spec_overview, -specimen), " - ", fill = "NA")) %>%
  select(specimen, condition)


## Read in experimental data and contatenate different sets ====================
input_data <- lapply(args$input, readRDS)
names(input_data) <- names(configs)
data_type_names <- names(input_data[[1]])
input_data <- lapply(1:length(input_data[[1]]), function(i){
  bind_rows(lapply(input_data, "[[", i), .id = "run.set")
})
names(input_data) <- data_type_names


message("Starting analysis.")
# Opening graphic --------------------------------------------------------------
graphic_order <- c("algnmts", "pile_up_algns", "paired_algns", "matched_algns")
graphic_data <- input_data[graphic_order]
graphic_grl <- GRangesList(lapply(
  graphic_data, 
  makeGRangesFromDataFrame, 
  seqinfo = seqinfo(ref_genome)))


# Specimen summary -------------------------------------------------------------
# Summarize components and append to specimen table
tbl_algn_counts <- input_data$algnmts %>% group_by(specimen)

if(umitag_option){
  tbl_algn_counts <- summarise(
    tbl_algn_counts, 
    Reads = sum(count), UMItags = sum(umitag), Alignments = n())
}else{
  tbl_algn_counts <- summarise(
    tbl_algn_counts, Reads = sum(count), Alignments = n())
}

spec_overview_join <- dplyr::left_join(
  spec_overview, tbl_algn_counts, by = "specimen") 


# On-target summary ------------------------------------------------------------
# Algnmts
tbl_ot_algn <- input_data$algnmts %>%
  mutate(specimen = factor(
    specimen, levels = sort(unique(sample_info$specimen)))) %>%
  group_by(specimen) %>%
  summarise(
    ot_algns = pNums(sum(as.integer(edit.site %in% on_targets))),
    ot_algns_pct = 100 * sum(as.integer(edit.site %in% on_targets))/n()) %>%
  ungroup() %>% as.data.frame()

# Probable edited sites
tbl_ot_prob <- input_data$probable_algns %>% 
  mutate(specimen = factor(
    specimen, levels = sort(unique(sample_info$specimen)))) %>%
  group_by(specimen) %>%
  summarise(
    ot_prob = pNums(sum(as.integer(edit.site %in% on_targets))),
    ot_prob_pct = 100 * sum(as.integer(edit.site %in% on_targets))/n()) %>%
  ungroup() %>% as.data.frame()

# Pile ups of read alignments
tbl_ot_pile <- input_data$pile_up_algns %>% 
  mutate(specimen = factor(
    specimen, levels = sort(unique(sample_info$specimen)))) %>%
  group_by(specimen) %>%
  summarise(
    ot_pile = pNums(sum(as.integer(edit.site %in% on_targets))),
    ot_pile_pct = 100 * sum(as.integer(edit.site %in% on_targets))/n()) %>%
  ungroup() %>% as.data.frame()

# Paired or flanking algnments
tbl_ot_pair <- input_data$paired_regions %>%
  mutate(
    specimen = factor(
      specimen, levels = sort(unique(sample_info$specimen))),
    on.off.target = factor(
      on.off.target, levels = c("On-target", "Off-target"))) %>%
  group_by(specimen, on.off.target) %>%
  summarise(cnt = sum(algns)) %>%
  group_by(specimen) %>%
  summarise(
    ot_pair = pNums(sum(ifelse(on.off.target == "On-target", cnt, 0))),
    ot_pair_pct = 100 * sum(ifelse(on.off.target == "On-target", cnt, 0)) /
      sum(cnt)) %>%
  ungroup() %>% as.data.frame()

# Guide RNA matched within 6 mismatches
tbl_ot_match <- input_data$matched_summary %>%
  mutate(specimen = factor(
    specimen, levels = sort(unique(sample_info$specimen)))) %>%
  group_by(specimen, on.off.target) %>%
  summarise(cnt = sum(algns)) %>%
  ungroup() %>% group_by(specimen) %>%
  summarise(
    ot_match = pNums(sum(ifelse(on.off.target == "On-target", cnt, 0))),
    ot_match_pct = 100 * sum(ifelse(on.off.target == "On-target", cnt, 0)) /
      sum(cnt)) %>%
  ungroup() %>% as.data.frame()

# Summary table
ot_tbl_summary <- mutate(treatment_df, specimen = factor(
  specimen, levels = sort(unique(sample_info$specimen))))

ot_tbl_summary <- Reduce(
  function(x,y){ dplyr::left_join(x, y, by = "specimen") },
  list(tbl_ot_algn[,c(1,3)], tbl_ot_pile[,c(1,3)],
       tbl_ot_pair[,c(1,3)], tbl_ot_match[,c(1,3)]),
  init = ot_tbl_summary)

names(ot_tbl_summary) <- c(
  "Specimen", "Treatment", "All\nAlign.", "Align.\nPileups", 
  "Flanking\nPairs", "gRNA\nMatched")


# On-target incorporation distribution -----------------------------------------
on_tar_dists <- input_data$matched_algns %>%
  filter(on.off.target == "On-target") %>%
  mutate(
    gRNA = str_extract(guideRNA.match, "[\\w]+"),
    pos = as.numeric(str_extract(edit.site, "[0-9]+$")),
    edit.site.dist = ifelse(strand == "+", start - pos, end - pos)) %>%
  dplyr::left_join(cond_overview, by = "specimen") %>%
  select(
    run.set, specimen, gRNA, condition, edit.site, edit.site.dist, strand) %>%
  group_by(condition, gRNA, edit.site.dist, strand) %>%
  summarise(cnt = n()) %>%
  ungroup() %>%
  mutate(
    strand.cnt = ifelse(
      strand == "+", log(cnt, base = 10), -log(cnt, base = 10)))

if(length(unique(on_tar_dists$condition)) == 1){
  on_tar_dists$condition <- " "
}

if(is.null(args$support)){
  sites_included <- on_tar_dists %>% group_by(gRNA)
}else{
  sites_included <- on_tar_dists %>% group_by(condition, gRNA)
}

sites_included <- summarise(
    sites_included,
    prop = 100 * sum(cnt[
      abs(edit.site.dist) <= upstream_dist & 
        abs(edit.site.dist) >= -downstream_dist]) / 
      sum(cnt),
    x_pos = upstream_dist,
    y_pos = 0.8 * min(strand.cnt[
      abs(edit.site.dist) <= upstream_dist & 
        abs(edit.site.dist) >= -downstream_dist])) %>%
  ungroup() %>%
  mutate(prop = paste0(pNums(prop, digits = 4), "%"))


# Off-target summary -----------------------------------------------------------
# All alignments
tbl_ft_algn <- input_data$algnmts %>%
  mutate(specimen = factor(
    specimen, levels = sort(unique(sample_info$specimen)))) %>%
  filter(!edit.site %in% on_targets) %>%
  group_by(specimen) %>%
  summarise(ft_algns = n_distinct(clus.ori)) %>%
  ungroup() %>% as.data.frame()

# Probable edit sites
tbl_ft_prob <- input_data$probable_algns %>%
  mutate(specimen = factor(
    specimen, levels = sort(unique(sample_info$specimen)))) %>%
  filter(on.off.target == "Off-target") %>%
  group_by(specimen) %>%
  summarise(ft_prob = n_distinct(clus.ori)) %>%
  ungroup() %>% as.data.frame()

# Pile ups
tbl_ft_pile <- input_data$pile_up_algns %>%
  mutate(specimen = factor(
    specimen, levels = sort(unique(sample_info$specimen)))) %>%
  filter(on.off.target == "Off-target") %>%
  group_by(specimen) %>%
  summarise(ft_pile = n_distinct(clus.ori)) %>%
  ungroup() %>% as.data.frame()

# Paired or flanked alignments
tbl_ft_pair <- input_data$paired_regions %>%
  mutate(specimen = factor(
    specimen, levels = sort(unique(sample_info$specimen)))) %>%
  filter(on.off.target == "Off-target") %>%
  group_by(specimen) %>%
  summarise(ft_pair = n()) %>%
  ungroup() %>% as.data.frame()

# gRNA sequence matched
tbl_ft_match <- input_data$matched_summary %>%
  mutate(specimen = factor(
    specimen, levels = sort(unique(sample_info$specimen)))) %>%
  filter(on.off.target == "Off-target") %>%
  group_by(specimen) %>%
  summarise(ft_match = n()) %>%
  ungroup() %>% as.data.frame()

# Summary table
ft_tbl_summary <- mutate(treatment_df, specimen = factor(
  specimen, levels = sort(unique(sample_info$specimen))))

ft_tbl_summary <- Reduce(
  function(x,y){ dplyr::left_join(x, y, by = "specimen") },
  list(tbl_ft_algn, tbl_ft_pile, tbl_ft_pair, tbl_ft_match),
  init = ft_tbl_summary)

names(ft_tbl_summary) <- c(
  "Specimen", "Treatment", "All\nAlign.", "Align.\nPileups", 
  "Flanking\nPairs", "gRNA\nMatched")


# Onco-gene enrichment analysis ------------------------------------------------
rand_sites <- selectRandomSites(
  num = nrow(input_data$paired_regions) + nrow(input_data$matched_summary), 
  refGenome = ref_genome, drop_extra_seqs = TRUE, setSeed = 714)

rand_sites$gene_id <- suppressMessages(assign_gene_id(
  seqnames(rand_sites), start(rand_sites), 
  reference = ref_genome, ref_genes = ref_genes, 
  onco_genes = onco_genes, special_genes = special_genes))

rand_df <- data.frame(
  condition = "Random Sites", 
  "total" = length(rand_sites), 
  "onco" = sum(stringr::str_detect(rand_sites$gene_id, "~")), 
  "special" = sum(stringr::str_detect(rand_sites$gene_id, "!")))

paired_list <- split(
  input_data$paired_regions, 
  cond_overview$condition[
    match(input_data$paired_regions$specimen, cond_overview$specimen)])

paired_df <- bind_rows(lapply(paired_list, function(df){
  data.frame(
    "total" = nrow(df), 
    "onco" = sum(stringr::str_detect(df$gene_id, "~")), 
    "special" = sum(stringr::str_detect(df$gene_id, "!")))
}), .id = "condition")

matched_list <- split(
  input_data$matched_summary, 
  cond_overview$condition[
    match(input_data$matched_summary$specimen, cond_overview$specimen)])

matched_df <- bind_rows(lapply(matched_list, function(df){
  data.frame(
    "total" = nrow(df), 
    "onco" = sum(stringr::str_detect(df$gene_id, "~")), 
    "special" = sum(stringr::str_detect(df$gene_id, "!")))
}), .id = "condition")

enrich_df <- bind_rows(list(
  "Reference" = rand_df, 
  "Flanking Pairs" = paired_df, 
  "gRNA Matched" = matched_df), .id = "origin")

enrich_df$onco.p.value <- p.adjust(sapply(seq_len(nrow(enrich_df)), function(i){
  ref <- enrich_df[1, c("total", "onco")]
  query <- enrich_df[i, c("total", "onco")]
  ref$diff <- abs(diff(as.numeric(ref)))
  query$diff <- abs(diff(as.numeric(query)))
  fisher.test(as.matrix(rbind(
    ref[,c("diff", "onco")], query[,c("diff", "onco")])))$p.value
}), method = "BH")

enrich_df$special.p.value <- p.adjust(sapply(seq_len(nrow(enrich_df)), function(i){
  ref <- enrich_df[1, c("total", "special")]
  query <- enrich_df[i, c("total", "special")]
  ref$diff <- abs(diff(as.numeric(ref)))
  query$diff <- abs(diff(as.numeric(query)))
  fisher.test(as.matrix(rbind(
    ref[,c("diff", "special")], query[,c("diff", "special")])))$p.value
}), method = "BH")

names(enrich_df) <- c(
  "Origin", "Condition", "Total", "Onco-Rel.", 
  "Spec. Genes", "Onco Enrich.", "Special Enrich.")


# Genomic Distribution of edited sites -----------------------------------------
genomic_grl <- GRangesList(lapply(graphic_data, function(x){
  y <- makeGRangesFromDataFrame(x, seqinfo = seqinfo(ref_genome))
  mcols(y) <- cond_overview[
    match(x$specimen, cond_overview$specimen), "condition", drop = FALSE]
  y
}))

num_conds <- max(length(unique(cond_overview$condition)), 1)

names(genomic_grl) <- c(
  "All Align.", "Pileup Align.", "Flanking Pairs", "gRNA Matched")


# Off-target sequence analysis -------------------------------------------------
ft_seqs <- input_data$matched_summary %>%
  select(
    specimen, aligned.sequence, guideRNA.match, edit.site,
    guideRNA.mismatch, on.off.target, algns, gene_id) %>% 
  dplyr::left_join(cond_overview, by = "specimen")

if(is.null(args$support)){
  ft_seqs <- group_by(
      ft_seqs, 
      guideRNA.match, edit.site, aligned.sequence, 
      guideRNA.mismatch, on.off.target, gene_id) %>%
    summarise(algns = sum(algns))
}else{
  ft_seqs <- group_by(
      ft_seqs,
      condition, guideRNA.match, edit.site, aligned.sequence, 
      guideRNA.mismatch, on.off.target, gene_id) %>%
    summarise(algns = sum(algns))
}

ft_seqs <- arrange(ft_seqs, desc(algns), guideRNA.mismatch) %>%
  ungroup() %>%
  mutate(on.off.target = stringr::str_extract(on.off.target, "[\\w]+")) %>%
  rename(
    on.off.target = "target", 
    guideRNA.mismatch = "mismatch", 
    guideRNA.match = "gRNA")

if(is.null(args$support)){
  ft_seqs_list <- split(ft_seqs, ft_seqs$gRNA)
}else{
  ft_seqs_list <- split(ft_seqs, paste0(ft_seqs$condition, " - ", ft_seqs$gRNA))
}

message("Analysis complete. Starting report generation.")

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

if(args$format == "html" & !str_detect(output_file, ".html$")){
  output_file <- paste0(output_file, ".html")
}

if(args$format == "pdf" & !str_detect(output_file, ".pdf$")){
  output_file <- paste0(output_file, ".pdf")
}

if(args$data){
  if(args$format == "html"){
    save.image(file = file.path(
      output_dir, stringr::str_replace(output_file, ".html$", ".RData"))) 
  }else if(args$format == "pdf"){
  save.image(file = file.path(
    output_dir, stringr::str_replace(output_file, ".pdf$", ".RData"))) 
  }
}

if(args$format == "html"){
  rmarkdown::render(
    input = file.path(code_dir, "iGUIDE_report_template.Rmd"),
    output_format = output_format, 
    output_file = output_file,
    output_dir = output_dir,
    output_options = list("css" = file.path(code_dir, "iguide.css")))
}else{
  rmarkdown::render(
    input = file.path(code_dir, "iGUIDE_report_template.Rmd"),
    output_format = output_format, 
    output_file = output_file,
    output_dir = output_dir)
}

q()
