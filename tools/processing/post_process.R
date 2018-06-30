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
panderOptions("table.style", "simple")
panderOptions("table.split.table", Inf)

# Set up and gather command line arguments -------------------------------------
parser <- ArgumentParser(
  description = "Post-processing script for iGUIDE.")
parser$add_argument(
  "uniqSites", nargs = 1, type = "character",
  help = "Unique sites output from blatCoupleR. The output from an entire run can be concatenated together as a single input.")
parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", 
  help = "Output file name in .rds format.")
parser$add_argument(
  "-c", "--config", nargs = 1, type = "character",
  help = "Run specific config file in yaml format.")
parser$add_argument(
  "-u", "--umitags", nargs = "+", type = "character",
  help = "Path(s) to associated fasta files containing read specific random captured sequences. Multiple file paths can be separated by a space.")
parser$add_argument(
  "-m", "--multihits", nargs = "+", type = "character",
  help = "Path(s) to associated multihit files (.rds) as produced by blatCoupleR. Multiple file paths can be separated by a space.")


args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(1:length(args), function(i){
    paste(args[[i]], collapse = ", ")}))
input_table <- input_table[
  match(c("uniqSites :", "output :", "config :", "umitags :", "multihits :"),
        input_table$Variables),]
pandoc.title("Post-processing Inputs")
pandoc.table(data.frame(input_table, row.names = NULL), 
             justify = c("left", "left"), 
             split.tables = Inf)

# Load dependancies ------------------------------------------------------------
add_packs <- c(
  "stringr", "magrittr", "dplyr", "Matrix", "igraph", "ShortRead",
  "Biostrings", "GenomicRanges", "BSgenome", "hiAnnotator")

add_packs_loaded <- suppressMessages(
  sapply(add_packs, require, character.only = TRUE))
if(!all(add_packs_loaded)){
  pandoc.table(data.frame(
    "R-Packages" = names(add_packs_loaded), 
    "Loaded" = add_packs_loaded, 
    row.names = NULL))
  stop("Check dependancies.")
}

# Source supporting functions --------------------------------------------------
code_dir <- dirname(
  sub("--file=", "", grep("--file=", 
                          commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(code_dir, "post_process_support.R"))

# Inputs and parameters --------------------------------------------------------
config <- yaml::yaml.load_file(args$config)
sample_info <- data.table::fread(
  file.path(config$Install_Directory, config$Sample_Info), data.table = FALSE)
submat <- banmat()

## Load reference genome =======================================================
if(grepl(".fa", config$RefGenome)){
  if(!file.exists(config$RefGenome)){
    stop("Specified reference genome file not found.")}
  ref_file_type <- ifelse(grepl(".fastq", config$RefGenome), "fastq", "fasta")
  ref_genome <- readDNAStringSet(config$RefGenome, format = ref_file_type)
}else{
  genome <- grep(
    config$RefGenome, unique(BSgenome::installed.genomes()), value = TRUE)
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

## Load refGenes and gene lists for annotation =================================
ref_genes <- loadRefFiles(
  config$refGenes, type = "GRanges", freeze = config$RefGenome)
onco_genes <- loadRefFiles(
  config$oncoGeneList, type = "gene.list", freeze = config$RefGenome)
special_genes <- loadRefFiles(
  config$specialGeneList, type = "gene.list", freeze = config$RefGenome)

## Incorporation site parameters ===============================================
upstream_dist <- config$upstreamDist
downstream_dist <- config$downstreamDist
maxGuideMismatch <- config$maxGuideMismatch
pile_up_min <- config$pileUpMin
on_target_sites <- config$On_Target_Sites 

## Load Guide RNAs and sample metadata =========================================
guide_rna_seqs <- lapply(config$Guide_RNA_Sequences, toupper)
pam_seq <- lapply(config$PAM_Sequence, toupper)
pam_mat <- matrix(unlist(lapply(pam_seq, function(pat){
    stringr::str_detect(guide_rna_seqs, paste0(pat,"$")) })), 
  ncol = length(pam_seq))
if(any(rowSums(pam_mat) > 1)){ 
  stop("Multiple PAM sequences detected on a single guide RNA.") }
rownames(pam_mat) <- names(guide_rna_seqs)
colnames(pam_mat) <- unlist(pam_seq)
gRNA_tbl <- data.frame(
  row.names = names(guide_rna_seqs),
  "Guide" = names(guide_rna_seqs),
  "gRNA" = sapply(1:length(guide_rna_seqs), function(i){
    stringr::str_replace(
      guide_rna_seqs[[i]], paste0(colnames(pam_mat)[pam_mat[i,]], "$"), "")}),
  "PAM" = colnames(pam_mat)[
    sapply(1:nrow(pam_mat), function(i) which(pam_mat[i,]))])

cat("Guide RNA Sequence Table:")
pandoc.table(
  gRNA_tbl, 
  row.names = FALSE,
  style = "simple", 
  split.table = Inf,
  justify = "center",
  emphasize.rownames = FALSE)

## Load data related to how samples were processed =============================
treatment <- config$Treatment
if(any(grepl("sampleInfo:", treatment[1]))){
  info_col <- match(str_extract(treatment[1], "[\\w]+$"), names(sample_info))
  if(length(info_col) != 1){
    stop("Cannot parse treatment data. Check config yaml and sampleInfo.")
  }
  treatment_df <- data.frame(
    sampleName = sample_info$sampleName, 
    treatment = sample_info[,info_col])
  treatment_df$specimen <- str_extract(treatment_df$sampleName, "[\\w]+")
  treatment_df <- unique(treatment_df[,c("specimen", "treatment")])
  treatment <- strsplit(treatment_df$treatment, ";")
  names(treatment) <- treatment_df$specimen
}else if(any(grepl("all", names(treatment)))){
  treatment_df <- data.frame(
    sampleName = sample_info$sampleName, 
    treatment = unique(unlist(treatment)))
  treatment_df$specimen <- str_extract(treatment_df$sampleName, "[\\w]+")
  treatment_df <- unique(treatment_df[,c("specimen", "treatment")])
  treatment <- strsplit(treatment_df$treatment, ";")
  names(treatment) <- treatment_df$specimen
}else{
  treatment_df <- data.frame(
    "specimen" = names(treatment), 
    "treatment" = sapply(treatment, paste, collapse = ";"))
}

cat("Sample Treatment Table:\n")
pandoc.table(
  t(as.data.frame(treatment)), 
  style = "simple", 
  split.table = Inf,
  justify = "center",
  emphasize.rownames = FALSE)


# Load input data --------------------------------------------------------------
## Unique sites ================================================================
reads <- data.table::fread(
  args$uniqSites, data.table = FALSE, stringsAsFactors = FALSE)

## Multihits if requested ======================================================
if(all(!is.null(args$multihits))){
  uniq_reads <- GenomicRanges::makeGRangesFromDataFrame(
    reads, keep.extra.columns = TRUE, seqinfo = seqinfo(ref_genome))
  
  multi_reads <- unlist(GRangesList(lapply(args$multihits, function(x){
    multi <- readRDS(x)
    seqinfo(multi$unclustered_multihits) <- seqinfo(ref_genome)
    if(length(multi$unclustered_multihits) > 0){
      mcols(multi$unclustered_multihits) <- mcols(multi$unclustered_multihits)[
        ,c(names(mcols(uniq_reads)))]
    }else{
      mcols(multi$unclustered_multihits) <- mcols(uniq_reads)[0
        ,c(names(mcols(uniq_reads)))]
    }
    multi$unclustered_multihits
  })))
  
  comb_reads <- c(uniq_reads, multi_reads)
  mcols(comb_reads)$type <- rep(
    c("uniq", "multi"), c(length(uniq_reads), length(multi_reads)))
  mcols(comb_reads)$clus.id <- pileupCluster(
    comb_reads, grouping = "sampleName", return = "ID")
  
  filt_multi_reads <- bind_rows(lapply(
    split(comb_reads, comb_reads$sampleName), function(x){
      uniq_id <- unique(x$clus.id[x$type == "uniq"])
      multi_id <- unique(x$clus.id[x$type == "multi"])
      y <- x[x$type == "multi" & x$clus.id %in% intersect(uniq_id, multi_id)]
      mcols(y)$clus.id <- NULL
      if(length(y) > 0){
        contrib_amt <- 1/table(mcols(y)$ID)
        mcols(y)$contrib <- as.numeric(contrib_amt[mcols(y)$ID])
      }
      as.data.frame(y, row.names = NULL) %>%
        mutate(seqnames = as.character(seqnames), strand = as.character(strand))
    }))
  
  reads <- mutate(reads, type = "uniq", contrib = 1) %>%
    bind_rows(., filt_multi_reads)
}else{
  reads <- mutate(reads, type = "uniq", contrib = 1)
}

# Print out stats during analysis.
cat("Tabulation of aligned reads per specimen:")
temp_table <- table(str_extract(reads$sampleName, "[\\w]+"))
pandoc.table(
  data.frame(
    "Specimen" = names(temp_table), 
    "Aligned_Reads" = format(as.numeric(temp_table), big.mark = ","),
    row.names = NULL),
  style = "simple",
  justify = "cr",
  emphasize.rownames = FALSE)
rm(temp_table)

## Umitags or captured random sequences ========================================
if(all(!is.null(args$umitags))){
  umitags <- lapply(args$umitags, ShortRead::readFasta)
  umitags <- serial_append_S4(umitags)
  reads$umitag <- as.character(ShortRead::sread(umitags))[
    match(reads$ID, as.character(ShortRead::id(umitags)))]
}

# Process input data -----------------------------------------------------------
## Format input alignments =====================================================
algnmts <- mutate(reads, specimen = str_extract(sampleName, "[\\w]+"))

if(config$UMItags){
  algnmts <- arrange(algnmts, desc(contrib)) %>%
    group_by(seqnames, start, end, strand, specimen, sampleName) %>%
    summarise(
      contrib = max(contrib),
      count = sum(contrib),
      umitag = as.integer(!duplicated(umitag)) * contrib) %>%
    ungroup() %>%
    select(
      seqnames, start, end, strand, specimen, 
      sampleName, contrib, umitag, count) %>%
    as.data.frame()
}else{
  algnmts <- group_by(
      algnmts, seqnames, start, end, strand, specimen, sampleName) %>%
    summarize(contrib = max(contrib), count = sum(contrib)) %>%
    ungroup() %>%
    select(
      seqnames, start, end, strand, specimen, sampleName, contrib, count) %>%
    as.data.frame()
}

sample_index <- ifelse(nrow(algnmts) > 10, 10, nrow(algnmts))
sample_index <- sample(1:nrow(algnmts), sample_index, replace = FALSE)

cat("Sample of aligned templates:")
pandoc.table(
  data.frame(algnmts[sample_index,], row.names = NULL),
  style = "simple", 
  caption = paste0("Number of templates: ", nrow(algnmts)),
  justify = "left")

rm(sample_index)

algnmts_gr <- GRanges(
  seqnames = algnmts$seqnames,
  ranges = IRanges(start = algnmts$start, end = algnmts$end),
  strand = algnmts$strand,
  seqinfo = seqinfo(ref_genome))

if(config$UMItags){
  mcols(algnmts_gr) <- dplyr::select(
    algnmts, specimen, sampleName, contrib, umitag, count)
}else{
  mcols(algnmts_gr) <- dplyr::select(
    algnmts, specimen, sampleName, contrib, count)
}

# Analyze alignments -----------------------------------------------------------
# Identify groups of alignments or pileups of aligned fragments
# These pileups give strong experimental evidence of directed incorporation of
# the dsODN into a region. Initially, pileups are identified and then checked 
# for pairing, or if there is another pileup on the opposite strand in close 
# proximity.
algnmts_gr$clus.ori <- pileupCluster(algnmts_gr, "specimen", 0L, "simple")
algnmts_gr$paired.algn <- identifyPairedAlgnmts(
  algnmts_gr, upstream_dist*2, "specimen")

# Create a GRange with only the unique cluster origins
split_clus_id <- stringr::str_split(
  unique(algnmts_gr$clus.ori), ":", simplify = TRUE)
algn_clusters <- GRanges(
  seqnames = split_clus_id[,1],
  ranges = IRanges(start = as.numeric(split_clus_id[,3]), width = 1),
  strand = split_clus_id[,2],
  seqinfo = seqinfo(ref_genome))

algn_clusters$clus.ori <- unique(algnmts_gr$clus.ori)
algn_clusters$clus.seq <- getSiteSeqs(
  algn_clusters, upstream_flank = upstream_dist, 
  downstream_flank = downstream_dist, ref_genome = ref_genome)

# Identify which guideRNAs potentially bind near clusters
algn_clusters <- compareGuideRNAs(
  algn_clusters, guideRNASeqs = guide_rna_seqs, 
  submat = submat, seq_col = "clus.seq", tolerance = maxGuideMismatch,
  upstream_flank = upstream_dist, downstream_flank = downstream_dist)

# Merge the guideRNA alignment information from the clusters back to all unique
# alignments
algnmts <- merge(
  as.data.frame(algnmts_gr), 
  mcols(algn_clusters)[,c(
    "clus.ori", "guideRNA.match", "guideRNA.mismatch", 
    "aligned.sequence", "edit.site")], 
  by = "clus.ori") %>%
  as.data.frame()

# Change guideRNA.match to No_Valid_Match if an inappropriate gRNA is annotated
algnmts$guideRNA.match <- filterInappropriateComparisons(
  algnmts$guideRNA.match, algnmts$specimen, treatment)

# Fragment pileups, paired clustering, and guideRNA alignments have been used to
# characterize the incorporation sites analyzed here. Each metric will be used 
# to create a list of incorporation sites that may be CRISPR cut sites. The 
# following identifies which alignments are associated with each of these 
# criteria.
tbl_clus_ori <- algnmts %>% 
  group_by(specimen, clus.ori) %>%
  filter(n() >= pile_up_min) %>%
  ungroup() %$%
  table(clus.ori)
idx_clus_ori <- which(algnmts$clus.ori %in% names(tbl_clus_ori))

tbl_paried_algn <- algnmts %>%
  filter(!is.na(paired.algn)) %$%
  table(paired.algn)
idx_paired_algn <- which(algnmts$paired.algn %in% names(tbl_paried_algn))

idx_matched <- which(algnmts$guideRNA.match != "No_valid_match")

idx_combined <- sort(unique(c(idx_clus_ori, idx_paired_algn, idx_matched)))

idx_df <- data.frame(
  "Type" = c("PileUp", "Paired", "gRNA_Matched", "Combined"),
  "Counts" = sapply(
    list(idx_clus_ori, idx_paired_algn, idx_matched, idx_combined), length))
cat("Table of uniquely aligned template counts:")
pandoc.table(
  idx_df,
  style = "simple", 
  caption = paste0("Total number of templates: ", nrow(algnmts)),
  justify = "left")

probable_algns <- algnmts[idx_combined,]
probable_algns$on.off.target <- ifelse(
  probable_algns$edit.site %in% unlist(on_target_sites), 
  "On-target", "Off-target")
cat("On / Off target alignment counts.")
print(table(probable_algns$on.off.target))

# Create summary and output formated object related to each of the criteria for
# edited site detection.
matched_algns <- probable_algns[
  probable_algns$guideRNA.match != "No_valid_match",]

matched_summary <- matched_algns %>%
  mutate(guideRNA.match = stringr::str_extract(guideRNA.match, "[\\w]+")) %>%
  group_by(
    specimen, edit.site, aligned.sequence, 
    guideRNA.match, guideRNA.mismatch)

if(config$UMItags){
  matched_summary <- summarise(
    matched_summary,
    on.off.target = paste(sort(unique(on.off.target)), collapse = ";"),
    paired.algn = paste(sort(unique(paired.algn)), collapse = ";"),
    umitag = sum(umitag),
    count = sum(count), 
    algns = sum(contrib),
    orient = paste(sort(unique(as.character(strand))), collapse = ";"))
}else{
  matched_summary <- summarise(
    matched_summary,
    on.off.target = paste(sort(unique(on.off.target)), collapse = ";"),
    paired.algn = paste(sort(unique(paired.algn)), collapse = ";"),
    count = sum(count), 
    algns = sum(contrib),
    orient = paste(sort(unique(as.character(strand))), collapse = ";"))
}

matched_summary <- ungroup(matched_summary) %>% 
  as.data.frame() %>%
  mutate(gene_id = assignGeneID( 
    seqnames = stringr::str_extract(edit.site, "[\\w]+"), 
    positions = as.numeric(stringr::str_extract(edit.site, "[\\w]+$")), 
    reference = ref_genome, 
    ref_genes = ref_genes, onco_genes = onco_genes, 
    special_genes = special_genes))

paired_algns <- probable_algns[
  probable_algns$paired.algn %in% names(tbl_paried_algn),]

paired_regions <- paired_algns %>%
  group_by(specimen, paired.algn, strand) %>%
  mutate(pos = ifelse(strand == "+", min(start), max(end))) %>%
  group_by(specimen, paired.algn)

if(config$UMItags){
  paired_regions <- summarise(
      paired_regions,
      seqnames = unique(seqnames),
      start = min(pos), end = max(pos), mid = start + (end-start)/2,
      strand = "*", width = end - start, algns = sum(contrib), 
      umitag = sum(umitag), count = sum(count)) %>%
    ungroup()
}else{
  paired_regions <- summarise(
      paired_regions,
      seqnames = unique(seqnames),
      start = min(pos), end = max(pos), mid = start + (end-start)/2,
      strand = "*", width = end - start, algns = sum(contrib), 
      count = sum(count)) %>%
    ungroup()
}

paired_regions <- mutate(
    paired_regions,     
    gene_id = assignGeneID(
      seqnames, mid, reference = ref_genome, 
      ref_genes = ref_genes, onco_genes = onco_genes, 
      special_genes = special_genes)) %>%
  group_by(specimen, paired.algn) %>%
  mutate(
    on.off.target = ifelse(
      any(sapply(
        unlist(on_target_sites[
          which(
            str_extract(names(on_target_sites), "[\\w\\-\\_]+") %in% 
              treatment[[specimen]])]),
        function(x, seq, st, en){
          seq == str_extract(x, "[\\w]+") &
            st <= as.numeric(str_extract(x, "[\\w]+$")) + downstream_dist &
            en >= as.numeric(str_extract(x, "[\\w]+$")) - downstream_dist
      }, seq = seqnames, st = start, en = end)), "On-target", "Off-target")) %>%
  ungroup() %>% 
  as.data.frame()
      
pile_up_algns <- probable_algns[
  probable_algns$clus.ori %in% names(tbl_clus_ori),]

# Output data composition ------------------------------------------------------
# rds file that can be read into reports or loaded
# into a data base with some additional scripting.
data_comp <- list(
  "algnmts" = algnmts, 
  "probable_algns" = probable_algns,
  "matched_algns" = matched_algns,
  "matched_summary" = matched_summary,
  "paired_algns" = paired_algns,
  "paired_regions" = paired_regions,
  "pile_up_algns" = pile_up_algns)

saveRDS(data_comp, file = args$output)
if(file.exists(args$output)){
  message("Successfully completed script.")
}else{
  message("Check output, not detected at the end of post-processing script.")
}
q()



# Previous code and development ------------------------------------------------

edit_sites <- split(probable_sites, probable_sites$edit.site)
mated_edit_sites <- edit_sites[
  sapply(edit_sites, function(x){
    all(c("+", "-") %in% names(table(as.character(strand(x)))))})]
mated_edit_sites <- unname(unlist(mated_edit_sites))
mated_edit_sites <- split(mated_edit_sites, mated_edit_sites$specimen)
mated_edit_sites <- lapply(
  mated_edit_sites, function(gr) split(gr, gr$edit.site))

mated_edit_loci_plots <- lapply(
  mated_edit_sites, function(x) lapply(
    x, plot_edit_sites, sampleName = "specimen", resolution = 1L))

edit_loci_plots <- lapply(edit_sites, plot_edit_sites, resolution = 1L)

crispr_sites <- GenomicRanges::as.data.frame(
  probable_sites, row.names = NULL) %>%
  select(specimen, count, guideRNA.match, guideRNA.mismatch, 
         aligned.sequence, edit.site) %>%
  mutate(guideRNA.match = str_extract(guideRNA.match, "[\\w\\.\\-]+")) %>%
  group_by(specimen, guideRNA.match, guideRNA.mismatch, 
           aligned.sequence, edit.site) %>%
  summarize(
    read_counts = sum(count), 
    template_counts = n()) %>%
  ungroup() %>%
  mutate(
    site.chr = sapply(strsplit(edit.site, ":"), "[[", 1),
    site.ort = sapply(strsplit(edit.site, ":"), "[[", 2),
    site.pos = as.numeric(sapply(strsplit(edit.site, ":"), "[[", 3))) %>%
  as.data.frame()
crispr_sites$gene_id <- assignGeneID(
  crispr_sites$site.chr, 
  positions = crispr_sites$site.pos, 
  reference = ref_genome, 
  ref_genes = ref_genes, 
  onco_genes = onco_genes, 
  special_genes = special_genes)
crispr_sites$gene_id <- sapply(1:nrow(crispr_sites), function(i, df, on){
  if(df$edit.site[i] %in% on){
    return(paste0(names(on)[match(df$edit.site[i], on)], "*"))
  }else{
    return(df$gene_id[i])
  }
}, df = crispr_sites, on = on_target_sites)

crispr_sites$on.off.target <- ifelse(
  crispr_sites$edit.site %in% on_target_sites, "On-target", "Off-target")

crispr_gr <- GRanges(
  seqnames = str_extract(crispr_sites$edit.site, "[\\w]+"),
  ranges = IRanges(
    start = as.numeric(str_extract(crispr_sites$edit.site, "[\\d]+$")),
    width = 1L),
  strand = str_extract(crispr_sites$edit.site, "[\\+\\-\\*]"))

crispr_gr <- doAnnotation(
  "within", crispr_gr, ref_genes, 
  colnam = "inGene", feature.colnam = "annot_sym")
crispr_gr <- doAnnotation(
  "within", crispr_gr, ref_exons, 
  colnam = "inExon", feature.colnam = "geneNames")
crispr_gr <- doAnnotation(
  "nearest", crispr_gr, ref_genes, 
  colnam = "nearestGene", feature.colnam = "annot_sym")
crispr_gr <- doAnnotation(
  "within", crispr_gr, oregs, 
  colnam = "inReg", feature.colnam = "id")
crispr_gr <- doAnnotation(
  "within", crispr_gr, oregs, 
  colnam = "inRegType", feature.colnam = "type")
crispr_gr <- doAnnotation(
  "nearest", crispr_gr, oregs, 
  colnam = "nearestReg", feature.colnam = "id")
crispr_gr <- doAnnotation(
  "nearest", crispr_gr, oregs, 
  colnam = "nearestRegType", feature.colnam = "type")
crispr_gr <- doAnnotation(
  "within", crispr_gr, oregs[oregs$type == "Transcription Factor Binding Site"],
  colnam = "inTFbs", feature.colnam = "TFbs")
crispr_gr <- doAnnotation(
  "nearest", crispr_gr, oregs[oregs$type == "Transcription Factor Binding Site"],
  colnam = "nearestTFbs", feature.colnam = "TFbs")
crispr_gr <- doAnnotation(
  "nearest", crispr_gr, ref_genes[ref_genes$annot_sym %in% onco_genes],
  colnam = "nearestOnco", feature.colnam = "annot_sym")

crispr_sites <- bind_cols(
  crispr_sites, 
  as.data.frame(mcols(crispr_gr)[,c(
    "inGene", "inExon", "nearestGene", "nearestGeneDist", "inReg", "inRegType",
    "nearestReg", "nearestRegDist", "nearestRegType", "inTFbs", "nearestTFbs",
    "nearestTFbsDist", "nearestOnco", "nearestOncoDist")]))




## Dev and notes ---------------------------------------------------------------
## Dev:
dev <- as.data.frame(algnmts_gr, row.names = NULL) %>%
  group_by(clusStart) %>%
  mutate(
    pos = as.numeric(stringr::str_extract(clusStart, "[0-9]+$")),
    dis = ifelse(strand == "+", start - pos, pos - end))

all_dev <- dev %>% ungroup() %>%
  mutate(grp = "all") %>% select(grp, dis)
on_tar <- dev %>% ungroup() %>% 
  filter(seqnames == "chr14", pos < 22547670, pos > 22547660) %>%
  mutate(grp = "on_tar") %>% select(grp, dis)
pile_ups <- dev %>% filter(n() > 1) %>% ungroup() %>%
  mutate(grp = "piles") %>% select(grp, dis)

bind_rows(all_dev, on_tar, pile_ups) %>%
  as.data.frame() %>%
  ggplot(aes(x = dis)) +
  geom_density(aes(color = factor(grp)), adjust = 2) +
  lims(x = c(-10, 100), y = c(0, 0.5)) +
  theme_bw()

# breakpoints may need to be refined prior to analysis
library(gintools)
g <- unname(unlist(GRangesList(lapply(
  split(algnmts_gr, algnmts_gr$sampleName), 
  refine_breakpoints, counts = "count"))))
gu <- g
gu$called.bp <- gu$adj.bp <- NULL
gu <- unique_granges(gu, sum.counts = TRUE, counts.col = "count")
