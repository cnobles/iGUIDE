# Global options ---------------------------------------------------------------
#' Adjust environmental options
options(stringsAsFactors = FALSE)
cores <- 6

# Load dependancies ------------------------------------------------------------
#' Load all required dependencies / packages for analysis
dependencies <- c("GenomicRanges", "stringr", "Biostrings", "Matrix", "ggplot2",
                  "parallel", "BSgenome", "Rsamtools", "dplyr", "pander", 
                  "knitr", "markdown", "igraph", "BSgenome.Hsapiens.UCSC.hg38",
                  "hiAnnotator")

null <- suppressMessages(sapply(dependencies, library, character.only = TRUE))
if(all(paste0("package:", dependencies) %in% search())){
  message("Dependencies loaded successfully.")
  rm(dependencies, null)
}else{
  pandoc.table(data.frame(
    "Dependencies" = dependencies,
    "Loaded" = paste0("package:", dependencies) %in% search()
  ))
  stop("Not all dependencies loaded. Check for issues.")
}

# Input parameters -------------------------------------------------------------
#' Input parameters, later will be moved to yaml or other input format.
## log table for guide sequences.
guideRNASeqs <- list(
  "NYCE.PD1-3" = c("ggcgccctggccagtcgtctngg"),
  "NYCE.TRAC-5" = c("tgtgctagacatgaggtctangg"),
  "NYCE.TRBC-4" = c("ggagaatgacgagtggacccngg"))

cat("Guide RNA Sequence Table:")
pandoc.table(
  t(as.data.frame(guideRNASeqs)), 
  style = "simple", 
  split.table = Inf,
  justify = "left",
  emphasize.rownames = FALSE)

## log table for which samples were treated with which guides.
treated_with <- list(
    "R1C4" = c("NYCE.PD1-3", "NYCE.TRAC-5", "NYCE.TRBC-4"),
    "R1C5" = c("NYCE.PD1-3", "NYCE.TRAC-5", "NYCE.TRBC-4"),
    "UninfectedControl" = c("NYCE.PD1-3", "NYCE.TRAC-5", "NYCE.TRBC-4"),
    "NoTemplateControl" = c("NYCE.PD1-3", "NYCE.TRAC-5", "NYCE.TRBC-4")
)

cat("Sample Treatment Table:")
pandoc.table(
  t(as.data.frame(treated_with)), 
  style = "simple", 
  split.table = Inf,
  justify = "left",
  emphasize.rownames = FALSE)

## log table for on target sites
## log table for on target sites
on_target_sites <- c(
  "NYCE:PD1" = "chr2:-:241858808",
  "NYCE:TRAC" = "chr14:+:22547664",
  "NYCE:TRBC1" = "chr7:+:142792020",
  "NYCE:TRBC2" = "chr7:+:142801367")

cat("On-target CRISPR Editing Genomic Locations:")
pandoc.table(
  data.frame(
    "Target_Gene" = names(on_target_sites),
    "Genomic_Locus" = on_target_sites,
    row.names = NULL),
  style = "simple",
  split.table = Inf,
  justify = "left",
  emphasize.rownames = FALSE)

guideRNASeqs <- unlist(DNAStringSetList(lapply(guideRNASeqs, DNAStringSet)))

upstream_dist <- 1200L
downstream_dist <- 30L

runDate <- "170907"
hg38 <- BSgenome.Hsapiens.UCSC.hg38
data_dir <- "~/data/projects/guideseq_analysis/170907_ENG02_clinical_run/processedData"
algnmt_file <- file.path(data_dir, "unique_sites.170907.csv")
code_dir <- dirname(
  sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(code_dir, "postProcessingSupport.R"))
submat <- banmat()
ref_dir <- "~/data/util_files"
ref_genes <- readRDS(file.path(ref_dir, "hg38.refSeq.rds"))
ref_exons <- extract_exons(
  ref_genes, ref_genes$exonStarts, ref_genes$exonEnds, ref_genes$name2)
onco_genes_data <- read.delim(
  file.path(ref_dir, "allOnco.human.v3.tsv"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE)
onco_genes <- unique(onco_genes_data[,"symbol"])
bad_actors <- read.delim(
  file.path(ref_dir, "humanLymph.v1.list"),
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE)[,1]
oregs <- readRDS(file.path(ref_dir, "oreganno.hg38.2017-06-09.rds"))

# Load input data --------------------------------------------------------------
#' Load all read alignments 
reads <- read.csv(algnmt_file)

# Print out stats during analysis.
cat("Tabulation of aligned templates per specimen:")
temp_table <- table(str_extract(reads$sampleName, "[\\w]+"))
pandoc.table(
  data.frame(
    "Specimen" = names(temp_table), 
    "Aligned_Templates" = as.numeric(temp_table),
    row.names = NULL),
  style = "simple",
  justify = "left",
  emphasize.rownames = FALSE)
rm(temp_table)

# Process input data -----------------------------------------------------------
# Format input alignments ======================================================
algnmts <- mutate(reads, specimen = str_extract(sampleName, "[\\w]+")) %>%
  mutate(seqOrt = ifelse(grepl("pos", sampleName), "pos", "neg")) %>%
  group_by(seqnames, start, end, strand, specimen, seqOrt) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  select(seqnames, start, end, strand, specimen, seqOrt, count) %>%
  as.data.frame()

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
    seqinfo = seqinfo(hg38))
mcols(algnmts_gr) <- select(
    algnmts, specimen, seqOrt, count)

# Analyze alignments -----------------------------------------------------------
algnmts_gr <- get_site_seqs(
  algnmts_gr, upstream_flank = upstream_dist, 
  downstream_flank = downstream_dist, ref_genome = hg38)

algnmts_gr <- compare_guideRNAs(
  algnmts_gr, guideRNASeqs = guideRNASeqs, 
  submat = submat, seq_col = "flanking.sequence",
  upstream_flank = upstream_dist, downstream_flank = downstream_dist)

sites <- filter_inappropriate_comparisons(
  algnmts_gr, spec_col = "specimen", treated_with)
probable_sites <- sites[sites$guideRNA.match != "No_valid_match",]
probable_sites$on_off_target <- ifelse(
  probable_sites$edit.site %in% on_target_sites, "On-target", "Off-target")
#probable_sites <- probable_sites[!grepl("N", probable_sites$aligned.sequence)]

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
crispr_sites$gene_id <- assign_gene_id(
  crispr_sites$site.chr, 
  positions = crispr_sites$site.pos, 
  reference = hg38, 
  ref_genes = ref_genes, 
  onco_genes = onco_genes, 
  bad_actors = bad_actors)
crispr_sites$gene_id <- sapply(1:nrow(crispr_sites), function(i, df, on){
  if(df$edit.site[i] %in% on){
    return(paste0(names(on)[match(df$edit.site[i], on)], "*"))
  }else{
    return(df$gene_id[i])
  }
}, df = crispr_sites, on = on_target_sites)

crispr_sites$on_off_target <- ifelse(
  crispr_sites$edit.site %in% on_target_sites, "On-target", "Off-target")

crispr_gr <- GRanges(
  seqnames = str_extract(crispr_sites$edit.site, "[\\w]+"),
  ranges = IRanges(
    start = as.numeric(str_extract(crispr_sites$edit.site, "[\\d]+$")),
    width = 1L),
  strand = str_extract(crispr_sites$edit.site, "[\\+\\-\\*]"))

crispr_gr <- doAnnotation(
  "within", crispr_gr, ref_genes, 
  colnam = "inGene", feature.colnam = "name2")
crispr_gr <- doAnnotation(
  "within", crispr_gr, ref_exons, 
  colnam = "inExon", feature.colnam = "geneNames")
crispr_gr <- doAnnotation(
  "nearest", crispr_gr, ref_genes, 
  colnam = "nearestGene", feature.colnam = "name2")
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
  "nearest", crispr_gr, ref_genes[ref_genes$name2 %in% onco_genes],
  colnam = "nearestOnco", feature.colnam = "name2")

crispr_sites <- bind_cols(
  crispr_sites, 
  as.data.frame(mcols(crispr_gr)[,c(
    "inGene", "inExon", "nearestGene", "nearestGeneDist", "inReg", "inRegType",
    "nearestReg", "nearestRegDist", "nearestRegType", "inTFbs", "nearestTFbs",
    "nearestTFbsDist", "nearestOnco", "nearestOncoDist")]))


# Write output data ------------------------------------------------------------
saveRDS(
  list(
    "edit_sites" = sites, 
    "mated_edit_sites" = mated_edit_sites, 
    "edit_loci_plots" = edit_loci_plots, 
    "cripsr_sites" = crispr_sites), 
  file = file.path(data_dir, paste0("edited_sites.", runDate, ".rds")))


q()

# One-off analyses -------------------------------------------------------------
#' Summary of information
summary <- data.frame(
  "Sample_Num" = 1:length(reduced_ranges_df),
  "SampleName" = sapply(reduced_ranges_df, function(x) unique(x$sampleName)),
  "Cas9mRNA" = c("Yes", "Yes", "No", "No", "Yes", "Yes", "No"),
  "guideRNA" = c("Yes (NG)", "Yes (NYCE)", "Yes (NYCE)", "No", "No", "Yes (NYCE)", "No"),
  "dsODN" = c("Yes", "Yes", "Yes", "Yes", "Yes", "No", "No"),
  "Total_reads" = sapply(reduced_ranges_df, nrow),
  "pct_PD1" = sapply(
    reduced_ranges_df, function(x){
      round(100 * length(grep("chr2_PD1", x$rname))/nrow(x), digits = 1)}),
  "pct_TRBC1" = sapply(
    reduced_ranges_df, function(x){
      round(100 * length(grep("chr7_TRBC1", x$rname))/nrow(x), digits = 1)}),
  "pct_TRBC2" = sapply(
    reduced_ranges_df, function(x){
      round(100 * length(grep("chr7_TRBC2", x$rname))/nrow(x), digits = 1)}),
  "pct_TRAC" = sapply(
    reduced_ranges_df, function(x){
      round(100 * length(grep("chr14_TRAC", x$rname))/nrow(x), digits = 1)}),
  row.names = NULL)

pandoc.table(summary, split.tables = Inf)
write.csv(
  summary, 
  file = "output_data_and_graphics/initial_samples/170411_on-target_summary_for_initial_samples.csv", 
  row.names = FALSE, 
  quote = TRUE)

summary <- select(summary, -SampleName) %>% as.data.frame()



knit(input = "output_data_and_graphics/initial_samples/170411_initial_on-target_summary_report.Rmd", 
     output = "output_data_and_graphics/initial_samples/170411_initial_on-target_summary_report.md")
markdownToHTML(file = "output_data_and_graphics/initial_samples/170411_initial_on-target_summary_report.md", 
               output = "output_data_and_graphics/initial_samples/170411_initial_on-target_summary_report.html", 
               extensions = c('tables'))
system("source activate guideseq && wkhtmltopdf 170411_initial_on_target_summary_report.html 170411_initial_on_target_summary_report.pdf")


render(input = "output_data_and_graphics/initial_samples/170411_initial_on-target_summary_report.Rmd", 
       ouput_format = )


lapply(flanking_templates, function(x) lapply(x, function(y) flank(reduce(y, min.gapwidth = 25L), -1)))

bounds <- lapply(flanking_templates, function(x) lapply(x, function(y) flank(reduce(y, min.gapwidth = 25L), -1)))
b <- lapply(bounds, function(x) lapply(x, function(y) GRanges(seqnames = unique(seqnames(y)), ranges = IRanges(start = min(start(y)), end = max(start(y))), strand = "+")))
lapply(guideRNASeqs, function(x) filter(data.frame(pct.id = round(pairwiseAlignment(getSeq(hg38, b), x, type = "overlap", substitutionMatrix = submat, scoreOnly = TRUE) / width(b) * 100), width = width(b), names = names(b)), width >= 5, pct.id >= 80))
