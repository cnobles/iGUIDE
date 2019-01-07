#' For anyone reviewing the code below, the following is a small style guide 
#' outlining the various formats for the code. 
#' 
#' Names with "_": objects, inlucding data.frames, GRanges, vectors, ...
#' Names in caMel format: functions or components of objects (i.e. columns 
#' within a data.frame).
#' Names with ".": arguments / options for functions

# Set Global options and load intiial packages ---------------------------------
options(stringsAsFactors = FALSE, scipen = 99)
suppressMessages(library("argparse"))
suppressMessages(library("pander"))
panderOptions("table.style", "simple")
panderOptions("table.split.table", Inf)

# Set up and gather command line arguments -------------------------------------
parser <- ArgumentParser(
  description = "Simulate incorporation site data to generate test data for iGUIDE software.")
parser$add_argument(
  "config", nargs = 1, type = "character", help = "Configuration file.")
parser$add_argument(
  "-o", "--outfolder", nargs = 1, type = "character", default = "Data",
  help = "Output folder. [Data] by default.")
parser$add_argument(
  "-s", "--seed", nargs = 1, type = "integer", default = 1,
  help = "Set random seed, 1 by default.")

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
set.seed(args$seed)

input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(1:length(args), function(i){
    paste(args[[i]], collapse = ", ")}))
input_table <- input_table[
  match(
    c("config :", "outfolder :", "seed :"), 
    input_table$Variables),]
pandoc.title("SimIncorpData Inputs")
cat("\n[", paste(Sys.time()), "]")
pandoc.table(data.frame(input_table, row.names = NULL), 
             justify = c("left", "left"), 
             split.tables = Inf)

## Create output directory if not currently available --------------------------
if(!file.exists(args$outfolder)){
  attempt <- try(system(paste0("mkdir ", args$outfolder)))
  if(attempt == 1) stop("Cannot create output folder.")
}

# Load additional packages -----------------------------------------------------
packs <- c(
  "ShortRead", "Biostrings", "GenomicRanges", "BSgenome", 
  "tidyverse", "magrittr", "yaml")
packs_loaded <- suppressMessages(
  sapply(packs, require, character.only = TRUE))
if(!all(packs_loaded)){
  pandoc.table(data.frame(
    "R-Packages" = names(packs_loaded),
    "Loaded" = packs_loaded, 
    row.names = NULL))
  stop("Check dependancies.")
}

# Supporting functions ---------------------------------------------------------
offTargetLoci <- function(gRNA, pam, ref, num, mismatch, 
                            cut_pos = -4, output = "ids"){
  # Checks
  stopifnot(output %in% c("ids", "granges"))
  
  # Standardize input sequences
  pam <- toupper(pam)
  gRNA <- stringr::str_remove(toupper(gRNA), paste0(pam, "$"))
  
  # Identify reference location where gRNA matches
  algns <- Biostrings::vmatchPattern(gRNA, ref, max.mismatch = mismatch)
  # Verify which locatoins contain adjacent pam sequence matches
  pam_flanks <- GenomicRanges::flank(algns, 3, start = FALSE)
  pam_seqs <- BSgenome::getSeq(ref_genome, pam_flanks)
  pam_pos_idx <- which(stringr::str_detect(
    as.character(pam_seqs), gsub("N", ".", pam)))
  algns_pos <- algns[pam_pos_idx]
  
  # Append additional information to ranges
  algns_pos$seq <- BSgenome::getSeq(ref_genome, algns_pos)
  algns_pos$mis <- nchar(gRNA) - Biostrings::pairwiseAlignment(
    algns_pos$seq, gRNA, 
    substitutionMatrix = Biostrings::nucleotideSubstitutionMatrix(
      baseOnly = TRUE), 
    scoreOnly = TRUE)
  
  # Select top 'num' of sites and identify the edit location
  algns_pos <- algns_pos[order(algns_pos$mis)]
  algns_pos <- algns_pos[algns_pos$mis != 0]
  top_pos <- head(algns_pos, n = num)
  top_off_sites <- GenomicRanges::flank(
    GenomicRanges::flank(top_pos, cut_pos, start = FALSE), -1, 
    start = TRUE)
  
  # Return data as either ids or the grange object
  if(output == "ids"){
    return(paste0(
      GenomicRanges::seqnames(top_off_sites), ":", 
      GenomicRanges::strand(top_off_sites), ":", 
      GenomicRanges::start(top_off_sites)))
  }else if(output == "granges"){
    return(top_off_sites)
  }
}

selectRandomSites <- function(num, ref.genome, drop.extra.seqs = TRUE, 
                              seq.names = NULL, set.seed = NULL){
  if(!is.null(set.seed)) set.seed(set.seed)
  if(is.null(GenomicRanges::seqinfo(ref.genome))) stop("Ref genome does not have seqinfo.")
  
  if(!is.null(seq.names)){
    seq_lengths <- ref.genome@seqinfo@seqlengths[
      match(seq.names, GenomicRanges::seqnames(ref.genome))]
    names(seq_lengths) <- seq.names
  }else{
    if(drop.extra.seqs){
      seq.names <- grep("chr[0-9XYM]+$", GenomicRanges::seqnames(ref.genome), value = TRUE)
      seq_lengths <- ref.genome@seqinfo@seqlengths[
        match(seq.names, GenomicRanges::seqnames(ref.genome))]
      names(seq_lengths) <- seq.names
    }else{
      seq.names <- GenomicRanges::seqnames(ref.genome)
      seq_lengths <- ref.genome@seqinfo@seqlengths
      names(seq_lengths) <- seq.names
    }
  }
  
  chrs <- sort(factor(
    sample(seq.names, num, replace = TRUE, prob = seq_lengths),
    levels = seq.names))
  strands <- sample(c("+", "-"), num, replace = TRUE)
  splt_chrs <- split(chrs, chrs)
  splt_chrs <- splt_chrs[sapply(splt_chrs, length) > 0]
  positions <- unname(unlist(sapply(
    splt_chrs, function(seq, seq_len){
      sample.int(
        n = seq_len[unique(as.character(seq))], size = length(seq), 
        replace = FALSE, useHash = FALSE)
    }, seq_len = seq_lengths)))
  
  gr <- GenomicRanges::GRanges(
    seqnames = chrs,
    ranges = IRanges::IRanges(start = positions, width = 1L),
    strand = strands,
    seqinfo = seqinfo(ref.genome))
  
  gr
}

genIncorpDist <- function(posid, num, sd = 41){
  chr <- stringr::str_extract(posid, "[\\w]+")
  pos <- as.numeric(stringr::str_extract(posid, "[0-9]+$"))
  dist <- round(rnorm(num, 0, sd = sd))
  strands <- ifelse(
    dist > 0, "+", 
    ifelse(dist < 0, "-", sample(c("+", "-"), num, replace = TRUE)))
  paste0(chr, ":", strands, ":", pos + dist)
}

fillRandomNucleotides <- function(seqs, code_map = Biostrings::IUPAC_CODE_MAP){
  code_map <- code_map[!names(code_map) %in% c("A", "T", "G", "C")]
  mod_seqs <- new.env()
  mod_seqs$seqs <- seqs
  null <- lapply(seq_along(code_map), function(i){
    anuc <- names(code_map[i])
    nucs <- unlist(strsplit(code_map[i], ""))
    null <- lapply(seq_along(mod_seqs$seqs), function(j){
      cnt <- str_count(mod_seqs$seqs[j], anuc)
      ops <- sample(nucs, cnt, replace = TRUE)
      null <- lapply(seq_len(cnt), function(k){
        mod_seqs$seqs[j] <- str_replace(mod_seqs$seqs[j], anuc, ops[k])
      })
    })
  })
  return(mod_seqs$seqs)
}

captureReadSequence <- function(templates, primer, cycles, RC = FALSE){
  if(RC) templates <- reverseComplement(templates)
  primer_anneal <- vmatchPattern(primer, templates, max.mismatch = 0)
  DNAStringSet(
    templates, start = unlist(primer_anneal@ends) + 1, width = cycles)
}

# Load configuration -----------------------------------------------------------
config <- yaml.load_file(args$config)
pandoc.title("Configuration")
print.simple.list(config)

# Load sampleInfo
sample_info <- read.csv(config$sample_info)
pandoc.title("Sample Information")
pander(sample_info)

# Load ref_genome --------------------------------------------------------------
genome <- grep(
  config$ref_genome, unique(BSgenome::installed.genomes()), value = TRUE)
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
cat("\n[", paste(Sys.time()), "] Loaded reference genome: ", genome)

# Identify off-target sites ----------------------------------------------------
cat(
  "\n[", paste(Sys.time()), 
  "] Identifying top potential off-target sites for gRNA(s): ", 
  paste(names(config$guide_rna_sequences), collapse = ", "))

off_target_sites <- lapply(
  config$guide_rna_sequences, offTargetLoci, 
  pam = config$pam_sequence, 
  ref = ref_genome, 
  num = config$off_target_params$num, 
  mismatch = config$off_target_params$mismatch)

# Generate random sites in a uniform distribution across the reference ---------
cat(
  "\n[", paste(Sys.time()), 
  "] Generating random pool of ", 
  nrow(sample_info) * config$random_background$num, 
  " sites.")
if(!is.null(config$random_background$chrs)){
  random_gr <- selectRandomSites(
    num = nrow(sample_info) * config$random_background$num, 
    ref.genome = ref_genome, 
    drop.extra.seqs = TRUE, 
    seq.names = config$random_background$chrs,
    set.seed = args$seed)
}else{
  random_gr <- selectRandomSites(
    num = nrow(sample_info) * config$random_background$num, 
    ref.genome = ref_genome, 
    drop.extra.seqs = TRUE, 
    set.seed = args$seed)
}

random_sites <- paste0(
  seqnames(random_gr), ":", strand(random_gr), ":", start(random_gr))

# Consolidate on / off / false / random sites together -------------------------
on_target_df <- data.frame(
  class = "on",
  gRNA = stringr::str_extract(names(config$on_target_sites), "[\\w]+"),
  posid = unlist(config$on_target_sites),
  abund = config$on_target_params$abund,
  row.names = NULL)

off_target_df <- data.frame(
  class = "off",
  gRNA = rep(names(off_target_sites), lengths(off_target_sites)),
  posid = unlist(off_target_sites),
  abund = config$off_target_params$abund[
    do.call(c, lapply(off_target_sites, seq_along))],
  row.names = NULL)

false_target_df <- data.frame(
  class = "false",
  gRNA = "all",
  posid = unlist(config$false_target_sites),
  abund = config$false_target_params$abund,
  row.names = NULL)

random_target_df <- data.frame(
  class = "random",
  gRNA = "all",
  posid = random_sites,
  abund = sample(
    seq_len(config$random_background$max_abund), 
    length(random_sites), 
    replace = TRUE),
  row.names = NULL)

sites_df <- bind_rows(
  on_target_df, off_target_df, false_target_df, random_target_df)

pandoc.title("Sample from target data set.")
pander(rbind(head(sites_df, n = 5), tail(sites_df, n = 5)), row.names = FALSE)

# Initialize specimen specific grouping and determine ranges for seqs ----------
cat(
  "\n[", paste(Sys.time()), 
  "] Generating incorporation positions from targets.")
spec_sites <- lapply(
  structure(sample_info$gRNA, names = sample_info$sampleName),
  function(sgRNA){
    sgRNA <- unlist(strsplit(sgRNA, ";"))
    on_sites <- filter(on_target_df, gRNA %in% sgRNA)
    off_sites <- filter(off_target_df, gRNA %in% sgRNA)
    ran_sites <- random_target_df[
      sample(
        seq_len(nrow(random_target_df)), 
        round(1.1 * config$random_background$num), 
        replace = FALSE),]
    sites <- bind_rows(on_sites, off_sites, false_target_df, ran_sites)
    incorp_inter <- sites %>% 
      mutate(site = seq_len(n())) %>%
      group_by(site) %>%
      mutate(incorp = list(genIncorpDist(posid, abund))) %>%
      ungroup() %>%
      as.data.frame()
    
    breakpoint_dist <- round(rnorm(10000, 70, 250))
    breakpoint_dist <- breakpoint_dist[
      breakpoint_dist >= config$min_length & breakpoint_dist <= 2500]
    
    incorp <- incorp_inter[
      rep(seq_len(nrow(incorp_inter)), lengths(incorp_inter$incorp)),] %>%
      mutate(incorp = unlist(incorp_inter$incorp)) %>%
      select(class, gRNA, posid, incorp) %>%
      group_by(class, gRNA) %>%
      mutate(
        pos = as.numeric(stringr::str_extract(incorp, "[0-9]+$")),
        len = sample(breakpoint_dist, n(), replace = TRUE),
        brk = ifelse(
          str_detect(incorp, fixed("+")), pos + len - 1, pos - len + 1)) %>%
      ungroup()
    
    incorp_gr <- GRanges(
      seqnames = stringr::str_extract(incorp$incorp, "[\\w]+"),
      strand = stringr::str_extract(incorp$incorp, "[+-]"),
      ranges = IRanges(
        start = ifelse(
          stringr::str_detect(incorp$incorp, fixed("+")), 
          incorp$pos, incorp$brk),
        end = ifelse(
          stringr::str_detect(incorp$incorp, fixed("+")), 
          incorp$brk, incorp$pos)),
      seqinfo = seqinfo(ref_genome))
    mcols(incorp_gr) <- select(incorp, class, gRNA, posid, incorp)
    incorp_gr
  })

sample_idx <- sample(seq_along(spec_sites[[1]]), 10, replace = FALSE)
pandoc.title("Sample from incorporation data set.")
pander(as.data.frame(mcols(spec_sites[[1]][sample_idx]), row.names = NULL))

# Get sequences for ranges identified flanking the incorporation sites ---------
cat(
  "\n[", paste(Sys.time()), 
  "] Retrieving genomic sequences for template construction.")
spec_seqs <- lapply(spec_sites, function(set) getSeq(ref_genome, set))

# Build template ---------------------------------------------------------------
cat(
  "\n[", paste(Sys.time()), 
  "] Building templates with accessory sequences.")

incorp_side_seq <- structure(paste0(
  config$seqs$P7, 
  as.character(reverseComplement(DNAStringSet(sample_info$barcode1))), 
  config$seqs$PI1, 
  config$seqs$ODNbit),
  names = sample_info$sampleName)

breakpoint_side_seq <- structure(paste0(
  config$seqs$P5, sample_info$barcode2, config$seqs$PI2), 
  names = sample_info$sampleName)

template_seqs <- structure(lapply(
  seq_along(spec_seqs), 
  function(i, spec_seqs, incorp_side_seq, breakpoint_side_seq){
    incorp_side <- incorp_side_seq[i]
    breakpoint_side <- breakpoint_side_seq[i]
    seqs <- spec_seqs[[i]]
    gr <- spec_sites[[i]]
    n_idx <- which(stringr::str_detect(as.character(seqs), "N"))
    ran_idx <- which(gr$class == "random")
    remove_idx <- n_idx[which(n_idx %in% ran_idx)]
    if(length(remove_idx) > 0){
      seqs <- seqs[-remove_idx]
      gr <- gr[-remove_idx]
    }
  
    # Reduce number of random sites to specified count
    random_posids <- sample(
      unique(gr$posid[gr$class == "random"]), 
      size = config$random_background$num, 
      replace = FALSE)
    keep_idx <- which(gr$class != "random" | gr$posid %in% random_posids)
    seqs <- seqs[keep_idx]
    gr <- gr[keep_idx]
  
    # Append incorporation side sequence
    i_seqs <- Biostrings::DNAStringSet(paste0(incorp_side, as.character(seqs)))
    
    # Flip and append breakpoint side sequence
    rc_seqs_i <- Biostrings::reverseComplement(i_seqs)
    breakpoint_side_filled <- fillRandomNucleotides(
      rep(breakpoint_side, length(i_seqs)))
    b_seqs_i <- Biostrings::DNAStringSet(
      paste0(breakpoint_side_filled, as.character(rc_seqs_i)))
  
    # Append poly-T to end to simulate sequencing into the flowcell
    A_b_seqs_i_T <- Biostrings::DNAStringSet(paste0(
      paste(rep("A", 50), collapse = ""),
      as.character(b_seqs_i),
      paste(rep("T", 50), collapse = "")))
  
    return(list(seqs = A_b_seqs_i_T, gr = gr))
    }, spec_seqs = spec_seqs, 
    incorp_side_seq = incorp_side_seq, 
    breakpoint_side_seq = breakpoint_side_seq), 
  names = names(spec_seqs))

# Identify read-specific sequences ---------------------------------------------
cat(
  "\n[", paste(Sys.time()), 
  "] Simulating ", sum(lengths(template_seqs)), 
  " reads for sequencing cycles: ", 
  paste(names(config$seq_dists), collapse = ", "))

R1_seqs <- unlist(DNAStringSetList(lapply(
  lapply(template_seqs, "[[", "seqs"), 
  captureReadSequence, 
  primer = config$seqs$R1, 
  cycles = config$seq_dists$R1)))

I1_seqs <- unlist(DNAStringSetList(lapply(
  lapply(template_seqs, "[[", "seqs"), 
  captureReadSequence, 
  primer = config$seqs$I1, 
  cycles = config$seq_dists$I1)))

I2_seqs <- unlist(DNAStringSetList(lapply(
  lapply(template_seqs, "[[", "seqs"), 
  captureReadSequence, 
  primer = config$seqs$I2, 
  cycles = config$seq_dists$I2,
  RC = config$instrument != "miseq")))

R2_seqs <- unlist(DNAStringSetList(lapply(
  lapply(template_seqs, "[[", "seqs"), 
  captureReadSequence, 
  primer = config$seqs$R2, 
  cycles = config$seq_dists$R2,
  RC = TRUE)))

output_gr <- unlist(GRangesList(lapply(template_seqs, "[[", "gr")))

cat(
  "\n[", paste(Sys.time()), 
  "] Constructed reads. Saving output to: ", args$outfolder)

# Format for output ------------------------------------------------------------
read_names <- paste0(
  names(output_gr), "-", output_gr$class, "-", output_gr$gRNA, "-", 
  seqnames(output_gr), ":", ifelse(strand(output_gr) == "+", "pos", "neg"), ":",
  ifelse(strand(output_gr) == "+", start(output_gr), end(output_gr)), "-", 
  width(output_gr), ":", seq_along(output_gr))
names(R1_seqs) <- names(I1_seqs) <- read_names
names(I2_seqs) <- names(R2_seqs) <- read_names
output_gr$sampleName <- names(output_gr)
output_gr$read.name <- read_names
output_df <- as.data.frame(output_gr, row.names = NULL)
output_bed <- select(output_df, seqnames, start, end, strand, sampleName) %>%
  mutate(type = paste0("sim-", sampleName), score = 500) %>%
  select(seqnames, start, end, type, score, strand)

# Write truth files
cat("\n[", paste(Sys.time()), "] Writing truth files.")

write.table(
  output_bed, 
  file = file.path(args$outfolder, "truth.bed"), 
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.csv(
  output_df, 
  file = file.path(args$outfolder, "truth.csv"), 
  quote = FALSE, row.names = FALSE)

# Write Seq files
cat("\n[", paste(Sys.time()), "] Writing fastq files for reads.")
null <- mapply(function(read_seqs, path, fill){
  
    reads <- ShortRead::ShortReadQ(
      sread = read_seqs, 
      quality = ShortRead::SFastqQuality(Biostrings::BStringSet(rep(
        paste0(rep(fill, unique(Biostrings::width(read_seqs))), collapse = ""), 
        length(read_seqs)))),
      id = Biostrings::BStringSet(names(read_seqs)))
  
    ShortRead::writeFastq(
      reads, file = file.path(args$outfolder, path), compress = TRUE)
  }, 
  read_seqs = list(R1_seqs, I1_seqs, I2_seqs, R2_seqs), 
  path = config$seq_files, 
  MoreArgs = list(fill = config$fill_qual))

cat("\n[", paste(Sys.time()), "] Script completed.\n")
q()
