# Supporting functions for postProcessing GuideSeq data
get_site_seqs <- function(gr, upstream_flank, downstream_flank, ref_genome){
  require(BSgenome)
  require(GenomicRanges)
  fl_starts <- flank(gr, width = -1, start = TRUE)
  seq_ranges <- trim(flank(
    shift(
      fl_starts, 
      shift = ifelse(
        strand(fl_starts) == "+", 
        downstream_flank, 
        -downstream_flank)),
    width = (upstream_flank + downstream_flank),
    start = TRUE))
  gr$flanking.sequence <- getSeq(ref_genome, seq_ranges)
  gr
}

get_edit_site_and_seqs <- function(grl_of_clusters, flanking_seq, ref_genome){
  require(BSgenome)
  require(GenomicRanges)
  sites <- unlist(GRangesList(lapply(grl_of_clusters, function(gr, flanking_seq){
    fl.red <- flank(reduce(gr), width = -1, start = TRUE)
    GRanges(
      seqname = unique(seqnames(gr)),
      ranges = IRanges(
        start = min(start(fl.red)) - flanking_seq,
        end = max(start(fl.red)) + flanking_seq),
      strand = "+",
      sampleName = unique(gr$sampleName),
      count = sum(gr$count))
  }, flanking_seq = flanking_seq)))
  sites$sequence <- getSeq(ref_genome, sites)
  sites
}

aln_guideRNAs <- function(seqs, guideRNASeqs, tolerance){
  require(IRanges)
  require(Biostrings)
  require(dplyr)
  
  if(is.null(names(seqs))) names(seqs) <- 1:length(seqs)
  
  nt_widths <- data.frame(
    names = names(seqs),
    nt_width = width(seqs))
  
  alns <- lapply(0:tolerance, function(tol, guides, seqs){
    lapply(guides, function(guide){
      unlist(vmatchPattern(guide, seqs, max.mismatch = tol, fixed = 'subject'))
    })
  }, guides = guideRNASeqs, seqs = seqs)
  
  df <- do.call(rbind, lapply(1:length(alns), function(i){
    do.call(rbind, lapply(1:length(alns[[i]]), function(j){
      d <- as.data.frame(alns[[i]][[j]])
      d$guideRNA <- names(alns[[i]][j])
      d$mismatches <- i-1
      d
    }))
  }))
  
  group_by(df, start, end, width, names, guideRNA) %>%
    mutate(guideRNA.mismatch = min(mismatches)) %>%
    ungroup() %>%
    select(names, guideRNA, guideRNA.mismatch, start, end, width) %>%
    distinct() %>%
    left_join(., nt_widths, by = "names") %>%
    mutate(start = ifelse(start <= 0, 1, start)) %>%
    mutate(end = ifelse(
      end > nt_width, nt_width, end)) %>%
    mutate(width = end - start + 1)
}

calc_cutSite <- function(sites, matched_seqs, upstream_flank, 
                         downstream_flank, PAM, offset_nt){
  df <- data.frame(
    "chr" = seqnames(sites),
    "strand" = strand(sites),
    "pos" = start(flank(sites, -1, start = TRUE)),
    "guideRNA" = matched_seqs$guideRNA,
    "guide.aln" = ifelse(!grepl("(rev)", matched_seqs$guideRNA), "+", "-"),
    "guide.start" = matched_seqs$start,
    "guide.end" = matched_seqs$end)
  
  df$true.ort <- ifelse(df$strand == df$guide.aln, "+", "-")
  
  df$edit.pos <- ifelse(
    df$strand == "+",
    ifelse(
      df$guide.aln == "+",
      df$pos - upstream_flank + df$guide.end - nchar(PAM) - offset_nt,
      df$pos - df$guide.end + downstream_flank + nchar(PAM) + offset_nt - 1),
    ifelse(
      df$guide.aln == "+",
      df$pos + upstream_flank - df$guide.end + nchar(PAM) + offset_nt,
      df$pos + df$guide.end - downstream_flank - nchar(PAM) - offset_nt + 1
    ))
  
  return(paste0(df$chr, ":", df$true.ort, ":", df$edit.pos))
}

compare_guideRNAs <- function(gr_with_sequences, seq_col, 
                              guideRNASeqs, submat = NULL, tolerance = 6L,
                              upstream_flank, downstream_flank,
                              PAM = "NGG", offset_nt = 4L){ #23 matching all, -6 mismatches in guide
  require(Biostrings)
  require(GenomicRanges)
  
  if(is.null(submat)) submat <- banmat()
  pam <- paste0(gsub("N", ".", PAM), "$")
  
  sites <- gr_with_sequences
  sites$siteID <- 1:length(sites)
  seq_col_match <- match(seq_col, names(mcols(sites)))
  if(length(seq_col_match) == 0){
    stop("Cannot find sequences in column.")
  }else{
    seqs <- mcols(sites)[,seq_col_match]
    names(seqs) <- sites$siteID
  }
  
  fwd_df <- aln_guideRNAs(seqs, guideRNASeqs, tolerance)
  
  rev_df <- aln_guideRNAs(reverseComplement(seqs), guideRNASeqs, tolerance)
  
  fwd_df$aln.seq <- as.character(DNAStringSet(
    seqs[fwd_df$names], 
    start = fwd_df$start, 
    end = fwd_df$end))
  rev_df$aln.seq <- as.character(DNAStringSet(
    reverseComplement(seqs[rev_df$names]), 
    start = rev_df$start, 
    end = rev_df$end))
  rev_df$guideRNA <- paste0(rev_df$guideRNA, " (rev)")
  
  matched_seqs <- rbind(fwd_df, rev_df)
  matched_seqs <- matched_seqs[grepl(pam, matched_seqs$aln.seq),]
  
  good_alns <- as.numeric(matched_seqs$names)
  
  if(length(good_alns) == 0){
    sites$guideRNA.match <- "No_valid_match"
    sites$guideRNA.score <- NA
    sites$aligned.sequence <- NA
    sites$edit.site <- NA
    return(sites)
  }
  
  if(length(good_alns) != length(sites)){
    non_probable_sites <- sites[-good_alns]
    non_probable_sites$guideRNA.match <- "No_valid_match"
    non_probable_sites$guideRNA.mismatch <- NA
    non_probable_sites$aligned.sequence <- NA
    non_probable_sites$edit.site <- NA
  }
  
  potential_sites <- sites[good_alns]
  
  potential_sites$guideRNA.match <- matched_seqs$guideRNA
  potential_sites$guideRNA.mismatch <- matched_seqs$guideRNA.mismatch
  potential_sites$aligned.sequence <- matched_seqs$aln.seq
  
  potential_sites$edit.site <- calc_cutSite(
    potential_sites, matched_seqs, upstream_flank, 
    downstream_flank, PAM, offset_nt)
  
  if(length(good_alns) != length(sites)){
    all_sites <- c(potential_sites, non_probable_sites)
  }else{
    all_sites <- potential_sites
  }
  all_sites <- all_sites[order(all_sites$siteID)]
  if(any(duplicated(all_sites$siteID))){
    message("Some alignments with multiple guideRNA matches. Total number of alignments expanded. See siteID.")
  }else{  
    all_sites$siteID <- NULL
  }
  all_sites
}

calc_coverage <- function(gr, resolution){ ###!!!!!!!!Add option for counting coverage by reads or uniq frags
  #Set up coverage gr
  strandless <- gr
  strand(strandless) <- "*"
  gr_ranges <- range(strandless)
  
  window_seqs <- lapply(gr_ranges, function(chr, res){
    seq(start(chr), end(chr), res)
  }, res = resolution)
  
  coverage_grl <- GRangesList(lapply(
    1:length(gr_ranges), function(i, gr_ranges, window_seqs){
      seqname <- seqnames(gr_ranges[i])
      window <- window_seqs[[i]]
      GRanges(
        seqnames = rep(seqname, length(window)),
        ranges = IRanges(
          start = window, width = rep(resolution, length(window))),
        strand = rep("*", length(window)))
    }, gr_ranges = gr_ranges, window_seqs = window_seqs))
  
  coverage_pos <- coverage_grl
  coverage_pos <- GRangesList(lapply(coverage_pos, function(x){
    strand(x) <- rep("+", length(x))
    x}))
  coverage_neg <- coverage_grl
  coverage_neg <- GRangesList(lapply(coverage_pos, function(x){
    strand(x) <- rep("-", length(x))
    x}))
  
  bind_rows(lapply(1:length(coverage_grl), function(i, gr){
    as.data.frame(coverage_grl[[i]], row.names = NULL) %>%
      select(seqnames, start, end, width) %>%
      mutate(
        readCountsPos = countOverlaps(coverage_pos[[i]], gr),
        readCountsNeg = countOverlaps(coverage_neg[[i]], gr)) %>%
      arrange(seqnames)
  }, gr = gr))
}

plot_coverage <- function(gr, resolution = 10L){
  df <- calc_coverage(gr, resolution)
  ggplot(df, aes(x = start)) + 
    geom_bar(
      aes(y = readCountsPos), 
      stat = "identity", fill = "blue", width = resolution) +
    geom_bar(
      aes(y = -readCountsNeg), 
      stat = "identity", fill = "red", width = resolution) +
    facet_grid(. ~ seqnames, scales = "free") +
    geom_abline(slope = 0, intercept = 0, color = "grey") +
    labs(x = "Genomic Loci Position", y = "Read Counts") +
    theme_bw() +
    theme(
      legend.position = "None",
      strip.background = element_rect(fill = "white"),
      panel.border = element_rect(color = "white"),
      panel.grid.major = element_line(color = "white"),
      panel.grid.minor = element_line(color = "white"),
      axis.text = element_text(color = "black"),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"))
}

plot_edit_sites <- function(gr, sampleName = NULL, resolution = 10L){
  stopifnot(length(unique(gr$edit.site)) == 1)
  if(!is.null(sampleName)){  
    isThere <- match(sampleName, names(mcols(gr)))
    if(length(isThere) == 0) stop("SampleName column not found.")
    sample <- unique(as.character(mcols(gr)[,isThere]))
  }else{
    sample <- "Unspecified"
  }
  guide <- unique(str_extract(gr$guideRNA.match, "[\\w\\-\\.]+"))
  edit_site <- as.character(unique(gr$edit.site))
  edit_site <- unlist(strsplit(edit_site, ":"))
  edit_pos <- as.numeric(edit_site[3])
  edit_site[3] <- format(edit_pos, big.mark = ",", scientific = FALSE)
  df <- calc_coverage(gr, resolution)
  df$edit.site <- paste0(
    "Sample: ", sample, 
    "\nGuide: ", guide,
    "\nEdit Site: ", paste(edit_site, collapse = ":"))
  ggplot(df, aes(x = start)) + 
    geom_bar(
      aes(y = readCountsPos), 
      stat = "identity", fill = "blue", width = resolution) +
    geom_bar(
      aes(y = -readCountsNeg), 
      stat = "identity", fill = "red", width = resolution) +
    facet_grid(. ~ edit.site, scales = "free") +
    geom_abline(slope = 0, intercept = 0, color = "black") +
    geom_vline(xintercept = edit_pos, color = "black") +
    labs(x = "Genomic Loci Position", y = "Template Counts") +
    theme_bw() +
    theme(
      legend.position = "None",
      strip.background = element_rect(fill = NA, linetype = 0),
      strip.text = element_text(size = 14, lineheight = 1.1),
      panel.border = element_rect(color = "white"),
      plot.background = element_rect(color = "white"),
      panel.grid.major = element_line(color = "white"),
      panel.grid.minor = element_line(color = "white"),
      axis.title = element_text(color = "black", size = 14),
      axis.text = element_text(color = "black", size = 12),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"))
}

parse_ambiguous_cut_sites <- function(gr){
  require(GenomicRanges)
  require(stringr)
  guideRNA_matches <- strsplit(gr$guideRNA.match, ";")
  aligned_seqs <- strsplit(gr$aligned.sequence, ";")
  edit_sites <- strsplit(gr$edit.site, ";")
  stopifnot(all(sapply(guideRNA_matches, length) == sapply(aligned_seqs, length)))
  
  exp_gr <- gr[Rle(
    values = 1:length(gr), 
    lengths = sapply(guideRNA_matches, length))]
  
  exp_gr$guideRNA.match <- unlist(guideRNA_matches)
  exp_gr$aligned.sequence <- unlist(aligned_seqs)
  exp_gr$edit.site <- unlist(edit_sites)
  split(exp_gr, str_extract(exp_gr$guideRNA.match, "[\\w\\.\\-]+"))
}

filter_inappropriate_comparisons <- function(sites, spec_col, treated_with){
  require(GenomicRanges)
  require(stringr)
  
  spec_col_match <- match(spec_col, names(mcols(sites)))
  if(length(spec_col_match) == 0){
    stop("Cannot find specimen in column.")
  }else{
    specimens <- mcols(sites)[,spec_col_match]
  }
  
  if(!all(specimens %in% names(treated_with))){
    stop("Not all specimens accounted for in 'treated_with' list.")
  }
  
  sites$ori.order <- 1:length(sites)
  
  grl <- split(sites, specimens)
  treated_with <- treated_with[unique(sites$specimen)]
  grl <- grl[names(treated_with)]
  
  edited_sites <- unname(unlist(GRangesList(mapply(function(gr, treated){
    gr$guideRNA.match <- ifelse(
      str_extract(gr$guideRNA.match, "[\\w\\-\\.]+") %in% treated, 
      gr$guideRNA.match, 
      "No_valid_match")
    gr
  }, gr = grl, treated = treated_with, SIMPLIFY = FALSE))))
  
  edited_sites[order(edited_sites$ori.order)]
  edited_sites$ori.order <- NULL
  edited_sites
}

assign_gene_id <- function(seqnames, positions, reference, ref_genes, 
                           onco_genes, bad_actors, annotations = TRUE){
  require(GenomicRanges)
  require(hiAnnotator)
  
  gr <- GRanges(
    seqnames = seqnames,
    ranges = IRanges(start = positions, width = 1),
    strand = "*",
    seqinfo = seqinfo(reference))
  
  # Annotate Sites with Gene names for within and nearest genes
  gr <- getSitesInFeature(
    gr, ref_genes, colnam = "in_gene", feature.colnam = "name2")
  gr <- getNearestFeature(
    gr, ref_genes, colnam = "nearest_gene", feature.colnam = "name2")
  
  ## Add gene marks ("*" for in_gene, "~" for onco gene, and "!" for bad_actors)
  gr$gene_id_wo_annot <- ifelse(
    gr$in_gene == "FALSE", gr$nearest_gene, gr$in_gene)
  gr$gene_id_wo_annot <- sapply(strsplit(gr$gene_id_wo_annot, ","), "[[", 1)
  
  gr$gene_id <- ifelse(
    gr$in_gene == "FALSE", gr$gene_id_wo_annot, paste0(gr$gene_id_wo_annot, "*"))
  
  gr$gene_id <- ifelse(
    gr$gene_id_wo_annot %in% onco_genes, paste0(gr$gene_id, "~"), gr$gene_id)
  
  gr$gene_id <- ifelse(
    gr$gene_id_wo_annot %in% bad_actors, paste0(gr$gene_id, "!"), gr$gene_id)
  
  if(annotations){
    return(gr$gene_id)
  }else{
    return(gr$gene_id_wo_annot)
  }
}

refine_edit_sites <- function(gr, seq_col, guideRNASeqs, offset_nt = 4L, 
                              PAM = "NGG", submat = NULL){
  require(BSgenome)
  require(GenomicRanges)
  require(Biostrings)
  require(stringr)
  
  if(is.null(submat)) submat <- banmat()
  
  seqs <- mcols(gr)[,match(seq_col, names(mcols(gr)))]
  revl <- grepl("(rev)", gr$guideRNA.match)
  
  alns <- mapply(function(seq, rl, guide, guideRNASeqs, submat){
    guide <- str_extract(guide, "[\\w\\-\\.]+")
    if(!rl){
      guideRNA <- guideRNASeqs[[guide]]
    }else{
      guideRNA <- reverseComplement(guideRNASeqs[[guide]])
    }
    pairwiseAlignment(
      seq, guideRNA, type = "overlap", substitutionMatrix = submat)},
    seq = seqs,
    rl = revl,
    guide = gr$guideRNA.match,
    MoreArgs = list(guideRNASeqs = guideRNASeqs, submat = submat),
    SIMPLIFY = FALSE)
  
  stopifnot(all(sapply(alns, score) == gr$guideRNA.score))
  
  gr$edit.sites <- sapply(
    1:length(gr), function(i, alns, gr, revl, PAM, offset_nt){
      aln <- alns[[i]]
      g <- gr[i]
      rl <- revl[i]
      if(as.logical(strand(g) == "+" & !rl)){
        pos <- start(g) + end(aln@pattern@range) - nchar(PAM) - offset_nt
        return(paste0(seqnames(g), ":+:", pos))
      }else if(as.logical(strand(g) == "+" & rl)){
        pos <- start(g) + start(aln@pattern@range) + nchar(PAM) + offset_nt
        return(paste0(seqnames(g), ":-:", pos))
      }else if(as.logical(strand(g) == "-" & rl)){
        pos <- end(g) - end(aln@pattern@range) + nchar(PAM) + offset_nt
        return(paste0(seqnames(g), ":-:", pos))
      }else if(as.logical(strand(g) == "-" & !rl)){
        pos <- end(g) - start(aln@pattern@range) - nchar(PAM) - offset_nt
        return(paste0(seqnames(g), ":+:", pos))
      }},
    alns = alns, gr = gr, revl = revl, PAM = PAM, offset_nt = offset_nt)
  gr
}

extract_exons <- function(refGenes, exonStarts, exonEnds, geneNames){
  stopifnot(
    length(refGenes) == length(exonStarts) &
    length(exonStarts) == length(exonEnds) &
    length(exonEnds) == length(geneNames))
  
  exonStarts <- gsub(",$", "", exonStarts)
  exonEnds <- gsub(",$", "", exonEnds)
  
  exonStarts <- strsplit(exonStarts, ",")
  exonEnds <- strsplit(exonEnds, ",")
  stopifnot(all(sapply(
    1:length(geneNames), 
    function(i) length(exonStarts[[i]]) == length(exonEnds[[i]]))))
  
  gr <- GRanges(
    seqnames = Rle(
      values = as.character(seqnames(refGenes)), 
      lengths = sapply(exonStarts, length)),
    ranges = IRanges(
      start = as.numeric(unlist(exonStarts)),
      end = as.numeric(unlist(exonEnds))),
    strand = Rle(
      values = as.character(strand(refGenes)),
      length = sapply(exonStarts, length)),
    seqinfo = seqinfo(refGenes))
  
  geneNames <- as.character(Rle(
    values = geneNames,
    lengths = sapply(exonStarts, length)))
  exonNums <- unlist(mapply(function(strand, exSt){
    if(strand == "+"){
      return(seq_along(exSt))
    }else{
      return(rev(seq_along(exSt)))
    }},
    strand = as.character(strand(refGenes)),
    exSt = exonStarts))
  
  geneNames <- paste0(geneNames, ":exon", sprintf("%03d", exonNums))
  
  gr$geneNames <- geneNames
  gr
}

#' A Binary Ambiguous Nucleotide scoring Matrix (BAN Mat)
#' 
#' Constructed based on NUC4.4.
#' 
#' Meant for comparing ambiguous sequences against "A", "T", "G", "C", and "N"
#' containing sequences. Currently matches between ambiuous nucleotides are 
#' considered mismatch.
#' 
#' @author Christopher Nobles, Ph.D.

banmat <- function(){
  matrix(c(
    1,0,0,0,0,1,1,0,0,1,0,1,1,1,1,
    0,1,0,0,0,1,0,1,1,0,1,0,1,1,1,
    0,0,1,0,1,0,1,0,1,0,1,1,0,1,1,
    0,0,0,1,1,0,0,1,0,1,1,1,1,0,1,
    0,0,1,1,1,0,0,0,0,0,0,0,0,0,1,
    1,1,0,0,0,1,0,0,0,0,0,0,0,0,1,
    1,0,1,0,0,0,1,0,0,0,0,0,0,0,1,
    0,1,0,1,0,0,0,1,0,0,0,0,0,0,1,
    0,1,1,0,0,0,0,0,1,0,0,0,0,0,1,
    1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,
    0,1,1,1,0,0,0,0,0,0,1,0,0,0,1,
    1,0,1,1,0,0,0,0,0,0,0,1,0,0,1,
    1,1,0,1,0,0,0,0,0,0,0,0,1,0,1,
    1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
    ncol = 15,
    nrow = 15,
    byrow = TRUE,
    dimnames = list(
      c("A", "T", "G", "C", "S", "W", "R", "Y", 
        "K", "M", "B", "V", "H", "D", "N"),
      c("A", "T", "G", "C", "S", "W", "R", "Y", 
        "K", "M", "B", "V", "H", "D", "N")))
}
