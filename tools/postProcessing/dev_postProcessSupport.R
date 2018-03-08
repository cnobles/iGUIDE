## GUIDE-seq matching method
## Pile up alignments and find the beginning of the most frequent alignment
## (based on fragment lengths) for the group. Look +/- 25 bp around the site for
## identifying on/off-target sites with the guideRNA.
## This needs to be made into a function for standardizing locations, solve the 
## bias problem by doing it twice on the data, since it's so fast, it shouldn't
## matter the extra round. Further, this pile up method would be good for
## identifying genomic regions to scan and would reduce the query load, just 
## increase the range.
pileupCluster <- function(gr, grouping = NULL, 
                              maxgap = 0L, return = "full"){
  if(!is.null(grouping)){
    group_vec <- mcols(gr)[,match(grouping, names(mcols(gr)))]
  }else{
    group_vec <- rep(1, length(gr))
  }
  
  # Process each group individually to keep clustering separate
  gr$grouping <- group_vec
  gr$order <- 1:length(gr)
  grl <- split(gr, gr$grouping)
  gr <- unname(unlist(GRangesList(lapply(grl, function(g){
    g_pos <- g[strand(g) == "+"]
    g_pos <- groupPileups(g_pos, "+", maxgap = maxgap)
    g_neg <- g[strand(g) == "-"]
    g_neg <- groupPileups(g_neg, "-", maxgap = maxgap)
    g_ret <- c(g_pos, g_neg)
    g_ret <- g_ret[order(g_ret$order)]
    g_ret$clusID <- as.integer(factor(g_ret$clus.ori))
    g_ret
  }))))
  
  gr$clusID <- as.integer(factor(gr$clus.ori))
  gr <- gr[order(gr$order)]
  gr$order <- gr$grouping <- NULL
  
  if(return == "full"){
    return(list("gr" = gr, "clusID" = gr$clusID, "clusOri" = gr$clus.ori))
  }else if(return == "simple"){
    return(gr$clus.ori)
  }
}

groupPileups <- function(gr, strand, maxgap = maxgap){
  # Implement axial cluster structure rather than all vs all.
  red.gr <- reduce(gr, min.gapwidth = maxgap, with.revmap = TRUE)
  axil_nodes <- unlist(Rle(
    values = unlist(red.gr$revmap)[start(red.gr$revmap@partitioning)], 
    lengths = width(red.gr$revmap@partitioning)))
  nodes <- unlist(red.gr$revmap)
  edgelist <- unique(matrix( c(axil_nodes, nodes), ncol = 2 ))
  pile_ups <- igraph::clusters(graph.edgelist(edgelist, directed = FALSE))
  gr$clusID <- membership(pile_ups)
  if(strand == "+"){
    pile_starts <- sort(split(Rle(start(gr)), gr$clusID))
    pile_starts <- runValue(pile_starts)[
      runLength(pile_starts) == max(runLength(pile_starts))]
    pile_starts <- unlist(pile_starts)[start(pile_starts@partitioning)]
  }else{
    pile_starts <- sort(split(Rle(end(gr)), gr$clusID))
    pile_starts <- runValue(pile_starts)[
      runLength(pile_starts) == max(runLength(pile_starts))]
    pile_starts <- unlist(pile_starts)[end(pile_starts@partitioning)]
  }
  gr <- gr[order(gr$clusID)]
  gr$clus.ori <- paste0(
    seqnames(gr), ":", strand(gr), ":", unlist(
      Rle(values = pile_starts, lengths = pile_ups$csize)))
  gr
}  

identifyPairedAlgnmts <- function(gr, maxgap, grouping = NULL){
  if(!is.null(grouping)){
    group_vec <- mcols(gr)[,match(grouping, names(mcols(gr)))]
  }else{
    group_vec <- rep(1, length(gr))
  }
  
  gr$ori.order <- 1:length(gr)
  grl <- split(gr, group_vec)
  
  pairs <- unname(unlist(lapply(1:length(grl), function(i, grl, maxgap){
    gs <- grl[[i]]
    red <- reduce(
      flank(gs, width = -1, start = TRUE), 
      min.gapwidth = 0L, with.revmap = TRUE)
  
    hits <- findOverlaps(red, maxgap = maxgap, ignore.strand = TRUE)
    hits <- as.data.frame(hits, row.names = NULL) %>%
      mutate(
        qStart = start(red[queryHits]),
        sStart = start(red[subjectHits]),
        qStrand = as.character(strand(red[queryHits])),
        sStrand = as.character(strand(red[subjectHits]))) %>%
      filter(
        queryHits != subjectHits, 
        qStrand != sStrand,
        qStrand == "+",
        qStart >= sStart) 
    
    if(nrow(hits) == 0){
      return(rep(NA, length(gs)))
    }else{
      el <- hits[,c("queryHits", "subjectHits")]
      g <- simplify(graph_from_edgelist(as.matrix(el), directed = FALSE))
      mem <- data.frame(
          idx = 1:vcount(g), 
          grp = as.numeric(membership(clusters(g)))) %>%
        group_by(grp) %>%
        filter(n() > 1) %>%
        ungroup() %>%
        mutate(grp = as.integer(as.factor(grp))) %>%
        as.data.frame()
      mem$revmap <- red$revmap[mem$idx]
      mem_exp <- data.frame(
        gs_idx = unlist(mem$revmap),
        grp = paste0(
          names(grl)[i], ".", 
          Rle(values = mem$grp, lengths = lengths(mem$revmap))))
      ret <- rep(NA, length(gs))
      ret[mem_exp$gs_idx] <- mem_exp$grp
      return(ret)
    }
  }, grl = grl, maxgap = maxgap)))
  
  if(is.null(grouping)){
    pairs <- as.numeric(stringr::str_extract(pairs, "[0-9]+$"))}
  
  un_grl <- unlist(grl)
  pairs[order(un_grl$ori.order)]
}

# Supporting functions
# Format number in tables with big.marks conveinently
pNums <- function(x, ...){
  format(x, big.mark = ",", ...)
}

# Convert matches to GRanges objects
mindexToGranges <- function(mindex, strand, ref = NULL){
  if(is.null(mindex@NAMES)) stop("NAMES column not found for seqnames.")
  ir <- unlist(mindex)
  if(is.null(ref)){
    return(unname(GRanges(
      seqnames = names(ir),
      ranges = ir,
      strand = rep(strand, length(ir)))))
  }else{
    return(unname(GRanges(
      seqnames = names(ir),
      ranges = ir,
      strand = rep(strand, length(ir)),
      seqinfo = seqinfo(ref))))
  }
}

# Generate random sites across the reference genome.
selectRandomSites <- function(num, refGenome, drop_extra_seqs = TRUE, 
                              seqNames = NULL, setSeed = NULL){
  require(GenomicRanges)
  if(!is.null(setSeed)) set.seed(setSeed)
  if(is.null(seqinfo(refGenome))) stop("Ref genome does not have seqinfo.")
  
  if(!is.null(seqNames)){
    seqLengths <- refGenome@seqinfo@seqlengths[
      match(seqNames, seqnames(refGenome))]
    names(seqLengths) <- seqNames
  }else{
    if(drop_extra_seqs){
      seqNames <- grep("chr[0-9XY]+$", seqnames(refGenome), value = TRUE)
      seqLengths <- refGenome@seqinfo@seqlengths[
        match(seqNames, seqnames(refGenome))]
      names(seqLengths) <- seqNames
    }else{
      seqNames <- seqnames(refGenome)
      seqLengths <- refGenome@seqinfo@seqlengths
      names(seqLengths) <- seqNames
    }
  }
  
  chrs <- sort(factor(
    sample(seqNames, num, replace = TRUE, prob = seqLengths),
    levels = seqNames))
  strands <- sample(c("+", "-"), num, replace = TRUE)
  spltChrs <- split(chrs, chrs)
  spltChrs <- spltChrs[sapply(spltChrs, length) > 0]
  positions <- unname(unlist(sapply(
    spltChrs, function(seq, seqLen){
      sample.int(
        n = seqLen[unique(as.character(seq))], size = length(seq), 
        replace = FALSE, useHash = FALSE)
    }, seqLen = seqLengths)))
  
  gr <- GRanges(
    seqnames = chrs,
    ranges = IRanges(start = positions, width = 1L),
    strand = strands,
    seqinfo = seqinfo(refGenome))
}

getSiteSeqs <- function(gr, upstream_flank, downstream_flank, ref_genome){
  require(BSgenome)
  require(GenomicRanges)
  if(class(gr) != "GRanges"){
    if(class(gr) == "character" & all(stringr::str_detect(gr, ":"))){
      clusID <- stringr::str_split(gr, ":", simplify = TRUE)
      gr <- GRanges(
        seqnames = clusID[,1], 
        ranges = IRanges(start = as.numeric(clusID[,3]), width = 1), 
        strand = clusID[,2], seqinfo = seqinfo(ref_genome))
    }else{
      stop("Initial input must either be GRanges or posIDs delimited by ':'.")
  }}

  fl_starts <- GenomicRanges::flank(gr, width = -1, start = TRUE)
  seq_ranges <- GenomicRanges::trim(GenomicRanges::flank(
    GenomicRanges::shift(
      fl_starts, 
      shift = ifelse(
        strand(fl_starts) == "+", 
        downstream_flank, 
        -downstream_flank)),
    width = (upstream_flank + downstream_flank),
    start = TRUE))
  uniq_ranges <- unique(seq_ranges)
  
  seqs <- BSgenome::getSeq(ref_genome, uniq_ranges)
  seqs[match(seq_ranges, uniq_ranges)]
}

compareGuideRNAs <- function(gr_with_sequences, seq_col, 
                             guideRNASeqs, submat = NULL, tolerance = 6L,
                             upstream_flank, downstream_flank,
                             PAM = "NGG", offset_nt = 4L){ 
  #23 matching all, -6 mismatches in guide
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
  
  fwd_df <- alnGuideRNAs(seqs, guideRNASeqs, tolerance)
  
  rev_df <- alnGuideRNAs(reverseComplement(seqs), guideRNASeqs, tolerance)
  
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
  
  potential_sites$edit.site <- calcCutSite(
    potential_sites, matched_seqs, upstream_flank, 
    downstream_flank, PAM, offset_nt)
  
  if(length(good_alns) != length(sites)){
    all_sites <- c(potential_sites, non_probable_sites)
  }else{
    all_sites <- potential_sites
  }
  all_sites <- all_sites[order(all_sites$siteID)]
  if(any(duplicated(all_sites$siteID))){
    message("Some alignments with multiple guideRNA matches. 
            Total number of alignments expanded. See siteID.")
  }else{  
    all_sites$siteID <- NULL
  }
  all_sites
}

alnGuideRNAs <- function(seqs, guideRNASeqs, tolerance){
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

calcCutSite <- function(sites, matched_seqs, upstream_flank, 
                         downstream_flank, PAM, offset_nt){
  df <- data.frame(
    "chr" = seqnames(sites),
    "strand" = strand(sites),
    "pos" = start(flank(sites, -1, start = TRUE)),
    "guideRNA" = matched_seqs$guideRNA,
    "guide.aln" = ifelse(!grepl("(rev)", matched_seqs$guideRNA), "+", "-"),
    "guide.start" = matched_seqs$start,
    "guide.end" = matched_seqs$end)
  
  df$true.ort <- ifelse(
    as.character(df$strand) == as.character(df$guide.aln), "+", "-")
  
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

filter_inappropriate_comparisons <- function(guideRNA.match, specimen, 
                                             treatment){
  require(stringr)
  
  stopifnot(length(guideRNA.match) == length(specimen))
  
  df <- data.frame(
    order = 1:length(specimen),
    specimen = specimen,
    guideRNA = guideRNA.match)
  
  dfl <- split(df, df$specimen)

  edited_df <- do.call(rbind, mapply(function(d, treat){
        d$guideRNA <- ifelse(
          str_extract(d$guideRNA, "[\\w\\-\\.]+") %in% as.character(treat), 
          d$guideRNA, "No_valid_match")
        d
      }, d = dfl, treat = treatment[names(dfl)], SIMPLIFY = FALSE)) %>%
    arrange(order)
  
  edited_df$guideRNA
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
    gr$in_gene == "FALSE", 
    paste0(gr$gene_id_wo_annot, " "), 
    paste0(gr$gene_id_wo_annot, " *"))
  
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
    1,0,0,0,0,1,1,0,0,1,0,1,1,1,1,0,
    0,1,0,0,0,1,0,1,1,0,1,0,1,1,1,0,
    0,0,1,0,1,0,1,0,1,0,1,1,0,1,1,0,
    0,0,0,1,1,0,0,1,0,1,1,1,1,0,1,0,
    0,0,1,1,1,0,0,0,0,0,0,0,0,0,1,0,
    1,1,0,0,0,1,0,0,0,0,0,0,0,0,1,0,
    1,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,
    0,1,0,1,0,0,0,1,0,0,0,0,0,0,1,0,
    0,1,1,0,0,0,0,0,1,0,0,0,0,0,1,0,
    1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,
    0,1,1,1,0,0,0,0,0,0,1,0,0,0,1,0,
    1,0,1,1,0,0,0,0,0,0,0,1,0,0,1,0,
    1,1,0,1,0,0,0,0,0,0,0,0,1,0,1,0,
    1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
    ncol = 16,
    nrow = 16,
    byrow = TRUE,
    dimnames = list(
      c("A", "T", "G", "C", "S", "W", "R", "Y", 
        "K", "M", "B", "V", "H", "D", "N", "?"),
      c("A", "T", "G", "C", "S", "W", "R", "Y", 
        "K", "M", "B", "V", "H", "D", "N", "?")))
}
