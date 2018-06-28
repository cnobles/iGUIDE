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
  gr$order <- seq_along(gr)
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
  
  if(length(gr) == 0){ gr$grouping <- gr$order <- seq_along(gr) }
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
  clus.ori <- paste0(
    seqnames(gr), ":", strand(gr), ":", unlist(
      Rle(values = pile_starts, lengths = pile_ups$csize)))
  gr$clus.ori <- clus.ori[seq_along(gr)]
  gr
}  

identifyPairedAlgnmts <- function(gr, maxgap, grouping = NULL, maxovlp = 10L){
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
        qStart - sStart >= -maxovlp)
    
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

assignLociID <- function(gr, pilegap = 0L, pairgap = 200L, 
                         maxovlp = 10L, grouping = NULL){
  if(!is.null(grouping)){
    group_vec <- mcols(gr)[,match(grouping, names(mcols(gr)))]
  }else{
    group_vec <- rep(1, length(gr))
  }
  gr <- granges(gr)
  gr$grouping <- group_vec
  gr$ori.order <- seq_along(gr)
  gr$pile.id <- pileupCluster(
    gr, maxgap = pilegap, grouping = "grouping", return = "simple")
  gr$pair.id <- identifyPairedAlgnmts(
    gr, maxgap = pairgap, maxovlp = maxovlp, grouping = "grouping")
  gr$pair.id[is.na(gr$pair.id)] <- paste0("NA.", 1:sum(is.na(gr$pair.id)))
  
  grl <- split(gr, group_vec)
  gr_mod <- unlist(GenomicRanges::GRangesList(lapply(grl, function(gs){
    #red <- reduce(gs, min.gapwidth = pilegap, with.revmap = TRUE)
    pp <- as.data.frame(gs) %>% 
      distinct(pile.id, pair.id)
    pp_gr <- GenomicRanges::GRanges(
      seqnames = "mock", 
      ranges = IRanges::IRanges(
        start = as.integer(factor(pp$pair.id)), width = 1),
      strand = "+")
    pp_grl <- split(pp_gr, pp$pile.id)
    pp_el <- as.matrix(GenomicRanges::findOverlaps(pp_grl))
    pp_el <- matrix(names(pp_grl)[pp_el], ncol = 2)

    # Construct graph and identify clusters.
    g <- igraph::simplify(igraph::graph_from_edgelist(pp_el, directed = FALSE))
    mem <- igraph::membership(igraph::clusters(g))
    gs$mem <- mem[gs$pile.id]
    gs
  })))
  
  gr_mod$mem <- as.integer(factor(paste0(gr$grouping, ":", gr_mod$mem)))
  gr_mod <- gr_mod[order(gr_mod$ori.order)]
  return(gr_mod$mem)
}

# Format number in tables with big.marks conveinently
pNums <- function(x, ...){
  format(x, big.mark = ",", ...)
}

load_ref_files <- function(ref, type = "gene.list", freeze = NULL){
  stopifnot(type %in% c("gene.list", "GRanges", "data.frame"))
  suppressMessages(require("data.table"))
  suppressMessages(require("GenomicRanges"))
  
  if(grepl(".rds$", ref$file)){
    ref_set <- readRDS(ref$file)
  }else if(grepl(".RData$", ref$file)){
    ref_env <- new.env()
    load(config$refGenes, envir = refs)
    ref_set <- ref_env[[ls(ref_env)]]
  }else if(file.exists(ref$file) | grepl("^http", ref$file)){
    ref_set <- data.table::fread(ref$file, data.table = FALSE)
  }else{
    stopifnot(require("hiAnnotator"))
    stopifnot(grepl(":", ref$file))
    trackTable <- unlist(strsplit(ref$file, ":"))
    ucsc_session <- makeUCSCsession(freeze)
    
    stopifnot(
      trackTable[1] %in% rtracklayer::trackNames(ucsc_session) &
      trackTable[2] %in% rtracklayer::tableNames(
        rtracklayer::ucscTableQuery(ucsc_session, track = trackTable[1])))
    
    ref_tbl <- rtracklayer::getTable(rtracklayer::ucscTableQuery(
      ucsc_session, track = trackTable[1], table = trackTable[2]))
    
    ref_set <- rtracklayer::track(
      ucsc_session, name = trackTable[1], table = trackTable[2])
    
    stopifnot(all(ref_tbl$name == ref_set$name))
    ref_set <- granges(ref_set)
    mcols(ref_set) <- ref_tbl
  }
  
  if(!class(ref_set) %in% c("data.frame", "GRanges")){
    stop("Import of reference data failed. Check input parameters.")
  }
  
  if(type == "GRanges"){
    if(class(ref_set) == "GRanges"){
      ref_set$annot_sym <- mcols(ref_set)[,ref$symbol]
      return(ref_set)
    }else{
      ref_set$annot_sym <- ref_set[,ref$symbol]
      return(GenomicRanges::makeGRangesFromDataFrame(
        ref_set, keep.extra.columns = TRUE))
    }
  }else if(type == "gene.list"){
    ref_set <- try(as.data.frame(ref_set), silent = TRUE)
    if(class(ref_set) == "data.frame"){
      return(ref_set[,ref$symbol])
    }else{
      stop("Cannot coerce to geneList. Check input parameters.")
    }
  }else if(type == "data.frame"){
    ref_set <- try(as.data.frame(ref_set), silent = TRUE)
    if(class(ref_set) == "data.frame"){
      ref_set$annot_sym <- ref_set[,ref$symbol]
      return(ref_set)
    }else{
      stop("Cannot coerce to data.frame. Check input parameters.")
    }
  }
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
  rev_df$guideRNA <- paste0(rev_df$guideRNA, rep(" (rev)", nrow(rev_df)))
  
  matched_seqs <- rbind(fwd_df, rev_df)
  matched_seqs <- matched_seqs[grepl(pam, matched_seqs$aln.seq),]
  
  good_alns <- as.numeric(matched_seqs$names)
  
  if(length(good_alns) == 0){
    sites$guideRNA.match <- "No_valid_match"
    sites$guideRNA.mismatch <- NA
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
      d$guideRNA <- rep(names(alns[[i]][j]), nrow(d))
      d$mismatches <- rep(i-1, nrow(d))
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
                           onco_genes, special_genes, annotations = TRUE){
  require(GenomicRanges)
  require(hiAnnotator)
  stopifnot(length(seqnames) == length(positions))
  
  gr <- GRanges(
    seqnames = seqnames,
    ranges = IRanges(start = positions, width = rep(1, length(positions))),
    strand = rep("+", length(positions)),
    seqinfo = seqinfo(reference))
  
  if(length(gr) == 0){
    return(character())
  }else{
    # Annotate Sites with Gene names for within and nearest genes
    gr <- getSitesInFeature(
      gr, ref_genes, colnam = "in_gene", feature.colnam = "annot_sym")
    gr <- getNearestFeature(
      gr, ref_genes, colnam = "nearest_gene", feature.colnam = "annot_sym")
  
    ## Add gene marks ("*" for in_gene, "~" for onco gene, and "!" for special_genes)
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
      gr$gene_id_wo_annot %in% special_genes, paste0(gr$gene_id, "!"), gr$gene_id)
  
    if(annotations){
      return(gr$gene_id)
    }else{
      return(gr$gene_id_wo_annot)
    }
  }
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

make_square <- function(p, dims, fudge=1) {
  dims <- heatmap_dims(p)
  p + ggplot2::theme(aspect.ratio = (dims$nrows/dims$ncols)*fudge)
}

#' Combine a list of ShortRead objects
#' 
#' @param split.seqs list of ShortRead objects
#' @author Christopher Nobles, Ph.D.
serial_append_S4 <- function(split.seqs){
  require("ShortRead")
  stopifnot(class(split.seqs) == "list")
  
  app_env <- new.env()
  app_env$seqs <- split.seqs[[1]]
  
  null <- lapply(2:(length(split.seqs)), function(i){
    app_env$seqs <- append(app_env$seqs, split.seqs[[i]])
  })
  
  return(app_env$seqs)
}