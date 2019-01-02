## GUIDE-seq matching method
## Pile up alignments and find the beginning of the most frequent alignment
## (based on fragment lengths) for the group. Look +/- 25 bp around the site for
## identifying on/off-target sites with the guideRNA.
## This needs to be made into a function for standardizing locations, solve the 
## bias problem by doing it twice on the data, since it's so fast, it shouldn't
## matter the extra round. Further, this pile up method would be good for
## identifying genomic regions to scan and would reduce the query load, just 
## increase the range.
pileupCluster <- function(gr, grouping = NULL, maxgap = 0L, return = "full"){
  
  stopifnot(return %in% c("full", "simple", "ID"))
  
  if( !is.null(grouping) ){
    group_vec <- GenomicRanges::mcols(gr)[
        ,match(grouping, names(GenomicRanges::mcols(gr)))
      ]
  }else{
    group_vec <- rep(1, length(gr))
  }
  
  # Process each group individually to keep clustering separate
  gr$grouping <- group_vec
  gr$order <- seq_along(gr)
  grl <- split(gr, gr$grouping)
  gr <- unname(unlist(GenomicRanges::GRangesList(lapply(
    grl, 
    function(g){
      
      g_pos <- g[strand(g) == "+"]
      g_pos <- groupPileups(g_pos, "+", maxgap = maxgap)
      g_neg <- g[strand(g) == "-"]
      g_neg <- groupPileups(g_neg, "-", maxgap = maxgap)
      g_ret <- c(g_pos, g_neg)
      g_ret <- g_ret[order(g_ret$order)]
      g_ret$clusID <- as.integer(factor(g_ret$clus.ori))
      g_ret
      
    }
  ))))

  if( length(gr) == 0 ) gr$grouping <- gr$order <- seq_along(gr)
  
  gr$clus.ID <- as.integer(factor(gr$clus.ori))
  gr <- gr[order(gr$order)]
  gr$order <- gr$grouping <- NULL
  
  if( return == "full" ){
    return(list("gr" = gr, "clusID" = gr$clus.ID, "clusOri" = gr$clus.ori))
  }else if( return == "simple" ){
    return(gr$clus.ori)
  }else if( return == "ID" ){
    return(gr$clus.ID)
  }
  
}

groupPileups <- function(gr, strand, maxgap = maxgap){
  
  # Implement axial cluster structure rather than all vs all.
  red.gr <- GenomicRanges::reduce(gr, min.gapwidth = maxgap, with.revmap = TRUE)
  
  axil_nodes <- unlist(S4Vectors::Rle(
    values = unlist(red.gr$revmap)[
      S4Vectors::start(red.gr$revmap@partitioning)
    ], 
    lengths = S4Vectors::width(red.gr$revmap@partitioning)
  ))
  
  nodes <- unlist(red.gr$revmap)
  edgelist <- unique(matrix( c(axil_nodes, nodes), ncol = 2 ))
  pile_ups <- igraph::clusters(igraph::graph.edgelist(
    el = edgelist, directed = FALSE
  ))
  gr$clusID <- igraph::membership(pile_ups)
  
  if( strand == "+" ){
    
    pile_starts <- sort(split(
      x = S4Vectors::Rle(GenomicRanges::start(gr)), 
      f = gr$clusID
    ))
    
    pile_starts <- S4Vectors::runValue(pile_starts)[
      S4Vectors::runLength(pile_starts) == 
        max(S4Vectors::runLength(pile_starts))
    ]
    
    pile_starts <- unlist(pile_starts)[
      S4Vectors::start(pile_starts@partitioning)
    ]
    
  }else{
    
    pile_starts <- sort(split(
      x = S4Vectors::Rle(GenomicRanges::end(gr)), 
      f = gr$clusID
    ))
    
    pile_starts <- S4Vectors::runValue(pile_starts)[
      S4Vectors::runLength(pile_starts) == 
        max(S4Vectors::runLength(pile_starts))
    ]
    
    pile_starts <- unlist(pile_starts)[
      S4Vectors::end(pile_starts@partitioning)
    ]
    
  }
  
  gr <- gr[order(gr$clusID)]
  
  clus.ori <- paste0(
    GenomicRanges::seqnames(gr), ":", GenomicRanges::strand(gr), ":", 
    unlist(S4Vectors::Rle(values = pile_starts, lengths = pile_ups$csize))
  )
  
  gr$clus.ori <- clus.ori[seq_along(gr)]
  
  gr
  
}  

identifyPairedAlgnmts <- function(gr, grouping = NULL, maxgap, maxovlp = 10L){
  
  if( !is.null(grouping) ){
    group_vec <- GenomicRanges::mcols(gr)[
      ,match(grouping, names(GenomicRanges::mcols(gr)))
    ]
  }else{
    group_vec <- rep(1, length(gr))
  }
  
  gr$ori.order <- seq_along(gr)
  grl <- split(x = gr, f = group_vec)
  
  pairs <- unname(unlist(lapply(
    seq_along(grl), 
    function(i, grl, maxgap){
      
      gs <- grl[[i]]
      
      red <- GenomicRanges::reduce(
        x = GenomicRanges::flank(gs, width = -1, start = TRUE), 
        min.gapwidth = 0L, 
        with.revmap = TRUE
      )
    
      hits <- GenomicRanges::findOverlaps(
        red, maxgap = maxgap, ignore.strand = TRUE
      )
      
      hits <- as.data.frame(hits, row.names = NULL) %>%
        dplyr::mutate(
          qStart = GenomicRanges::start(red[queryHits]),
          sStart = GenomicRanges::start(red[subjectHits]),
          qStrand = as.character(GenomicRanges::strand(red[queryHits])),
          sStrand = as.character(GenomicRanges::strand(red[subjectHits]))
        ) %>%
        dplyr::filter(
          queryHits != subjectHits, 
          qStrand != sStrand,
          qStrand == "+",
          qStart - sStart >= -maxovlp
        )
      
      if( nrow(hits) == 0 ){
        
        return(rep(NA, length(gs)))
        
      }else{
        
        el <- hits[,c("queryHits", "subjectHits")]
        
        g <- igraph::simplify(igraph::graph_from_edgelist(
          el = as.matrix(el), 
          directed = FALSE
        ))
        
        mem <- data.frame(
            idx = seq_len(igraph::vcount(g)), 
            grp = as.numeric(igraph::membership(igraph::clusters(g)))) %>%
          dplyr::group_by(grp) %>%
          dplyr::filter(n() > 1) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(grp = as.integer(as.factor(grp))) %>%
          as.data.frame()
        
        mem$revmap <- red$revmap[mem$idx]
        
        mem_exp <- data.frame(
          gs_idx = unlist(mem$revmap),
          grp = paste0(
            names(grl)[i], ".", 
            S4Vectors::Rle(values = mem$grp, lengths = lengths(mem$revmap))
          )
        )
        
        ret <- rep(NA, length(gs))
        ret[mem_exp$gs_idx] <- mem_exp$grp
        return(ret)
        
      }
    }, 
    grl = grl, 
    maxgap = maxgap
  )))
  
  if( is.null(grouping) ){
    pairs <- as.numeric(stringr::str_extract(pairs, "[0-9]+$"))
  }
  
  un_grl <- unlist(grl)
  pairs[order(un_grl$ori.order)]
  
}

assignLociID <- function(gr, pilegap = 0L, pairgap = 200L, maxovlp = 10L, 
                         grouping = NULL){
  
  if( !is.null(grouping) ){
    group_vec <- GenomicRanges::mcols(gr)[
      ,match(grouping, names(GenomicRanges::mcols(gr)))
    ]
  }else{
    group_vec <- rep(1, length(gr))
  }
  
  gr <- GenomicRanges::granges(gr)
  gr$grouping <- group_vec
  gr$ori.order <- seq_along(gr)
  
  gr$pile.id <- pileupCluster(
    gr, maxgap = pilegap, grouping = "grouping", return = "simple"
  )
  
  gr$pair.id <- identifyPairedAlgnmts(
    gr, maxgap = pairgap, maxovlp = maxovlp, grouping = "grouping"
  )
  
  gr$pair.id[is.na(gr$pair.id)] <- paste0(
    "NA.", seq_len(sum(is.na(gr$pair.id)))
  )
  
  grl <- split(x = gr, f = group_vec)
  gr_mod <- unlist(GenomicRanges::GRangesList(lapply(
    grl, 
    function(gs){
    
      pp <- as.data.frame(gs) %>% 
        dplyr::distinct(pile.id, pair.id)
      
      pp_gr <- GenomicRanges::GRanges(
        seqnames = "mock", 
        ranges = IRanges::IRanges(
          start = as.integer(factor(pp$pair.id)), width = 1
        ),
        strand = "+"
      )
      
      pp_grl <- split(x = pp_gr, f = pp$pile.id)
      pp_el <- as.matrix(GenomicRanges::findOverlaps(pp_grl))
      pp_el <- matrix(names(pp_grl)[pp_el], ncol = 2)
  
      # Construct graph and identify clusters.
      g <- igraph::simplify(igraph::graph_from_edgelist(
        el = pp_el, directed = FALSE
      ))
      
      mem <- igraph::membership(igraph::clusters(g))
      gs$mem <- mem[gs$pile.id]
      gs
      
    }
  )))

  gr_mod$mem <- as.integer(factor(paste0(gr$grouping, ":", gr_mod$mem)))
  gr_mod <- gr_mod[order(gr_mod$ori.order)]
  return(gr_mod$mem)
  
}

# Format number in tables with big.marks conveinently
pNums <- function(x, ...){
  format(x, big.mark = ",", ...)
}

loadRefFiles <- function(ref, type = "gene.list", freeze = NULL){
  
  stopifnot(type %in% c("gene.list", "GRanges", "data.frame"))

  if( grepl(".rds$", ref$file) ){
    
    ref_set <- readRDS(ref$file)
    
  }else if( grepl(".RData$", ref$file) ){
    
    ref_env <- new.env()
    load(config$refGenes, envir = refs)
    ref_set <- ref_env[[ls(ref_env)]]
    
  }else if( file.exists(ref$file) | grepl("^http", ref$file) ){
    
    ref_set <- data.table::fread(ref$file, data.table = FALSE)
    
  }else{
    
    stopifnot( "hiAnnotator" %in% row.names(installed.packages()) )
    stopifnot( grepl(":", ref$file) )
    trackTable <- unlist(strsplit(ref$file, ":"))
    ucsc_session <- hiAnnotator::makeUCSCsession(freeze)
    
    stopifnot(
      trackTable[1] %in% rtracklayer::trackNames(ucsc_session) &
      trackTable[2] %in% rtracklayer::tableNames(
        rtracklayer::ucscTableQuery(ucsc_session, track = trackTable[1])
      )
    )
    
    ref_tbl <- rtracklayer::getTable(rtracklayer::ucscTableQuery(
      ucsc_session, track = trackTable[1], table = trackTable[2]
    ))
    
    ref_set <- rtracklayer::track(
      ucsc_session, name = trackTable[1], table = trackTable[2]
    )
    
    stopifnot( all(ref_tbl$name == ref_set$name) )
    ref_set <- GenomicRanges::granges(ref_set)
    GenomicRanges::mcols(ref_set) <- ref_tbl
    
  }
  
  if( !class(ref_set) %in% c("data.frame", "GRanges") ){
    stop("Import of reference data failed. Check input parameters.")
  }
  
  if( type == "GRanges" ){
    
    if( class(ref_set) == "GRanges" ){
      ref_set$annot_sym <- GenomicRanges::mcols(ref_set)[,ref$symbol]
      return(ref_set)
    }else{
      ref_set$annot_sym <- ref_set[,ref$symbol]
      return(GenomicRanges::makeGRangesFromDataFrame(
        ref_set, keep.extra.columns = TRUE
      ))
    }
    
  }else if( type == "gene.list" ){
    
    ref_set <- try(as.data.frame(ref_set), silent = TRUE)
    
    if( class(ref_set) == "data.frame" ){
      return(ref_set[,ref$symbol])
    }else{
      stop("Cannot coerce to geneList. Check input parameters.")
    }
    
  }else if( type == "data.frame" ){
    
    ref_set <- try(as.data.frame(ref_set), silent = TRUE)
    
    if( class(ref_set) == "data.frame" ){
      ref_set$annot_sym <- ref_set[,ref$symbol]
      return(ref_set)
      
    }else{
      
      stop("Cannot coerce to data.frame. Check input parameters.")
      
    }
    
  }
  
}

# Convert matches to GRanges objects
mindexToGranges <- function(mindex, strand, ref = NULL){
  
  if( is.null(mindex@NAMES) ) stop("NAMES column not found for seqnames.")
  ir <- unlist(mindex)
  
  if( is.null(ref) ){
    
    return(unname(GenomicRanges::GRanges(
      seqnames = names(ir),
      ranges = ir,
      strand = rep(strand, length(ir))
    )))
    
  }else{
    
    return(unname(GenomicRanges::GRanges(
      seqnames = names(ir),
      ranges = ir,
      strand = rep(strand, length(ir)),
      seqinfo = GenomeInfoDb::seqinfo(ref)
    )))
    
  }
  
}

# Generate random sites across the reference genome.
selectRandomSites <- function(num, ref.genome, drop.extra.seqs = TRUE, 
                              seq.names = NULL, rnd.seed = NULL){

  if( !is.null(rnd.seed) ) set.seed(rnd.seed)
  
  if( is.null(GenomeInfoDb::seqinfo(ref.genome)) ){
    stop("Ref genome does not have seqinfo.")
  }
  
  if( !is.null(seq.names) ){
    
    seq_names <- seq.names
    
    seq_lengths <- ref.genome@seqinfo@seqlengths[
      match(seq.names, GenomeInfoDb::seqnames(ref.genome))
    ]
    
    names(seq_lengths) <- seq.names
    
  }else{
    
    if( drop.extra.seqs ){
      
      seq_names <- grep(
        pattern = "chr[0-9XY]+$", 
        x = GenomeInfoDb::seqnames(ref.genome), 
        value = TRUE
      )
      
      seq_lengths <- ref.genome@seqinfo@seqlengths[
        match(seq_names, seqnames(ref.genome))
      ]
      
      names(seq_lengths) <- seq_names
      
    }else{
      
      seq_names <- GenomeInfoDb::seqnames(ref.genome)
      seq_lengths <- ref.genome@seqinfo@seqlengths
      names(seq_lengths) <- seq_names
      
    }
    
  }
  
  chrs <- sort(factor(
    sample(seq_names, num, replace = TRUE, prob = seq_lengths),
    levels = seq_names
  ))
  
  strands <- sample(c("+", "-"), num, replace = TRUE)
  splt_chrs <- split(chrs, chrs)
  splt_chrs <- splt_chrs[ sapply(splt_chrs, length) > 0 ]
  
  positions <- unname(unlist(sapply(
    splt_chrs, 
    function(seq, seq_leng){
      
      sample.int(
        n = seq_leng[unique(as.character(seq))], 
        size = length(seq), 
        replace = FALSE, 
        useHash = FALSE
      )
      
    }, 
    seq_leng = seq_lengths
  )))
  
  GenomicRanges::GRanges(
    seqnames = chrs,
    ranges = IRanges::IRanges(start = positions, width = 1L),
    strand = strands,
    seqinfo = GenomeInfoDb::seqinfo(ref.genome)
  )
  
}

getSiteSeqs <- function(gr, upstream.flank, downstream.flank, ref.genome){

  if( class(gr) != "GRanges" ){
    
    if( class(gr) == "character" & all(stringr::str_detect(gr, ":")) ){
      
      clusID <- stringr::str_split(gr, ":", simplify = TRUE)
      
      gr <- GenomicRanges::GRanges(
        seqnames = clusID[,1], 
        ranges = IRanges::IRanges(start = as.numeric(clusID[,3]), width = 1), 
        strand = clusID[,2], 
        seqinfo = GenomicRanges::seqinfo(ref.genome)
      )
      
    }else{
      
      stop("Initial input must either be GRanges or posIDs delimited by ':'.")
      
    }
    
  }

  fl_starts <- GenomicRanges::flank(gr, width = -1, start = TRUE)
  
  seq_ranges <- GenomicRanges::trim(GenomicRanges::flank(
    GenomicRanges::shift(
      fl_starts, 
      shift = ifelse(
        GenomicRanges::strand(fl_starts) == "+", 
        downstream.flank, 
        -downstream.flank
      )
    ),
    width = upstream.flank + downstream.flank,
    start = TRUE
  ))
  
  uniq_ranges <- unique(seq_ranges)
  
  seqs <- BSgenome::getSeq(ref.genome, uniq_ranges)
  seqs[match(seq_ranges, uniq_ranges)]
  
}

compareGuideRNAs <- function(gr.with.sequences, seq.col, 
                             guide.rna.seqs, submat = NULL, tolerance = 6L,
                             upstream.flank, downstream.flank,
                             PAM = "NGG", offset.nt = 4L){ 
  
  #23 matching all, -6 mismatches in guide

  if( is.null(submat) ) submat <- banmat()
  
  pam <- paste0(gsub("N", ".", PAM), "$")
  
  sites <- gr.with.sequences
  sites$siteID <- seq_len(length(sites))
  seq_col_match <- match(seq.col, names(GenomicRanges::mcols(sites)))
  
  if( length(seq_col_match) == 0 ){
    
    stop("Cannot find sequences in column.")
    
  }else{
    
    seqs <- GenomicRanges::mcols(sites)[,seq_col_match]
    names(seqs) <- sites$siteID
    
  }
  
  fwd_df <- alnGuideRNAs(seqs, guide.rna.seqs, tolerance)
  
  rev_df <- alnGuideRNAs(reverseComplement(seqs), guide.rna.seqs, tolerance)
  
  fwd_df$aln.seq <- as.character(Biostrings::DNAStringSet(
    x = seqs[fwd_df$names], 
    start = fwd_df$start, 
    end = fwd_df$end
  ))
  
  rev_df$aln.seq <- as.character(Biostrings::DNAStringSet(
    x = Biostrings::reverseComplement(seqs[rev_df$names]), 
    start = rev_df$start, 
    end = rev_df$end
  ))
  
  rev_df$guideRNA <- paste0(rev_df$guideRNA, rep(" (rev)", nrow(rev_df)))
  
  matched_seqs <- rbind(fwd_df, rev_df)
  matched_seqs <- matched_seqs[ grepl(pam, matched_seqs$aln.seq), ]
  
  good_alns <- as.numeric(matched_seqs$names)
  
  if( length(good_alns) == 0 ){
    
    sites$guideRNA.match <- "No_valid_match"
    sites$guideRNA.mismatch <- NA
    sites$guideRNA.score <- NA
    sites$aligned.sequence <- NA
    sites$edit.site <- NA
    return(sites)
    
  }
  
  if( length(good_alns) != length(sites) ){
    
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
    potential_sites, matched_seqs, upstream.flank, 
    downstream.flank, PAM, offset.nt)
  
  if( length(good_alns) != length(sites) ){
    all_sites <- c(potential_sites, non_probable_sites)
  }else{
    all_sites <- potential_sites
  }
  
  all_sites <- all_sites[order(all_sites$siteID)]
  
  if( any(duplicated(all_sites$siteID)) ){
    message("Some alignments with multiple guideRNA matches. 
            Total number of alignments expanded. See siteID.")
  }else{  
    all_sites$siteID <- NULL
  }
  
  all_sites
  
}

alnGuideRNAs <- function(seqs, guide.rna.seqs, tolerance){

  if( is.null(names(seqs)) ) names(seqs) <- seq_along(seqs)
  
  nt_widths <- data.frame(
    names = names(seqs),
    nt_width = Biostrings::width(seqs)
  )
  
  alns <- lapply(
    0:tolerance, 
    function(tol, guides, seqs){
      
      lapply(
        guides, 
        function(guide){
          
          unlist(Biostrings::vmatchPattern(
            pattern = guide, 
            subject = seqs, 
            max.mismatch = tol, 
            fixed = 'subject'
          ))
          
        }
      )
      
    }, 
    guides = guide.rna.seqs, 
    seqs = seqs
  )
  
  df <- do.call(
    rbind, 
    lapply(
      seq_along(alns), 
      function(i){
        
        do.call(
          rbind, 
          lapply(
            seq_along(alns[[i]]), 
            function(j){
              
              d <- as.data.frame(alns[[i]][[j]])
              d$guideRNA <- rep(names(alns[[i]][j]), nrow(d))
              d$mismatches <- rep(i-1, nrow(d))
              d
              
            }
          )
        )
        
      }
    )
  )
        

  dplyr::group_by(df, start, end, width, names, guideRNA) %>%
    dplyr::mutate(guideRNA.mismatch = min(mismatches)) %>%
    dplyr::ungroup() %>%
    dplyr::select(names, guideRNA, guideRNA.mismatch, start, end, width) %>%
    dplyr::distinct() %>%
    dplyr::left_join(., nt_widths, by = "names") %>%
    dplyr::mutate(
      start = ifelse(start <= 0, 1, start),
      end = ifelse(end > nt_width, nt_width, end),
      width = end - start + 1
    )
  
}

calcCutSite <- function(sites, matched.seqs, upstream.flank, 
                         downstream.flank, PAM, offset.nt){
  df <- data.frame(
    "chr" = GenomicRanges::seqnames(sites),
    "strand" = GenomicRanges::strand(sites),
    "pos" = GenomicRanges::start(GenomicRanges::flank(sites, -1, start = TRUE)),
    "guideRNA" = matched.seqs$guideRNA,
    "guide.aln" = ifelse(!grepl("(rev)", matched.seqs$guideRNA), "+", "-"),
    "guide.start" = matched.seqs$start,
    "guide.end" = matched.seqs$end
  )
  
  df$true.ort <- ifelse(
    as.character(df$strand) == as.character(df$guide.aln), "+", "-"
  )
  
  df$edit.pos <- ifelse(
    df$strand == "+",
    ifelse(
      df$guide.aln == "+",
      df$pos - upstream.flank + df$guide.end - nchar(PAM) - offset.nt,
      df$pos - df$guide.end + downstream.flank + nchar(PAM) + offset.nt - 1),
    ifelse(
      df$guide.aln == "+",
      df$pos + upstream.flank - df$guide.end + nchar(PAM) + offset.nt,
      df$pos + df$guide.end - downstream.flank - nchar(PAM) - offset.nt + 1
    )
  )
  
  return(paste0(df$chr, ":", df$true.ort, ":", df$edit.pos))
  
}

filterInappropriateComparisons <- function(guideRNA.match, specimen, treatment){

  stopifnot( length(guideRNA.match) == length(specimen) )
  
  df <- data.frame(
    order = seq_along(specimen),
    specimen = specimen,
    guideRNA = guideRNA.match
  )
  
  dfl <- split(x = df, f = df$specimen)

  edited_df <- dplyr::bind_rows(
      mapply(
        function(d, treat){
          
          d$guideRNA <- ifelse(
            stringr::str_extract(d$guideRNA, "[\\w\\-\\.]+") %in% 
              as.character(treat), 
            d$guideRNA, 
            "No_valid_match"
          )
          d
          
        }, 
        d = dfl, 
        treat = treatment[names(dfl)], 
        SIMPLIFY = FALSE
      )
    ) %>%
    dplyr::arrange(order)
  
  edited_df$guideRNA
  
}

assignGeneID <- function(seqnames, positions, reference, ref.genes, onco.genes,
                         special.genes, annotations = TRUE){

  stopifnot( length(seqnames) == length(positions) )
  
  gr <- GenomicRanges::GRanges(
    seqnames = seqnames,
    ranges = IRanges::IRanges(
      start = positions, width = rep(1, length(positions))
    ),
    strand = rep("+", length(positions)),
    seqinfo = GenomicRanges::seqinfo(reference)
  )
  
  if( length(gr) == 0 ){
    
    return(character())
    
  }else{
    
    # Annotate Sites with Gene names for within and nearest genes
    gr <- hiAnnotator::getSitesInFeature(
      gr, ref.genes, colnam = "in_gene", feature.colnam = "annot_sym"
    )
    
    gr <- hiAnnotator::getNearestFeature(
      gr, ref.genes, colnam = "nearest_gene", feature.colnam = "annot_sym"
    )
  
    ## Add gene marks 
    ## ("*" for in_gene, "~" for onco gene, and "!" for special.genes)
    
    gr$gene_id_wo_annot <- ifelse(
      gr$in_gene == "FALSE", gr$nearest_gene, gr$in_gene
    )
    
    gr$gene_id_wo_annot <- sapply(strsplit(gr$gene_id_wo_annot, ","), "[[", 1)
  
    gr$gene_id <- ifelse(
      gr$in_gene == "FALSE", 
      gr$gene_id_wo_annot, 
      paste0(gr$gene_id_wo_annot, "*")
    )
  
    gr$gene_id <- ifelse(
      gr$gene_id_wo_annot %in% onco.genes, 
      paste0(gr$gene_id, "~"), 
      gr$gene_id
    )
  
    gr$gene_id <- ifelse(
      gr$gene_id_wo_annot %in% special.genes, 
      paste0(gr$gene_id, "!"), 
      gr$gene_id
    )
  
    if( annotations ){
      return(gr$gene_id)
    }else{
      return(gr$gene_id_wo_annot)
    }
    
  }
  
}

calcCoverage <- function(gr, resolution){ 
  ###!!!!!!!!Add option for counting coverage by reads or uniq frags
  #Set up coverage gr
  strandless <- gr
  GenomicRanges::strand(strandless) <- "*"
  gr_ranges <- range(strandless)
  
  window_seqs <- lapply(
    gr_ranges, 
    function(chr, res){
      seq(GenomicRanges::start(chr), GenomicRanges::end(chr), res)
    }, 
    res = resolution
  )
  
  coverage_grl <- GenomicRanges::GRangesList(lapply(
    seq_along(gr_ranges), 
    function(i, gr_ranges, window_seqs){
      
      seqname <- GenomicRanges::seqnames(gr_ranges[i])
      window <- window_seqs[[i]]
      
      GenomicRanges::GRanges(
        seqnames = rep(seqname, length(window)),
        ranges = IRanges::IRanges(
          start = window, width = rep(resolution, length(window))),
        strand = rep("*", length(window)))
      
    }, 
    gr_ranges = gr_ranges, 
    window_seqs = window_seqs
  ))
  
  coverage_pos <- coverage_grl
  
  coverage_pos <- GenomicRanges::GRangesList(lapply(
    coverage_pos, 
    function(x){
      strand(x) <- rep("+", length(x))
      x
    }
  ))
  
  coverage_neg <- coverage_grl
  
  coverage_neg <- GenomicRanges::GRangesList(lapply(
    coverage_pos, 
    function(x){
      strand(x) <- rep("-", length(x))
      x
    }
  ))
  
  dplyr::bind_rows(lapply(
    seq_along(coverage_grl), 
    function(i, gr){
      
      as.data.frame(coverage_grl[[i]], row.names = NULL) %>%
        select(seqnames, start, end, width) %>%
        mutate(
          readCountsPos = countOverlaps(coverage_pos[[i]], gr),
          readCountsNeg = countOverlaps(coverage_neg[[i]], gr)
        ) %>%
        arrange(seqnames)
      
    }, 
    gr = gr
  ))
  
}

plotCoverage <- function(gr, resolution = 10L){
  
  df <- calcCoverage(gr, resolution)
  
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

plotEditSites <- function(gr, sampleName = NULL, resolution = 10L){
  
  stopifnot( length(unique(gr$edit.site)) == 1 )
  
  if( !is.null(sampleName) ){
    
    isThere <- match(sampleName, names(GenomicRanges::mcols(gr)))
    if( length(isThere) == 0 ) stop("SampleName column not found.")
    sample <- unique(as.character(GenomicRanges::mcols(gr)[,isThere]))
    
  }else{
    
    sample <- "Unspecified"
    
  }
  
  guide <- unique(stringr::str_extract(gr$guideRNA.match, "[\\w\\-\\.]+"))
  edit_site <- as.character(unique(gr$edit.site))
  edit_site <- unlist(strsplit(edit_site, ":"))
  edit_pos <- as.numeric(edit_site[3])
  edit_site[3] <- format(edit_pos, big.mark = ",", scientific = FALSE)
  
  df <- calcCoverage(gr, resolution)
  
  df$edit.site <- paste0(
    "Sample: ", sample, 
    "\nGuide: ", guide,
    "\nEdit Site: ", paste(edit_site, collapse = ":")
  )
  
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

make_square <- function(p, dims, fudge=1){
  
  dims <- heatmap_dims(p)
  p + ggplot2::theme(aspect.ratio = (dims$nrows / dims$ncols) * fudge)
  
}

#' Combine a list of ShortRead objects
#' 
#' @param split.seqs list of ShortRead objects
#' @author Christopher Nobles, Ph.D.
serialAppendS4 <- function(split.seqs){

  stopifnot(class(split.seqs) == "list")
  
  app_env <- new.env()
  app_env$seqs <- split.seqs[[1]]
  
  null <- lapply(
    2:(length(split.seqs)), 
    function(i) app_env$seqs <- ShortRead::append(app_env$seqs, split.seqs[[i]])
  )
  
  return(app_env$seqs)
}

#' Cluster or group key - value(s) pairs
#' 
#' @usage clusterKV(key, val)
#' @usage clusterKV(key, val, return = "data.frame")
#' 
#' @description Function for clustering keys based on value content. Clustering
#' or grouping is solely based on presence / absence of values across keys. For
#' example, if we have the key - values of A = c(1, 2, 3), B = c(3, 4, 5), and
#' C = c(6, 7, 8), then keys A and B will be clustered together because they 
#' share the value 3, but C will not be clustered with either A or B because it
#' does not share any values with the respective keys. This type of clustering
#' can be helpful for grouping keys based on unique IDs, such as readnames or
#' character strings representing unique alignment locations.
#' 
#' @param key,val vector coercible into a factor vector. Both key and val 
#' vectors need to be equal length. Output will be of equal length and in same 
#' input order.
#' @param return options for output returned. "standard" will return a numeric 
#' vector of grouping IDs and is the default. "data.frame" will return a 
#' data.frame with key, val, and clus columns. "simple" will return a numeric
#' vector of grouping IDs with the names associated with unique keys. Lastly,
#' "graph" will return a simplifed graph with the keys as nodes and edges 
#' indicating which keys share values.
#' 
#' @author Christopher Nobles, Ph.D.
#' 
#' @importFrom magrittr %>%
#' 
clusterKV <- function(key, val, return = "standard"){
  
  # Check inputs
  stopifnot( "package:magrittr" %in% search())
  stopifnot( return %in% c("standard", "data.frame", "simple", "graph") )
  stopifnot( length(key) == length(val) )
  
  # Factorize keys and values for consistancy of indexing
  key_fac <- factor(key)
  val_fac <- factor(val)
  
  # Construct mock GRangesList where positions represent indices
  grl <- GenomicRanges::GRanges(
      seqnames = "mock", 
      ranges = IRanges::IRanges(start = as.integer(val_fac), width = 1), 
      strand = "*"
    ) %>%
    GenomicRanges::split(key) %>%
    GenomicRanges::reduce()
  
  # Determine which keys overlap based on associated values
  g <- GenomicRanges::findOverlaps(grl) %>%
    as.matrix() %>%
    igraph::graph.edgelist(directed = FALSE) %>%
    igraph::simplify()
  
  clus <- igraph::clusters(g)
  
  # Return data in different formats
  if( return == "standard" ){
    
    return(igraph::membership(clus)[as.integer(key_fac)])
    
  }else if( return == "data.frame" ){
    
    return(data.frame(
      key = key, 
      val = val, 
      clus = igraph::membership(clus)[as.integer(key_fac)]
    ))
    
  }else if( return == "simple" ){
    
    return(structure(igraph::membership(clus), names = levels(key_fac)))
    
  }else if( return == "graph" ){
    
    return(g)
    
  }
  
}

#' Predict Edit Site Probablility based on incorporation distance given an
#' On-target incorporation density
#' 
#' @usage predictESProb(x, density)
#' 
#' @param x integer position within range of density object indicating the 
#' distance from the predicted edit site.
#' @param density a density object constructed from the incorporation site 
#' distribution around associated On-target editing site(s).
#' @param range a numeric vector indicating the range of data to consider.
#' 
#' @author Christopher Nobles, Ph.D.
#' 
predictESProb <- function(x, density, range = NULL){
  
  if( is.null(range) & class(density) == "density" ) range <- range(density$x)
  
  if(suppressWarnings(
    x >= min(range) & x <= max(range) & class(density) == "density"
  )){
    
      idx <- match(x, round(density$x))
      return(
        with(density, 1 - sum(y[seq_len(idx)] * rep(diff(x)[1], idx)))
      )
      
  }else{
    
    return(NA)
    
  }
  
}
