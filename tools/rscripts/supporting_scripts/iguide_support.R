#' Cluster ranges into pile up clusters
#' 
#' @usage pileupCluster(gr)
#' @usage pileupCluster(gr, grouping = NULL, maxgap = 0L, return = "full")
#' 
#' @param gr a GRanges object to identify clusters with. If the grouping param
#' is included, the GRanges object should contain an column in the metadata 
#' columns with a name that matches the grouping param.
#' 
#' @param grouping a character vector of column names used for grouping ranges
#' together, i.e. clustering within specimens by use `grouping = "specimen"` 
#' where `specimen` is a metadata column in the input `gr`.
#' 
#' @param maxgap an integer indicating the distance allowed between ranges to 
#' consider them as overlapping. Follows the GenomicRanges::findOverlaps
#' notation, where 0L means adjacent ranges (not overlapping) will be identified
#' as overlapping. To only select ranges which overlap, use -1L.
#' 
#' @param return a character string of either "full" (default), "simple", or 
#' "ID". Full will return a listed object containing the input gr (`$gr`) with
#' appended data columns with information about the clusters. A "simple" return
#' returns a character vector of position identifiers for the considered origin 
#' of the cluster, while "ID" returns an integer vector of the cluster 
#' membership.
#' 
#' @description Pile up ranges and find the beginning of the most frequent 
#' range (based on fragment lengths) for the group. 
#' 
pileupCluster <- function(gr, grouping = NULL, maxgap = -1L, return = "full"){
  
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

#' groupPileup ranges, companion function to pileupCluster
#' 
#' @usage groupPileups(gr, strand, maxgap)
#' 
#' @param gr a GRanges object to identify clusters with. If the grouping param
#' is included, the GRanges object should contain an column in the metadata 
#' columns with a name that matches the grouping param.
#' 
#' @param strand character vector of singular length specifying either "+" or 
#' "-" strand to identify the grouping on. 
#' 
#' @param maxgap integer indicating the distance allowed between ranges to 
#' consider them as overlapping. Follows the GenomicRanges::findOverlaps
#' notation, where 0L means adjacent ranges (not overlapping) will be identified
#' as overlapping. To only select ranges which overlap, use -1L.
#' 
#' @description Identify groups of ranges that 'pileup' or overlap together.
#' 
groupPileups <- function(gr, strand, maxgap){
  
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

#' Identify paired alignments within a set of ranges
#' 
#' @usage identifyPairedAlgnmts(gr, grouping = NULL, maxgap, maxovlp = 10L)
#' 
#' @param gr a GRanges object to identify paired alignments. If the grouping 
#' param is included, the GRanges object should contain an column in the 
#' metadata columns with a name that matches the grouping param.
#' 
#' @param grouping a character vector of column names used for grouping ranges
#' together, i.e. clustering within specimens by use `grouping = "specimen"` 
#' where `specimen` is a metadata column in the input `gr`.
#' 
#' @param maxgap integer indicating the distance allowed between ranges to 
#' consider them as overlapping. Follows the GenomicRanges::findOverlaps
#' notation, where 0L means adjacent ranges (not overlapping) will be identified
#' as overlapping. To only select ranges which overlap, use -1L.
#'
#' @param maxovlp integer indicating the maximum allowed overlap before breaking
#' the paired criteria, this is implemented to give some control over alignments
#' that overlap more then would be expected by biology.
#' 
#' @description Given a set of ranges, identify which ranges match with the 
#' paired criteria. This means that the ranges are outward facing from eachother
#' and are found on opposite strands.
#' 
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

#' Assign locus IDs to alignments
#' 
#' @usage assignLociID(gr, pilegap = 0L, pairgap = 200L, maxovlp = 10L, grouping = NULL)
#' 
#' @param gr a GRanges object to identify pileup and paired alignments. If the 
#' grouping param is included, the GRanges object should contain an column in 
#' the metadata columns with a name that matches the grouping param.
#' 
#' @param pilegap,pairgap an integer indicating the distance allowed between 
#' ranges to consider them as overlapping. `pilegap` is passed to 
#' `pileupCluster` while `pairgap` is passed to `identifyPairedAlgnmts`.Follows 
#' the GenomicRanges::findOverlaps notation, where 0L means adjacent ranges 
#' (not overlapping) will be identified as overlapping. To only select ranges 
#' which overlap, use -1L.
#' 
#' @param pairgap integer indicating the distance allowed between ranges to 
#' consider them as overlapping. Follows the GenomicRanges::findOverlaps
#' notation, where 0L means adjacent ranges (not overlapping) will be identified
#' as overlapping. To only select ranges which overlap, use -1L.
#' 
#' @param maxovlp integer indicating the maximum allowed overlap before breaking
#' the paired criteria, this is implemented to give some control over alignments
#' that overlap more then would be expected by biology.
#' 
#' @param grouping a character vector of column names used for grouping ranges
#' together, i.e. clustering within specimens by use `grouping = "specimen"` 
#' where `specimen` is a metadata column in the input `gr`.
#' 
#' @description Using the pileup and paired identification, give each set of 
#' alignments a unique identifier that groups them together for analysis.
#' 
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

#' Format number in tables with big.marks conveinently
#' 
#' @usage pNums(x, ...)
#' 
#' @param x vector of numbers to format. Standard big.mark = ",".
#' @param ... arguments further passed to the `format` function.
#' 
pNums <- function(x, ...){
  format(x, big.mark = ",", ...)
}

#' Load reference files into the R environment
#' 
#' @usage loadRefFiles(ref, type = "gene.list", freeze = NULL, root = NULL)
#' 
#' @param ref reference object passed from config file.
#' @param type character value of either "gene.list", "GRanges", "data.frame",
#' which determines the output format.
#' @param freeze character value of the reference genome to use when loading
#' the reference file, if the reference needs it. Used if reference is 
#' downloading from UCSC.
#' @param root file path to the root or install directory of iGUIDE.
#' 
#' @description A convenience function that loads several different types of 
#' reference files into the R environment. Allows for input flexibility in 
#' config files when specifying reference datasets.
#' 
loadRefFiles <- function(ref, type = "gene.list", freeze = NULL, root = NULL){
  
  stopifnot(type %in% c("gene.list", "GRanges", "data.frame"))
  
  if( !is.null(root) ) stopifnot(dir.exists(root))

  if( grepl(".rds$", ref$file) ){
    
    if( file.exists(file.path(root, ref$file)) ){
      ref_set <- readRDS(file.path(root, ref$file))
    }else if( file.exists(ref$file) ){
      ref_set <- readRDS(ref$file)
    }else{
      stop("\n  Cannot find file: ", ref$file, ".\n")
    }
    
  }else if( grepl(".RData$", ref$file) ){
    
    if( file.exists(file.path(root, ref$file)) ){
      rdata_path <- file.path(root, ref$file)
    }else if( file.exists(ref$file) ){
      rdata_path <- ref$file
    }else{
      stop("\n  Cannot find file: ", ref$file, ".\n")
    }
    
    ref_env <- new.env()
    load(rdata_path, envir = refs)
    ref_set <- ref_env[[ls(ref_env)]]
    
  }else if( 
    file.exists(file.path(root, ref$file)) | 
      file.exists(ref$file) | 
        grepl("^http", ref$file) 
    ){
    
    if( file.exists(file.path(root, ref$file)) ){
      ref_set <- data.table::fread(
        file.path(root, ref$file), data.table = FALSE
      )
    }else{
      ref_set <- data.table::fread(ref$file, data.table = FALSE)
    }
    
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

#' Convert matches to GRanges objects
#' 
#' @usage mindexToGranges(mindex, strand, ref = NULL)
#' 
#' @param mindex an MIndex object, such as produced by 
#' Biostrings::vmatchPattern.
#' 
#' @param strand character to specify the strand of the alignment, either "+" or
#' "-" or "*".
#' 
#' @param ref reference genome object, sequence information will be included 
#' with output GRanges object if this parameter is specified.
#' 
#' @description Converts MIndex objects into GRange objects. Helpful conversion
#' if using the Biostrings::vmatchPattern function and genome objects.
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

#' Generate random sites across the reference genome.
#' 
#' @usage selectRandomSites(
#'   num, ref.genome, drop.extra.seqs = TRUE, seq.names = NULL, rnd.seed = NULL
#' )
#' 
#' @param num integer the number of random sites to choose.
#' 
#' @param ref.genome BSgenome object of the reference genome.
#' 
#' @param drop.extra.seqs logical, if TRUE then non-standard sequences will be
#' dropped from the possibilities of selecting random sites. FALSE will include
#' the non-standard sequences. Standard sequences are defined by chr1-22, X, Y,
#' and M.
#' @param seq.names character vector of sequence names to select sites from.
#' 
#' @param rnd.seed integer value that functions as the random seed for the 
#' selection process.
#' 
#' @description Generate any number of uniformily distributed random sites 
#' across a reference genome. Output in a GRanges object with the sequence info
#' of the reference genome.
#' 
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

#' Get sequences from a reference genome at a specific site.
#' 
#' @usage getSiteSeqs(gr, upstream.flank, downstream.flank, ref.genome)
#' 
#' @param gr GRanges object or character vector of positions following the 
#' pattern "seqname:strand:position". These objects specify the specific 
#' location in the reference genome and the strand, sense ("+") or anti-sense 
#' ("-"). Only the anchored positions of the ranges will be used. If gr is a 
#' GRanges object, it will pass through `flank(gr, -1, start = TRUE)` to
#' identify these positions.
#' 
#' @param upstream.flank integer value for upstream sequence to return.
#' 
#' @param downstream.flank integer value for downstream sequence to return.
#' 
#' @param ref.genome BSgenome object of the reference genome.
#' 
#' @description Given a list of positions or GRanges object of genomic
#' locations, return the sequences around the positions from a reference genome.
#' 
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
        downstream.flank + 1, 
        -downstream.flank - 1
      )
    ),
    width = upstream.flank + downstream.flank + 1,
    start = TRUE
  ))
  
  uniq_ranges <- unique(seq_ranges)
  
  seqs <- BSgenome::getSeq(ref.genome, uniq_ranges)
  seqs[match(seq_ranges, uniq_ranges)]
  
}

#' Compare target sequences with those found flanking incorporation sites
#' 
#' @usage compareTargetSeqs(
#'   gr.with.sequences, seq.col, target.seqs, tolerance = 6L, 
#'   nuc.profile = NULL, submat = NULL, upstream.flank, downstream.flank
#' )
#' 
#' @param gr.with.sequences GRanges object with a metadata column, specified by
#' `seq.col`, that contains sequences for the range.
#'  
#' @param seq.col character string indicating the metadata column containing 
#' sequence information for the range.
#' 
#' @param target.seqs a list of target sequences to compare against the ranges.
#' 
#' @param tolerance integer specifying the number of mismatches tolerated for a
#' target sequence alignment.
#' 
#' @param nuc.profile a listed object with the following named objects: "PAM" is
#' a character string specifying the PAM sequence (ambiguous nucleotides 
#' allowed, FALSE if no PAM sequence), "PAM_Loc" or PAM location with respect to
#' the target sequences ('5p', '3p', or FALSE), "PAM_Tol" is a positive integer
#' value specifying the allowed mismatch in the PAM sequence alignment (will be
#' ignored if PAM is FALSE), "Editing" is a character string indicating if 
#' editing by the nuclease is "upstream", "downstream", or "internal" with 
#' repect to the PAM sequence (or target sequence if PAM is FALSE), "Cut_Offset"
#' is an integer value specifying the nucleotide from the 5' nucleotide of the 
#' PAM sequence (or target sequence if PAM is FALSE) to call the predicted edit 
#' site (positive indicates downstream and negative indicates upstream, can also
#' specify "mid_insert" to indicate the middle of two paired alignments), 
#' "Insert_size" may be an integer value (15) or a string that will be converted
#' into a range ("15:20") indicating the distance between two paired alignments 
#' for editing (relates to TALEN and other dual-nickase systems). 
#' 
#' @param submat matrix that functions as the substitution matrix for sequence
#' alignment scores. Typically binary to calculate exact mismatches.
#' 
#' 
#' @param upstream.flank,downstream.flank integer value specifying the distance 
#' that was previously used to capture sequences flanking the incorporation 
#' sites.
#' 
#' @description Scan through sequences flanking incorporation sites for matches 
#' to target sequences and PAM motifs (if appliciable). This will create an all 
#' by all comparison and return a GRanges object with the input information and 
#' additional metadata cols with matching information.
#' 
compareTargetSeqs <- function(gr.with.sequences, seq.col, 
                              target.seqs, tolerance = 6L,
                              nuc.profile = NULL, submat = NULL,
                              upstream.flank, downstream.flank
                              ){ 
  
  if( is.null(submat) ) submat <- banmat()
  
  if( is.null(nuc.profile) ){
    
    nuc.profile <- list(
      "PAM" = FALSE, "PAM_Loc" = FALSE, "PAM_Tol" = FALSE,
      "Cut_Offset" = 0, "Insert_size" = FALSE
    )
    
  }
  
  sites <- gr.with.sequences
  sites$siteID <- seq_len(length(sites))
  seq_col_match <- match(seq.col, names(GenomicRanges::mcols(sites)))
  
  if( length(seq_col_match) == 0 ){
    
    stop("Cannot find sequences in column.")
    
  }else{
    
    seqs <- GenomicRanges::mcols(sites)[,seq_col_match]
    names(seqs) <- sites$siteID
    
  }
  
  fwd_df <- alnTargetSeqs(seqs, target.seqs, tolerance)
  
  rev_df <- alnTargetSeqs(
    Biostrings::reverseComplement(seqs), target.seqs, tolerance
  )

  if( nuc.profile$PAM != FALSE ){
    if( nuc.profile$PAM_Loc == "5p" ){
      
      fwd_df$start <- fwd_df$start - nchar(nuc.profile$PAM)
      rev_df$start <- rev_df$start - nchar(nuc.profile$PAM)
      
    }else if( nuc.profile$PAM_Loc == "3p" ){
      
      fwd_df$end <- fwd_df$end + nchar(nuc.profile$PAM)
      rev_df$end <- rev_df$end + nchar(nuc.profile$PAM)
      
    }else{
      
      stop(
        "\n  Config input error in Nuclease profile:\n",
        "    'PAM_Loc' must be [5p]rime or [3p]rime of the target ",
        "sequence if 'PAM' is not FALSE.\n",
        "    Acceptable inputs: '5p', '3p', FALSE."
      )
      
    }
  }
  
  fwd_df <- dplyr::filter(
    fwd_df, start >= 1, end <= upstream.flank + downstream.flank + 1
  )
  
  fwd_df$aln.seq <- as.character(Biostrings::DNAStringSet(
    x = seqs[fwd_df$names], 
    start = fwd_df$start, 
    end = fwd_df$end
  ))
  
  rev_df <- dplyr::filter(
    rev_df, start >= 1, end <= upstream.flank + downstream.flank + 1
  )
  
  rev_df$aln.seq <- as.character(Biostrings::DNAStringSet(
    x = Biostrings::reverseComplement(seqs[rev_df$names]), 
    start = rev_df$start, 
    end = rev_df$end
  ))
  
  fwd_df$target <- paste0(fwd_df$target, rep(":(sense)", nrow(fwd_df)))
  rev_df$target <- paste0(rev_df$target, rep(":(antisense)", nrow(rev_df)))
  
  matched_seqs <- rbind(fwd_df, rev_df)
  
  # Filter by PAM match if PAM present in nuclease profile.
  if( nuc.profile$PAM != FALSE ){
    
    pam_matched <- alnTargetSeqs(
      seqs = matched_seqs$aln.seq, 
      target.seqs = nuc.profile$PAM, 
      tolerance = nuc.profile$PAM_Tol, 
      fixed = FALSE
    )
    
    if( nuc.profile$PAM_Loc == "5p" ){
      pam_matched <- dplyr::filter(pam_matched, start == 1)
    }else if( nuc.profile$PAM_Loc == "3p" ){
      pam_matched <- dplyr::filter(pam_matched, end == nt_width)
    }
    
    matched_seqs <- matched_seqs[ unique(as.numeric(pam_matched$names)), ]
    
  }
  
  good_alns <- as.numeric(matched_seqs$names)
  
  if( length(good_alns) == 0 ){
    
    sites$target.match <- "No_valid_match"
    sites$target.mismatch <- NA
    sites$target.score <- NA
    sites$aligned.sequence <- NA
    sites$edit.site <- NA
    return(sites)
    
  }
  
  if( length(good_alns) != length(sites) ){
    
    non_probable_sites <- sites[!sites$siteID %in% good_alns]
    non_probable_sites$target.match <- "No_valid_match"
    non_probable_sites$target.mismatch <- NA
    non_probable_sites$aligned.sequence <- NA
    non_probable_sites$edit.site <- NA
    
  }
  
  potential_sites <- sites[good_alns]
  
  potential_sites$target.match <- matched_seqs$target
  potential_sites$target.mismatch <- matched_seqs$target.mismatch
  potential_sites$aligned.sequence <- matched_seqs$aln.seq
  
  potential_sites$edit.site <- calcCutSite(
    potential_sites, matched_seqs, upstream.flank, 
    downstream.flank, nuc.profile
  )
  
  if( length(good_alns) != length(sites) ){
    all_sites <- c(potential_sites, non_probable_sites)
  }else{
    all_sites <- potential_sites
  }
  
  all_sites <- all_sites[order(all_sites$siteID)]
  
  if( any(duplicated(all_sites$siteID)) ){
    
    dup_ids <- names(table(all_sites$siteID))[table(all_sites$siteID) > 1]
    dup_sites <- all_sites[all_sites$siteID %in% dup_ids]
    
    adj_sites <- unlist(GenomicRanges::GRangesList( 
      lapply(
        split(dup_sites, dup_sites$siteID), 
        function(gr){
          
          if( length(unique(gr$target.mismatch)) > 1){
            gr <- gr[gr$target.mismatch == min(gr$target.mismatch)]
          }
          
          if( length(gr) > 1){
            dists <- abs(start(gr) - 
              as.numeric(stringr::str_extract(gr$edit.site, "[\\d]+$")))
            gr <- gr[dists == min(dists)]
          }
          
          gr
          
        }
      )
    ))
    
    all_sites <- c(all_sites[!all_sites$siteID %in% dup_ids], adj_sites)
    all_sites <- all_sites[order(all_sites$siteID)]
    
  }
  
  if( any(duplicated(all_sites$siteID)) ){
    cat(
      "  Note: Some alignments with multiple target matches. Total number of", 
      "alignments expanded. See siteID.\n"
    )
  }else{  
    all_sites$siteID <- NULL
  }
  
  all_sites
  
}

#' Align target sequences against flanking sequences around incorporation sites
#' 
#' @usage alnTargetSeqs(seqs, target.seqs, tolerance)
#' 
#' @param seqs sequences in either Biostrings object or character vector.
#' 
#' @param target.seqs a list of target sequences to compare against the ranges.
#' 
#' @param tolerance integer specifying the number of mismatches tolerated for an
#' alignment.
#' 
#' @param fixed string following the behavior of 
#' Biostrings::vmatchPattern(fixed = ...). Options include logical where TRUE
#' indicates that ambiguous nucleotides can only match the same ambiguous nts,
#' FALSE refers to allowing for IUPAC ambiguity matching, 'subject' and 
#' 'pattern' can be used to fix one or the other where subject is the 'seqs' 
#' object and 'pattern' is the 'target.seqs' object.
#' 
#' @description Companion function to `compareTargetSeqs` which aligns target 
#' sequences to input sequences to determine potential locations of matched 
#' alignments.
#' 
alnTargetSeqs <- function(seqs, target.seqs, tolerance, fixed = 'subject'){
  
  if( is.null(names(seqs)) ) names(seqs) <- seq_along(seqs)
  if( class(seqs) != "DNAStringSet") seqs <- Biostrings::DNAStringSet(seqs)
  if( is.null(names(target.seqs)) ) names(target.seqs) <- target.seqs
  
  nt_widths <- data.frame(
    names = names(seqs),
    nt_width = Biostrings::width(seqs)
  )
  
  alns <- lapply(
    0:tolerance, 
    function(tol, targets, seqs){
      
      lapply(
        targets, 
        function(target){
          
          unlist(Biostrings::vmatchPattern(
            pattern = target, 
            subject = seqs, 
            max.mismatch = tol, 
            fixed = fixed
          ))
          
        }
      )
      
    }, 
    targets = target.seqs, 
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
              d$target <- rep(names(alns[[i]][j]), nrow(d))
              d$mismatches <- rep(i-1, nrow(d))
              d
              
            }
          )
        )
        
      }
    )
  )
  
  dplyr::group_by(df, start, end, width, names, target) %>%
    dplyr::mutate(target.mismatch = min(mismatches)) %>%
    dplyr::ungroup() %>%
    dplyr::select(names, target, target.mismatch, start, end, width) %>%
    dplyr::distinct() %>%
    dplyr::left_join(., nt_widths, by = "names") %>%
    dplyr::mutate(
      start = ifelse(start <= 0, 1, start),
      end = ifelse(end > nt_width, nt_width, end),
      width = end - start + 1
    )
  
}

#' Determine the editing site location
#' 
#' @usage calcCutSite(
#'   sites, matched.seqs, upstream.flank, downstream.flank, PAM, offset.nt
#' )
#' 
#' @param sites GRanges object representing the incorporation sites
#'  
#' @param matched.seqs character vector of sequences matching gRNAs, 
#' corresponding to the `sites`.
#' 
#' @param upstream.flank,downstream.flank integer value specifying the distance 
#' that was previously used to capture sequences flanking the incorporation 
#' sites.
#' 
#' @param PAM character string indicating the patter to recognize the PAM
#' sequence.
#' 
#' @param offset.nt integer value specifying the base in the gRNA sequence where
#' editing occurs most often. Numbered bases from the 3' end of the gRNA
#' sequences, not including the PAM.
#' 
#' @description Companion function to `compareGuideRNAs` which determines 
#' expected locations of nuclease editing based on matched sequences to the gRNA
#' sequences and known upstream and downstream flanking distances.
#' 
calcCutSite <- function(sites, matched.seqs, upstream.flank, 
                         downstream.flank, nuc.profile){
  
  # Include PAM if present
  if( nuc.profile$PAM != FALSE ){
    
    PAM <- nuc.profile$PAM
    pam_present <- TRUE
    pam_loc <- nuc.profile$PAM_Loc
    
  }else{
    
    PAM <- ""
    pam_present <- FALSE
    pam_loc <- FALSE
    
  }
  
  # Determine offset from PAM or target sequence
  if( is.numeric(nuc.profile$Cut_Offset) ){
    
    offset_nt <- nuc.profile$Cut_Offset
    
  }else if( nuc.profile$Cut_Offset == "mid_insert" ){
    
    insert_size <- nuc.profile$Insert_size
    
    if( grepl("\\:", insert_size) & length(insert_size) == 1 ){
      
      insert_size <- seq(
        unlist(strsplit(insert_size, ":"))[1],
        unlist(strsplit(insert_size, ":"))[2]
      )
      
    }
    
    offset_nt <- unname(round(quantile(insert_size, probs = 0.5)))
    
  }else{
    
    stop(
      "\n  Config input error in Nuclease profile:\n",
      "    'Cut_Offset' must be an integer value specifying the editing\n",
      "      distance in nucleotides from the 5' end of the PAM sequence\n", 
      "      or target sequence (if no PAM is present).\n",
      "      Acceptable inputs: positive and negative integer values, and\n",
      "      'mid_insert'.\n",
      "    'Insert_Size' should also be an integer value or a range,\n", 
      "      indicated by 'min:max'. FALSE can be provided if not \n", 
      "      applicable. Refer to documentation for more information.\n"
    )
    
  }
  
  # Combine alignment information and determine predicted cutsites
  df <- data.frame(
      "chr" = GenomicRanges::seqnames(sites),
      "strand" = GenomicRanges::strand(sites),
      "pos" = GenomicRanges::start(
        GenomicRanges::flank(sites, -1, start = TRUE)
      ),
      "target" = matched.seqs$target,
      "target.aln" = ifelse(
        !grepl(":(antisense)", matched.seqs$target, fixed = TRUE), "+", "-"
      ),
      "target.start" = matched.seqs$start,
      "target.end" = matched.seqs$end
    ) %>%
    dplyr::mutate(
      true.ort = ifelse(
        as.character(strand) == as.character(target.aln), "+", "-"
      ),
      flank.start = ifelse(
        strand == "+", pos - upstream.flank, pos - downstream.flank + 1
      ),
      flank.end = ifelse(
        strand == "+", pos + downstream.flank, pos + upstream.flank - 1
      ),
      tar.start = ifelse(
        strand == "+", 
        ifelse(
          target.aln == "+", 
          flank.start + target.start - 1,
          flank.end - target.end + 1
        ),
        ifelse(
          target.aln == "+",
          flank.end - target.end + 2,
          flank.start + target.start - 2
        )
      ),
      tar.end = tar.start + (target.end - target.start)
    )
  
  if( pam_present ){
    if( pam_loc == "5p" ){
      
      df <- df %>%
        dplyr::mutate(
          pam.start = ifelse(
            true.ort == "+", tar.start, tar.end - nchar(PAM) + 1
          ),
          pam.end = ifelse(
            true.ort == "+", tar.start + nchar(PAM) - 1, tar.end
          ),
          edit.pos = ifelse(
            true.ort == "+", pam.start + offset_nt, pam.end - offset_nt
          )
        )
      
    }else if( pam_loc == "3p" ){
      
      df <- df %>%
        dplyr::mutate(
          pam.start = ifelse(
            true.ort == "+", tar.end - nchar(PAM) + 1, tar.start
          ),
          pam.end = ifelse(
            true.ort == "+", tar.end, tar.start + nchar(PAM) - 1
          ),
          edit.pos = ifelse(
            true.ort == "+", pam.start + offset_nt, pam.end - offset_nt
          )
        )
      
    }else{
      
      stop(
        "\n  Config input error in Nuclease profile:\n",
        "    'PAM_Loc' must be [5p]rime or [3p]rime of the target ",
        "sequence if 'PAM' is not FALSE.\n",
        "    Acceptable inputs: '5p', '3p', FALSE."
      )
      
    }
    
  }else{
    
    df <- df %>%
      dplyr::mutate(
        edit.pos = ifelse(
          true.ort == "+", tar.start + offset_nt, tar.end - offset_nt
        )
      )
    
  }
  
  return(paste0(df$chr, ":", df$true.ort, ":", df$edit.pos))
  
}

#' Expand position range notation to position identifier strings
#' 
#' @usage expandPosStr(posStr, delim = ":")
#' 
#' @param pos.str character vector of position strings, in the format of seqname,
#' orientation, and position. For position, if two numbers are given, then they 
#' specify a range that will be expanded. For orientation, "+" and "-" indicates
#' the possible orientation, but if "*" is used, then the posStr will be 
#' expanded to "+" and "-".
#' 
#' @param delim character string indicating the symbol delimiting the 
#' information within `posStr`. A range should always be delimited by a "-" 
#' between two numbers.
#' 
#' @param return character string of either "vector" (default) or "list", 
#' indicating the type of output that is returned from the function. The vector
#' option will return a single vector of all possible position strings for the 
#' given input, while the list option will return a listed object where each 
#' object in the list is related to the expansion of the input. If the 
#' input object is named, then the vecotr / list names will be the same as the 
#' input names, otherwise the vector / list will be unnamed. 
#' 
#' @description Given position strings and string ranges, for example 
#' "chr2:+:392" and "chr2:+:392-394", expand to all possible position strings
#' of the format seqnames, orientation ("+"/"-"), and position delimited by the
#' symbol included as `delim`. Orientation will also be expanded. For example,
#' "chr2:+:392-394" will be expanded to "chr2:+:392", "chr2:+:393", and 
#' "chr2:+:394". Additionally, "chr2:*:392" will be expanded to "chr2:+:392" and
#' "chr2:-:392".
#' 
expandPosStr <- function(pos.str, delim = ":", return = "vector"){
  
  # Check
  if( !return %in% c("vector", "list") ){
    stop("\n  Output from expandPosStr must be either 'vector' or 'list'.\n")
  }
  
  # Split input into a matrix
  pos_mat <- stringr::str_split(
    pos.str, pattern = delim, n = 3, simplify = TRUE
  )
  
  # Expand matrix based on orientation
  exp_ort <- lapply(pos_mat[, 2, drop = TRUE], function(x){
    if( x %in% c("+", "-")){
      return(x)
    }else{
      return(c("+", "-"))
    }
  })
  
  pos_mod <- pos_mat[rep(seq_along(pos.str), lengths(exp_ort)), , drop = FALSE]
  pos_mod[,2] <- unlist(exp_ort)
  
  # Expand the matrix based on position ranges
  exp_pos <- stringr::str_split(
    pos_mod[, 3, drop = TRUE], pattern = "-", simplify = FALSE
  )
  
  exp_pos <- lapply(exp_pos, function(x){
    x <- as.numeric(x)
    if( length(x) == 1 ) x <- c(x, x)
    seq(from = x[1], to = x[2])
  })
  
  pos_mod <- pos_mod[
    rep(seq_len(nrow(pos_mod)), lengths(exp_pos)), , drop = FALSE
  ]
  
  pos_mod[,3] <- unlist(exp_pos)
  
  # Expand names if needed and return character vector of position strings
  if( !is.null(names(pos.str)) ){
    
    nam <- names(pos.str)[rep(seq_along(pos.str), lengths(exp_ort))]
    nam <- nam[rep(seq_along(nam), lengths(exp_pos))]
    
    if( return == "vector" ){
      return(structure(vcollapse(pos_mod, sep = delim), names = nam))
    }else if( return == "list" ){
      return(split(vcollapse(pos_mod, sep = delim), nam))
    }
    
  }else{
    
    idx <- rep(seq_along(pos.str), lengths(exp_ort))
    idx <- idx[rep(seq_along(idx), lengths(exp_pos))]
    
    if( return == "vector" ){
      return(vcollapse(pos_mod, sep = delim))
    }else if( return == "list" ){
      return(unname(split(vcollapse(pos_mod, sep = delim), idx)))
    }
    
  }
  
}

#' Filter out inappropriate comparisons based on experimental conditions
#' 
#' @usage filterInappropriateComparisons(guideRNA.match, specimen, treatment)
#' 
#' @param guideRNA.match character vector of gRNA names, indicating a gRNA that
#' a single incorporation is matched too. Needs to be the same length as 
#' `specimen`.
#' 
#' @param specimen character vector of specimen names. Needs to be the same 
#' length as `guideRNA.match`.
#' 
#' @param treatment list containing character vectors matching gRNA names and 
#' the names of each component of the list are the specimen names.
#' 
#' @description After comparing gRNA sequences to regions identified by 
#' incorporations, some comparisons may not be valid depending on how the
#' samples were prepared. For example, if sample A was treated with gRNA1 and 
#' sample B was treated with gRNA2, then it makes little sense to assume
#' incorporations in sample B with sequence related to gRNA1 are part of the 
#' nuclease specific activity. This function removes any inappropriately matched
#' annotations given speciment and treatment data.
#' 
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

#' Assign a Gene ID to a given position
#' 
#' @usage assignGeneID(
#'   seqnames, positions, reference, ref.genes, onco.genes, special.genes, 
#'   annotations = TRUE
#' )
#' 
#' @param seqnames character vector of seqnames present in `reference`.
#' 
#' @param positions integer vector of positions on the seqnames to assign the
#' gene ID, must be the same length as `seqnames`.
#' 
#' @param reference BSgenome object of the reference genome.
#' 
#' @param ref.genes GRanges object with transcription unit references. Names 
#' used in the assignment should be in a metadata column named "annot_sym".
#' 
#' @param onco.genes A list of gene names, present in the `ref.genes` GRanges
#' object, that indicate cancer-associated genes.
#' 
#' @param special.genes A list gene names for special annotation, should be 
#' present in the `ref.genes` GRanges object.
#' 
#' @param annotations logical indicating if annotations should be included in 
#' the output character vector.
#' 
#' @description For a given genomic location, assign it a gene ID which is the 
#' name of the nearest, or within, transcription unit. Additionally, annotations
#' are provided that specify the following: "*" means the position is within the
#' transcription unit for the specified gene, "~" means that the gene is on the
#' list provided by `onco.genes`, and lastly "!" means the gene is present on
#' the `special.genes` list.
#' 
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

#' Calculate the coverage of given ranges
#' 
#' @usage calcCoverage(gr, resolution)
#' 
#' @param gr GRanges object.
#' 
#' @param resolution integer value specifying the width for which to caclulate
#' coverage. 1L will determine the coverage at every position within the given
#' range, while 10L will determine the coverage in 10 bp chuncks.
#' 
#' @description Convert an input GRanges object into a coverage GRanges object
#' which indicates the number of alignments overlapping at different positions.
#' This is commonly known as the coverage, or how many times did you observe a 
#' region of the genome given GRange information. This is a companion function
#' to `plotCoverage` which does the work of generating the data object to plot.
#' 
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

#' Coverage plot
#' 
#' @usage plotCoverage(gr, resolution = 10L)
#' 
#' @param gr GRanges object.
#' 
#' @param resolution integer value specifying the width for which to caclulate
#' coverage. 1L will determine the coverage at every position within the given
#' range, while 10L will determine the coverage in 10 bp chuncks.
#' 
#' @description Plot the amount of coverage given a set of alignment ranges.
#' This function will convert the input GRanges object into a coverage object 
#' using the `calcCoverage` function and then plot it using a specific format.
#' 
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

#' Plot edit site coverage
#' 
#' @usage plotEditSites(gr, sampleName = NULL, resolution = 10L)
#' 
#' @param gr GRanges object with metadata column 'edit.site' which functions as
#' a grouping vector for the output plots. Additionally, the column 
#' "guideRNA.match" will be used in the title data.
#' 
#' @param sampleName character string specifying the name of metadata column 
#' with sample designation.
#' 
#' @param resolution integer value specifying the width for which to caclulate
#' coverage. 1L will determine the coverage at every position within the given
#' range, while 10L will determine the coverage in 10 bp chuncks.
#' 
#' @description Another way to plot coverage is to focus on called edit sites
#' and generate an individual plot for each. This plot will take a GRanges 
#' object with a metadata column 'edit.site' and generate a coverage plot for
#' each specified unique edit.site. This is a visual way to scan over the 
#' results and determine if edit sites that are being called really appear to be
#' nuclease dependent sites.
#' 
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

#' Modify tile-based plots to have square tiles
#' 
#' @usage makeSquare(p, dims, fudge=1)
#' 
#' @param p ggplot object using geom tile or that you would like to fix the 
#' aspect ratio
#' 
#' @param dims dimenstion object of the input matrix to a heatmap. A list of
#' integer vectors of length 2, with the names "ncols" and "nrows".
#' 
#' @description Use to fix an aspect ratio of a tile-based plot to be square.
#' 
makeSquare <- function(p, dims, fudge=1){
  
  dims <- heatmap_dims(p)
  p + ggplot2::theme(aspect.ratio = (dims$nrows / dims$ncols) * fudge)
  
}

#' Combine a list of ShortRead objects
#' 
#' @param split.seqs list of ShortRead objects
#' 
#' @description Utility function that combines ShortRead objects.
#' 
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
#' 
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
#' @usage predictESProb(x, density, range = NULL)
#' 
#' @param x integer position within range of density object indicating the 
#' distance from the predicted edit site.
#' 
#' @param density a density object constructed from the incorporation site 
#' distribution around associated On-target editing site(s).
#' 
#' @param range a numeric vector indicating the range of data to consider. NULL
#' will default to the range of input for the denisity object.
#' 
#' @description 
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

#' Generate a GRanges object that covers the entire genome at specified 
#' resolution
#' 
#' @usage generateGenomicRegions(ref, res, drop.alt.chr = TRUE)
#' 
#' @param ref reference genome as a BSgenome object.
#' 
#' @param res integer value to specify the width of each range. The number of
#' ranges generated will depend on the resolution to cover the entire genome.
#' 
#' @param drop.alt.chr logical specifying if alternative chromosomes 
#' (seqnames != chr1:22, X, Y, M) should be removed from the output GRanges
#' object.
#' 
#' @description Given a referene BSgenome object, generate a GRanges object with
#' ranges that span the entire genome with `res` resolution. This can then be
#' used as a scanning reference range for the density or counts of hits within 
#' the regions.
#' 
generateGenomicRegions <- function(ref, res, drop.alt.chr = TRUE){
  
  if( class(ref) == "BSgenome" ) ref <- GenomicRanges::seqinfo(ref)
  
  if( drop.alt.chr ){
    ref <- ref[names(ref)[which(!stringr::str_detect(names(ref), "_"))]]
  }
  
  ref_len <- GenomeInfoDb::seqlengths(ref)
  
  unlist(GenomicRanges::GRangesList(lapply(
    seq_along(ref_len), 
    function(i){
      
      seq <- names(ref_len)[i]
      starts <- seq(1, ref_len[seq], res)
      
      if( res > ref_len[seq] ){
        ends <- ref_len[seq]
      }else{
        ends <- seq(res, ref_len[seq], res)
      }
      
      if( length(starts) - length(ends) == 1 ) ends <- c(ends, ref_len[seq])
      
      GenomicRanges::GRanges(
        seqnames = seq, 
        ranges = IRanges::IRanges(start = starts, end = ends), 
        strand = "*",
        seqinfo = ref)
      
    }
  )))
  
}

#' Determine the genomic density of input GRanges across the reference genome.
#' 
#' @usage genomicDensity(gr, res, cutoff = 2, adj = 1, drop.alt.chr = TRUE)
#' 
#' @param gr GRanges object which will be used to calculate the genomic density.
#' 
#' @param res integer value to specify the width of each range. The number of
#' ranges generated will depend on the resolution to cover the entire reference
#' genome.
#' 
#' @param cutoff integer value specifying minimum counts of input ranges in `gr`
#' required in an output region of the genome to be included. Any regions with
#' counts below this cutoff will be dropped from the output.
#' 
#' @param adj integer value to add to every region as an adjustment factor. If 
#' using the log.count or norm.log.count outputs, adding 1 will still include 
#' all ranges rather than lead to -Inf values.
#' 
#' @param drop.alt.chr logical specifying if alternative chromosomes 
#' (seqnames != chr1:22, X, Y, M) should be removed from the output GRanges
#' object.
#' 
#' @description This function can be used to determine the genomic density 
#' of the input GRanges object. The output GRanges object will contain several
#' columns with count information in the metadata columns. The columns include:
#' "count" - the number of ranges within the region, "log.count" - the log 
#' transformation of the count plus the adjustment, and "norm.log.count" - where
#' the "log.count" has been normalized such that the highest value is 1.0. 
#' 
genomicDensity <- function(gr, res, cutoff = 2, adj = 1, drop.alt.chr = TRUE){
  
  stopifnot(class(gr) == "GRanges")
  
  if(
    is.na(sum(is.numeric(GenomeInfoDb::seqlengths(GenomicRanges::seqinfo(gr)))))
  ){
    stop("SeqInfo should be present for input GRanges.")
  }
  
  ref_regions <- generateGenomicRegions(
    GenomicRanges::seqinfo(gr), res, drop.alt.chr = drop.alt.chr
  )
  
  ref_regions$count <- GenomicRanges::countOverlaps(ref_regions, gr)
  ref_regions$log.count <- log(ref_regions$count + adj, base = 10)
  ref_regions$norm.log.count <- ref_regions$log.count/max(ref_regions$log.count)
  ref_regions <- ref_regions[ref_regions$count >= cutoff]
  ref_regions
  
}

#' Collapse row contents of a data.frame or matrix into single vector.
#'
#' \code{vcollapse} returns a single vector from input data.frame or matrix 
#' where row contents have been combined or collapsed, as with 
#' `paste(..., collaspe = "")`.
#'
#' @description Similar to python zip, `vzip` takes input vectors and merges
#' them together by their input order and index. A simple example is two numeric
#' vectors, A = c(1,1,1) and B = c(2,2,2). The output of vzip(A,B) would simply
#' be a single vector of c(1,2,1,2,1,2). Any number of vectors can be input, but
#' each input vector must be of the same length. Output vector class depends on
#' input vector consensus.
#'
#' @usage
#' vcollapse(d)
#' vcollapse(d, sep = "-", fill = "NA")
#'
#' @param d data.frame or matrix or object coercible to a matrix. Row contents
#' will be combined into a single output vector.
#' 
#' @param sep character used to separate collapsed contents.
#' 
#' @param fill character used to fill empty values within the coerced object.
#'
#' @examples
#' df <- data.frame(
#'   "A" = letters[1:5],
#'   "B" = 3:7,
#'   "C" = LETTERS[2:6])
#' vcollapse(df)
#' vcollapse(df, sep = "-")
#' vcollapse(df, sep = "-", fill = "z")
#'
vcollapse <- function(d, sep, fill = "NA"){
  
  if( is.vector(d) ){
    stop(
      "Function vcollapse() is not used on vectors, use paste(collapse = ...)."
    )
  }
  
  if( any(sapply(seq_len(ncol(d)), function(i) class(d[,i])) == "factor") ){
    
    fct_idx <- which(
      sapply(seq_len(ncol(d)), function(i) class(d[,i])) == "factor"
    )
    
    mod_env <- new.env()
    mod_env$d <- d
    
    null <- lapply(
      fct_idx, 
      function(i){
        mod_env$d[,i] <- as.character(d[,i])
      }
    )
    
    d <- mod_env$d
    
  }
  
  if( class(d) != "matrix" ) d <- as.matrix(d)
  
  mat <- d
  
  if( !is.null(fill) ) mat <- ifelse(is.na(mat), fill, mat)
  
  mat <- do.call(
    cbind, 
    lapply(
      seq_len(ncol(mat)), 
      function(i){
        
        if( i < ncol(mat) ){
          cbind(mat[,i], rep(sep, nrow(mat)))
        }else{
          mat[,i]
        } 
        
      }
    )
  )
  
  mat <- cbind(mat, rep(">|<", nrow(mat)))
  div_str <- stringr::str_c(t(mat), collapse = "")
  unlist(strsplit(div_str, split = ">\\|<"))
  
}

#' Generate a genomic density plot
#' 
#' @usage plotGenomicDensity(
#'   grl, res = 1E7, grp.col = NULL, cutoff = 2, drop.alt.chr = TRUE, 
#'   clean = FALSE
#' )
#' 
#' @param grl GRangesList object of length 4, representing all alignments, 
#' pileup alignments, paired alignments, and gRNA matched alignments.
#' 
#' @param res integer value specifying the resolution on the genome to display
#' the genomic density, "norm.log.count", see `genomicDensity`.
#' 
#' @param grp.col character string specifying the metadata column name in each
#' component of the list with grouping information. 
#' 
#' @param cutoff integer value specifying minimum counts of input ranges in `gr`
#' required in an output region of the genome to be included. Any regions with
#' counts below this cutoff will be dropped from the output.
#' 
#' @param abund.method character string specifying the abundance to plot. Types
#' include "count", "log.count", "norm.log.count". Refer to `genomicDensity` for
#' descriptions of the abundances.
#' 
#' @param drop.alt.chr logical specifying if alternative chromosomes 
#' (seqnames != chr1:22, X, Y, M) should be removed from the output GRanges
#' object.
#' 
#' @param clean logical specifying if the plot should be void of axis and legend
#' labels (or clean = TRUE). 
#' 
#' @description A circular plot that puts each chromosome end to end and 
#' displays the normalized log density of each level of alignment. Used to get 
#' a genomic overview / summary of detected double strand breaks.
#' 
plotGenomicDensity <- function(grl, res = 1E7, grp.col = NULL, cutoff = 2, 
                               abund.method = "norm.log.count", 
                               drop.alt.chr = TRUE, clean = FALSE){
  
  if( class(grl) == "GRanges" ) grl <- GenomicRanges::GRangesList(grl)
  
  if( !is.null(grp.col) ){
    
    stopifnot(all(unlist(lapply(
      grl, 
      function(x) grp.col %in% names(GenomicRanges::mcols(x))
    ))))
    
  }
  
  ref_len <- res * ceiling(
    GenomeInfoDb::seqlengths(GenomicRanges::seqinfo(grl)) / res
  )
  
  ref_cum_len <- structure(
    c(0, cumsum(as.numeric(ref_len))), names = c("start", names(ref_len))
  )
  
  gen_densities <- lapply(
    grl, 
    function(x){
      
      if( is.null(grp.col) ){
        x <- list(x)
      }else{
        x <- split(x, GenomicRanges::mcols(x)[,grp.col, drop = FALSE])
      }
      
      dplyr::bind_rows(lapply(
        x, 
        function(y){
          y <- genomicDensity(
            y, res = res, cutoff = cutoff, drop.alt.chr = drop.alt.chr
          )
          
          df <- as.data.frame(y, row.names = NULL) 
          
          cum_adj_pos <- ref_cum_len[
            match(df$seqnames, names(ref_cum_len)) - 1
            ]
          
          df$end <- res * ceiling(df$end / res)
          df$width <- df$end - df$start + 1
          
          df$adj.start <- cum_adj_pos + df$start
          df$adj.end <- cum_adj_pos + df$end
          
          df
          
        }), 
        .id = "cond")
      
    }
  )
  
  gen_den <- dplyr::bind_rows(gen_densities, .id = "grp")
  
  
  if( !is.null(names(grl)) ){
    
    gen_den$type <- factor(
      names(grl)[match(gen_den$grp, names(grl))], 
      levels = names(grl)
    )
    
  }else{
    
    gen_den$type <- factor(" ")
    
  }
  
  
  if( is.factor(GenomicRanges::mcols(grl[[1]])[,grp.col]) ){
    gen_den$cond <- factor(
      gen_den$cond, 
      levels = levels(GenomicRanges::mcols(grl[[1]])[,grp.col])
    )
  }
  
  max_scores <- sapply(
    split(gen_den[,match(abund.method, names(gen_den))], gen_den$type), max
  )
  
  cum_max_scores <- sapply(
    levels(gen_den$type), function(x){
      sum(max_scores[seq_len(match(x, names(max_scores))-1)])
  }) + 1
  
  gen_den$score <- cum_max_scores[gen_den$type] + 
    gen_den[,match(abund.method, names(gen_den))]
  
  # Grid layout
  x_breaks <- ref_cum_len[
    seq_len(max(match(unique(gen_den$seqnames), names(ref_cum_len))))
    ]
  
  x_lab_pos <- structure(
    sapply(
      2:length(x_breaks), 
      function(i) mean(c(x_breaks[i-1], x_breaks[i])) 
    ),
    names = names(x_breaks)[2:length(x_breaks)]
  )
  
  y_breaks <- cum_max_scores
  
  p <- ggplot2::ggplot(gen_den) 
  
  if( !clean ){
    
    p <- p +
      ggplot2::geom_hline(yintercept = y_breaks, color = "grey90") +
      ggplot2::geom_vline(xintercept = x_breaks, color = "grey90") +
      ggplot2::scale_x_continuous(
        breaks = x_lab_pos,
        labels = gsub("chr", "", names(x_lab_pos))
      )
    
  }else{
    
    p <- p + ggplot2::scale_x_continuous(labels = NULL)
    
  }
  
  p <- p + 
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = adj.start, xmax = adj.end, 
        ymin = cum_max_scores[type], ymax = score, 
        fill = type
      )
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, max(gen_den$score)), breaks = cum_max_scores, labels = NULL
    ) +
    ggplot2::scale_fill_brewer(
      type = "qual", palette = "Set1", direction = -1
    ) +
    ggplot2::labs(x = "Chromosome", fill = "Levels") + 
    ggplot2::coord_polar() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(color = "white"),
      panel.border = ggplot2::element_rect(color = "white"),
      panel.grid = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
  
  if( clean ){
    p + ggplot2::theme(
      legend.position = "none", 
      axis.title = ggplot2::element_blank()
    )
  }else{
    p
  }
  
}

#' Remove matching seuqences from input vector and replace with character
#' 
#' @usage divSeq(seqs, ref, match.chr = ".")
#' 
#' @param seqs character vector of sequences to compare against the reference.
#' 
#' @param ref character string of a sequence that functions as the reference.
#' 
#' @param match.chr character used to replace characters in `seqs` that match in
#' the same position in `ref`.
#' 
#' @param fill character string specifying if sequences of length less than the
#' reference should be filled in. Fill should be either "left" (insert blank 
#' spaces on the left) or "right" (for inserting blank spaces on the right). 
#' 
#' @description Given a vector of character strings, compare them to the
#' reference. Any character that matches the reference, in the same position,
#' will be converted into the `match.chr` while any sequence diverging from the
#' reference will be left unchanged. This allows the user to easily identify the
#' differences between sequence strings.
#' 
divSeq <- function(seqs, ref, match.chr = ".", fill = "left"){
  
  seqs <- as.character(seqs)
  ref <- as.character(ref)
  
  if( !all(nchar(as.character(seqs)) == nchar(ref)) ){
    
    fill_idx <- which(nchar(seqs) != nchar(ref))
    fill_width <- nchar(ref) - nchar(seqs)[fill_idx]
    
    if( fill == "left" ){
      
      seqs[fill_idx] <- sapply(
        seq_along(fill_idx), 
        function(i){
          paste0(
            paste(rep("N", fill_width[i]), collapse = ""), seqs[fill_idx[i]])
        }
      )
      
    }else if( fill == "right" ){
      
      seqs[fill_idx] <- sapply(
        seq_along(fill_idx), 
        function(i){
          paste0(
            seqs[fill_idx[i]], paste(rep("N", fill_width[i]), collapse = ""))
        }
      )
      
    }else{
      
      stop("fill parameter must be either left or right.")
      
    }
    
  }
  
  stopifnot( all(nchar(as.character(seqs)) == nchar(ref)) )
  
  seq_mat <- stringr::str_split(seqs, pattern = "", simplify = TRUE)
  ref_splt <- unlist(stringr::str_split(ref, pattern = ""))
  
  div_mat <- sapply(
    seq_along(ref_splt), 
    function(i) ifelse(seq_mat[,i] == ref_splt[i], match.chr, seq_mat[,i])
  )
  
  div_str <- stringr::str_c(t(div_mat), collapse = "")
  
  unlist(strsplit(
    div_str, split = paste0("(?<=.{", nchar(ref), "})"), perl = TRUE
  ))
  
}

#' Sequence diverge plot
#' 
#' @usage plotSeqDiverge(
#'   df, ref, nuc.col = NULL, padding = 4, text.size = 2, convert.seq = TRUE, 
#'   force.sq = FALSE, font.family = "Courier", font.face = "bold", 
#'   fill = "left"
#' )
#' 
#' @param df data.frame containing sequence and additional information to be 
#' displayed.
#' 
#' @param ref character string of a sequence that functions as the reference.
#' 
#' @param nuc.col character string specifying the name of the nucleotide 
#' sequence column to display in the input data.frame. If not provided, the 
#' function will assume it is the first column in the data.frame.
#' 
#' @param padding integer value specifying the space to insert between columns 
#' and the table in the output plot.
#' 
#' @param text.size integer value specifying the text size for the plot.
#' 
#' @param convert.seq logical to convert sequences to only show diverging 
#' sequences from reference or to display whole strings.
#' 
#' @param force.sq logical specifying if the output plot should have a fixed 
#' ratio.
#' 
#' @param font.family character string specifying the font family to use for the
#' plot. Recommended to use a monospace font family.
#' 
#' @param font.face character string specifying the font face to use in the 
#' plot.
#' 
#' @param fill character string specifying if sequences of length less than the
#' reference should be filled in. Fill should be either "left" (insert blank 
#' spaces on the left) or "right" (for inserting blank spaces on the right). 
#' 
#' @description A plot the displays the nucleotide divergence from the reference
#' sequence. This plot is easy to identify differences between the reference 
#' sequence and the input sequences. Additionally, other columns in the input
#' data.frame will be displayed in a table adjacent to the plot.
#' 
plotSeqDiverge <- function(df, ref, nuc.col = NULL, padding = 4, 
                             text.size = 2, convert.seq = TRUE, 
                             force.sq = FALSE, font.family = "Courier",
                             font.face = "bold", fill = "left"){
  
  if( is.null(nuc.col) ) nuc.col <- names(df)[1]
  
  seqs <- dplyr::pull(df, var = match(nuc.col, names(df)))
  
  if( !all(nchar(as.character(seqs)) == nchar(ref)) ){
    
    fill_idx <- which(nchar(seqs) != nchar(ref))
    fill_width <- nchar(ref) - nchar(seqs)[fill_idx]
    
    if( fill == "left" ){
      
      seqs[fill_idx] <- sapply(
        seq_along(fill_idx), 
        function(i){
          paste0(
            paste(rep("N", fill_width[i]), collapse = ""), seqs[fill_idx[i]])
        }
      )
      
    }else if( fill == "right" ){
      
      seqs[fill_idx] <- sapply(
        seq_along(fill_idx), 
        function(i){
          paste0(
            seqs[fill_idx[i]], paste(rep("N", fill_width[i]), collapse = ""))
        }
      )
      
    }else{
      
      stop("fill parameter must be either left or right.")
      
    }
    
  }
  
  nuc_len <- nchar(ref)
  
  # Convert seqs
  if( convert.seq ) seqs <- divSeq(seqs, ref, fill = fill)
  
  # Nucleotide color
  nucleotide_levels <- c("A", "T", "G", "C", ".", "N")
  nucleotide_colors <- RColorBrewer::brewer.pal(6, "Set1")
  nucleotide_colors <- c(nucleotide_colors[c(1,3,6,2)], "#FFFFFF", "#DCDCDC")
  names(nucleotide_colors) <- nucleotide_levels
  
  # Sequence matrix
  N_pos <- which(unlist(strsplit(ref, "")) == "N")
  
  nuc_melt <- stringr::str_split(
      string = c(ref, paste(rep(" ", nuc_len), collapse = ""), seqs), 
      pattern = "", 
      simplify = TRUE
    ) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    dplyr::mutate(pos.y = -(seq_len(n()))) %>%
    tidyr::gather(key = "var", value = "value", -pos.y) %>%
    dplyr::mutate(
      pos.x = as.numeric(stringr::str_extract(var, "[0-9]+$")),
      color = nucleotide_colors[value],      
      color = ifelse(
        pos.x %in% N_pos, rep(nucleotide_colors["N"], n()), color),
      color = ifelse(value == " ", "#FFFFFF", color)
    ) %>%
    dplyr::select(pos.x, pos.y, value, color)
  
  # Format remaining cols of input
  sup_df <- df[,-match(nuc.col, names(df))] %>%
    dplyr::mutate_all(format, big.mark = ",", justify = "centre")
  
  sup_names <- names(sup_df)
  if( length(sup_names) > 0 ){
    
    sup_df <- dplyr::bind_rows(
        as.data.frame(
          t(matrix(
            c(names(sup_df), rep(" ", ncol(sup_df))), 
            ncol = 2, dimnames = list(names(sup_df))))),
        sup_df
      ) %>%
      dplyr::mutate_all(format, justify = "centre") %>%
      dplyr::mutate(pos.y = -(seq_len(n())))
    
    sup_melt <- tidyr::gather(sup_df, key = "var", value = "value", -pos.y) %>%
      dplyr::mutate(
        pos.x = nuc_len + (match(var, names(sup_df))) * padding - padding * 0.25,
        color = "#FFFFFF"
      ) %>%
      dplyr::select(pos.x, pos.y, value, color) %>%
      dplyr::bind_rows(
        data.frame(
          pos.x = max(.$pos.x) + padding,
          pos.y = -1,
          value = " ",
          color = "#FFFFFF"
        )
      )
  }else{
    
    sup_melt <- data.frame()
    
  }
  
  plot_melt <- dplyr::bind_rows(nuc_melt, sup_melt)
  
  plot_colors <- structure(
    unique(plot_melt$color), 
    names = unique(plot_melt$color)
  )
  
  p <- ggplot2::ggplot(plot_melt, ggplot2::aes(x = pos.x, y = pos.y)) +
    ggplot2::geom_tile(ggplot2::aes(fill = color)) +
    ggplot2::geom_text(
      ggplot2::aes(label = value), size = text.size, 
      family = font.family, fontface = font.face
    ) +
    ggplot2::scale_fill_manual(values = plot_colors) +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      legend.position = "none"
    )
  
  if( force.sq ){
    p <- p + 
      ggplot2::theme(aspect.ratio = with(plot_melt, max(abs(pos.y))/max(pos.x)))
  }
  
  p
  
}

#' Load Sample Info files given a run / project config file and root directory.
#'
#' \code{loadSampleInfo} returns a data.frame of the sampleInfo for the specific
#' run or project configuration file.
#'
#' @description Given a configuration file and root directory, the function
#' loads the sampleInfo file specified by the configuration file into the 
#' environment as a data.frame. This is accomplished with the data.table::fread
#' function to allow for variable input formats.
#'
#' @usage
#' loadSampleInfo(config, root.dir)
#'
#' @param config an imported yaml object from an iGUIDE run / project 
#' configuration.
#' 
#' @param root.dir path to iGUIDE install directory.
#' 
loadSampleInfo <- function(config, root.dir){
  
  # Check for required component config object
  stopifnot("Sample_Info" %in% names(config))
  
  # Sample_Info path
  si_path <- config[["Sample_Info"]]
  
  # Load file into environment through relative or absolute paths
  if( file.exists(file.path(root.dir, si_path)) ){
    
    return(data.table::fread(file.path(root.dir, si_path), data.table = FALSE))
    
  }else if( file.exists(si_path) ){
    
    return(data.table::fread(si_path, data.table = FALSE))
    
  }else{
    
    stop("\n  Cannot find Sample_Info: ", si_path, ".\n")
    
  }
  
}

#' Return an object specifying how specimens were experimentally treated.
#'
#' \code{getTreatmentInfo} returns a data.frame or list object of indicating
#' the experimental treatment for each specimen.
#'
#' @description Given a configuration file and root directory, the function will
#' return either a data.frame (default, \code{return.obj = "data.frame"}) or 
#' list (\code{return.obj = "list"}) indicating the associated experimental 
#' treatment stategies used for each specimen.
#'
#' @usage
#' getTreatmentInfo(config, root.dir)
#' getTreatmentInfo(config, root.dir, return.obj = "data.frame")
#'
#' @param config an imported yaml object from an iGUIDE run / project 
#' configuration.
#' 
#' @param root.dir path to iGUIDE install directory.
#' 
#' @param return.obj character either "data.frame" or "list" indicating the 
#' prefered object to return.
#'
getTreatmentInfo <- function(config, root.dir, return.obj = "data.frame"){
  
  # Check for required components within config
  stopifnot(all(
    c("Treatment", "Sample_Info", "Sample_Name_Column") %in% names(config)
  ))
  
  # Check for correct input for return object
  stopifnot( return.obj %in% c("data.frame", "list") )
  
  # Identify treatment input data
  treatments <- config[["Treatment"]]
  
  sample.info <- loadSampleInfo(config, root.dir)
  
  all_specimens <- unique(stringr::str_extract(
    string = sample.info[,config$Sample_Name_Column], pattern = "[\\w]+"
  ))
  
  # Determine input format: sampleInfo / all / list / incorrect
  if( any(grepl("sampleInfo:", treatments)) ){
    
    if( tolower(names(treatments)) != "all" ){
      stop("Config info should be specified for 'all' specimens if in sampleInfo.")
    }
    
    info_col <- match(
      stringr::str_extract(string = treatments, pattern = "[\\w]+$"), 
      names(sample.info)
    )
    
    if( length(info_col) != 1 ){
      stop("\n  Cannot parse treatment data. Check config yaml and sampleInfo.\n")
    }
    
    treatment_df <- data.frame(
      #run_set = sample_info$run_set,
      sampleName = sample.info[,sample_name_col], 
      treatment = sample.info[,info_col]
    )
    
  }else if( any(grepl("all", names(treatments))) ){
    
    stopifnot(length(treatments) == 1)
    
    treatment_df <- data.frame(
      sampleName = sample.info[,config$Sample_Name_Column], 
      treatment = unname(unlist(treatments))
    )
    
  }else if( all(names(treatments) %in% all_specimens) ){
    
    if( any(lengths(treatments) > 1) ){
      treatments <- lapply(treatments, paste, collapse = ";")
    }
    
    treatment_df <- data.frame(
      sampleName = names(treatments),
      treatment = unlist(treatments)
    )
    
  }else{
    
    stop(
      "\n  Treatment information not accurately parsed from config(s).\n", 
      "  Check config(s) formating.\n"
    )
    
  }
  
  # Format output objects
  treatment_df <- treatment_df %>%
    dplyr::mutate(
      specimen = stringr::str_extract(string = sampleName, pattern = "[\\w]+")
    ) %>%
    dplyr::distinct(specimen, treatment) %>%
    dplyr::mutate(treatment = ifelse(is.na(treatment), "Mock", treatment)) %>%
    dplyr::select(specimen, treatment)
  
  treatment <- strsplit(treatment_df$treatment, ";")
  names(treatment) <- treatment_df$specimen
  
  if( return.obj == "data.frame" ){
    return(treatment_df)
  }else{
    return(treatment)
  }
  
}

#' Return an object specifying how specimens which nucleases were used 
#' experimentally.
#'
#' \code{getNucleaseInfo} returns a data.frame or list object of indicating
#' the associated nuclease for each specimen.
#'
#' @description Given a configuration file and root directory, the function will
#' return either a data.frame (default, \code{return.obj = "data.frame"}) or 
#' list (\code{return.obj = "list"}) indicating the associated experimental 
#' treatment stategies used for each specimen.
#'
#' @usage
#' getNucleaseInfo(config, root.dir)
#' getNucleaseInfo(config, root.dir, return.obj = "data.frame")
#'
#' @param config an imported yaml object from an iGUIDE run / project 
#' configuration.
#' 
#' @param root.dir path to iGUIDE install directory.
#' 
#' @param return.obj character either "data.frame" or "list" indicating the 
#' prefered object to return.
#'
getNucleaseInfo <- function(config, root.dir, return.obj = "data.frame"){
  
  # Check for required components within config
  stopifnot(all(
    c("Nuclease", "Sample_Info", "Sample_Name_Column") %in% names(config)
  ))
  
  # Check for correct input for return object
  stopifnot( return.obj %in% c("data.frame", "list") )
  
  # Identify treatment input data
  nucleases <- config[["Nuclease"]]
  
  sample.info <- loadSampleInfo(config, root.dir)
  
  all_specimens <- unique(stringr::str_extract(
    string = sample.info[,config$Sample_Name_Column], pattern = "[\\w]+"
  ))
  
  # Determine input format: sampleInfo / all / list / incorrect
  if( any(grepl("sampleInfo:", nucleases)) ){
    
    if( tolower(names(nucleases)) != "all" ){
      stop("Config info should be specified for 'all' specimens if in sampleInfo.")
    }
    
    info_col <- match(
      stringr::str_extract(string = nucleases, pattern = "[\\w]+$"), 
      names(sample.info)
    )
    
    if( length(info_col) != 1 ){
      stop("\n  Cannot parse nucleases data. Check config yaml and sampleInfo.\n")
    }
    
    nuclease_df <- data.frame(
      sampleName = sample.info[,sample_name_col], 
      nuclease = sample.info[,info_col]
    )
    
  }else if( any(grepl("all", names(nucleases))) ){
    
    stopifnot(length(nucleases) == 1)
    
    nuclease_df <- data.frame(
      sampleName = sample.info[,config$Sample_Name_Column], 
      nuclease = unname(unlist(nucleases))
    )
    
  }else if( all(names(nucleases) %in% all_specimens) ){
    
    if( any(lengths(nucleases) > 1) ){
      nucleases <- lapply(nucleases, paste, collapse = ";")
    }
    
    nuclease_df <- data.frame(
      sampleName = names(nucleases),
      nuclease = unlist(nucleases)
    )
    
  }else{
    
    stop(
      "\n  Treatment information not accurately parsed from config(s).\n", 
      "  Check config(s) formating.\n"
    )
    
  }
  
  # Format output objects
  nuclease_df <- nuclease_df %>%
    dplyr::mutate(
      specimen = stringr::str_extract(string = sampleName, pattern = "[\\w]+")
    ) %>%
    dplyr::distinct(specimen, nuclease) %>%
    dplyr::mutate(nuclease = ifelse(is.na(nuclease), "Mock", nuclease)) %>%
    dplyr::select(specimen, nuclease)
  
  nuclease <- strsplit(nuclease_df$nuclease, ";")
  names(nuclease) <- nuclease_df$specimen
  
  if( return.obj == "data.frame" ){
    return(nuclease_df)
  }else{
    return(nuclease)
  }
  
}
