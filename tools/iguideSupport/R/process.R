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
#' @author Christopher Nobles, Ph.D.
#' @export

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
  grl <- GenomicRanges::split(gr, gr$grouping)
  gr <- unname(BiocGenerics::unlist(GenomicRanges::GRangesList(lapply(
    grl,
    function(g){

      g_pos <- g[GenomicRanges::strand(g) == "+"]
      g_pos <- groupPileups(g_pos, "+", maxgap = maxgap)
      g_neg <- g[GenomicRanges::strand(g) == "-"]
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
#' @author Christopher Nobles, Ph.D.
#' @export

groupPileups <- function(gr, strand, maxgap){

  # Implement axial cluster structure rather than all vs all.
  red.gr <- GenomicRanges::reduce(gr, min.gapwidth = maxgap, with.revmap = TRUE)

  axil_nodes <- BiocGenerics::unlist(S4Vectors::Rle(
    values = BiocGenerics::unlist(red.gr$revmap)[
      S4Vectors::start(red.gr$revmap@partitioning)
    ],
    lengths = S4Vectors::width(red.gr$revmap@partitioning)
  ))

  nodes <- BiocGenerics::unlist(red.gr$revmap)
  edgelist <- unique(matrix( c(axil_nodes, nodes), ncol = 2 ))

  pile_ups <- igraph::clusters(igraph::graph.edgelist(
    el = edgelist, directed = FALSE
  ))

  gr$clusID <- igraph::membership(pile_ups)

  if( strand == "+" ){

    pile_starts <- GenomicRanges::sort(GenomicRanges::split(
      x = S4Vectors::Rle(GenomicRanges::start(gr)),
      f = gr$clusID
    ))

    pile_starts <- S4Vectors::runValue(pile_starts)[
      S4Vectors::runLength(pile_starts) ==
        max(S4Vectors::runLength(pile_starts))
    ]

    pile_starts <- BiocGenerics::unlist(pile_starts)[
      S4Vectors::start(pile_starts@partitioning)
    ]

  }else{

    pile_starts <- GenomicRanges::sort(GenomicRanges::split(
      x = S4Vectors::Rle(GenomicRanges::end(gr)),
      f = gr$clusID
    ))

    pile_starts <- S4Vectors::runValue(pile_starts)[
      S4Vectors::runLength(pile_starts) ==
        max(S4Vectors::runLength(pile_starts))
    ]

    pile_starts <- BiocGenerics::unlist(pile_starts)[
      S4Vectors::end(pile_starts@partitioning)
    ]

  }

  gr <- gr[order(gr$clusID)]

  clus.ori <- paste0(
    GenomicRanges::seqnames(gr), ":", GenomicRanges::strand(gr), ":",
    BiocGenerics::unlist(
      S4Vectors::Rle(values = pile_starts, lengths = pile_ups$csize)
    )
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
#' @author Christopher Nobles, Ph.D.
#' @export
#' @importFrom magrittr %>%

identifyPairedAlgnmts <- function(gr, grouping = NULL, maxgap, maxovlp = 10L){

  if( !is.null(grouping) ){
    group_vec <- GenomicRanges::mcols(gr)[
      ,match(grouping, names(GenomicRanges::mcols(gr)))
    ]
  }else{
    group_vec <- rep(1, length(gr))
  }

  gr$ori.order <- seq_along(gr)
  grl <- GenomicRanges::split(x = gr, f = group_vec)

  pairs <- unname(BiocGenerics::unlist(lapply(
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
          dplyr::filter(dplyr::n() > 1) %>%
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
#' @author Christopher Nobles, Ph.D.
#' @export
#' @importFrom magrittr %>%

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

  grl <- GenomicRanges::split(x = gr, f = group_vec)
  gr_mod <- BiocGenerics::unlist(GenomicRanges::GRangesList(lapply(
    grl,
    function(gs){

      pp <- GenomicRanges::as.data.frame(gs) %>%
        dplyr::distinct(pile.id, pair.id)

      pp_gr <- GenomicRanges::GRanges(
        seqnames = "mock",
        ranges = IRanges::IRanges(
          start = as.integer(factor(pp$pair.id)), width = 1
        ),
        strand = "+"
      )

      pp_grl <- GenomicRanges::split(x = pp_gr, f = pp$pile.id)
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
#' @author Christopher Nobles, Ph.D.
#' @export

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
#' @author Christopher Nobles, Ph.D.
#' @export
#' @importFrom magrittr %>%

alnTargetSeqs <- function(seqs, target.seqs, tolerance, fixed = 'subject'){

  if( is.null(names(seqs)) ) names(seqs) <- seq_along(seqs)
  if( class(seqs) != "DNAStringSet") seqs <- Biostrings::DNAStringSet(seqs)
  if( is.null(names(target.seqs)) ) names(target.seqs) <- target.seqs

  nt_widths <- data.frame(
    names = names(seqs),
    nt_width = Biostrings::width(seqs),
    stringsAsFactors = FALSE
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

  df <- dplyr::bind_rows(
    lapply(
      seq_along(alns),
      function(i){

        dplyr::bind_rows(
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
#' @author Christopher Nobles, Ph.D.
#' @export
#' @importFrom magrittr %>%

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
#' @author Christopher Nobles, Ph.D.
#' @export
#' @importFrom magrittr %>%

filterInappropriateComparisons <- function(guideRNA.match, specimen, treatment){

  stopifnot( length(guideRNA.match) == length(specimen) )

  df <- data.frame(
    order = seq_along(specimen),
    specimen = specimen,
    guideRNA = guideRNA.match,
    stringsAsFactors = FALSE
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
#' @author Christopher Nobles, Ph.D.
#' @export

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
#' @author Christopher Nobles, Ph.D.
#' @export

predictESProb <- function(z, density, range = NULL){

  if( is.null(range) & class(density) == "density" ) range <- range(density$x)

  if( class(density) == "density"){
    
    idx <- match(z, round(density$x))
    
    return(vapply(idx, function(i){
        if(is.na(i)){
          return(NA)
        }else{
          return(with(density, 1 - sum(y[seq_len(i)] * rep(diff(x)[1], i))))
        }
      },
      numeric(1)
    ))
    
  }else{
    
    return(rep(NA, length(z)))
    
  }

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
#' @author Christopher Nobles, Ph.D.
#' @export

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
#' @author Christopher Nobles, Ph.D.
#' @export
#' @importFrom magrittr %>%

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
#' @author Christopher Nobles, Ph.D.
#' @export
#' @importFrom magrittr %>%

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
