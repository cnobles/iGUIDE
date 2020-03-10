#' Format number in tables with big.marks conveinently
#'
#' @usage pNums(x, ...)
#'
#' @param x vector of numbers to format. Standard big.mark = ",".
#' @param ... arguments further passed to the `format` function.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

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
#' @author Christopher Nobles, Ph.D.
#' @export

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
#' @author Christopher Nobles, Ph.D.
#' @export

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
  seqs[GenomicRanges::match(seq_ranges, uniq_ranges)]

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
#' @author Christopher Nobles, Ph.D.
#' @export

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
        match(seq_names, GenomeInfoDb::seqnames(ref.genome))
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
#' @author Christopher Nobles, Ph.D.
#' @export

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
#' @author Christopher Nobles, Ph.D.
#' @export
#' @importFrom magrittr %>%

calcCoverage <- function(gr, resolution){

  #Set up coverage gr
  strandless <- gr
  GenomicRanges::strand(strandless) <- "*"
  gr_ranges <- range(strandless)
  grl_ranges <- GenomicRanges::split(
    gr_ranges, GenomeInfoDb::seqnames(gr_ranges), drop = TRUE
  )

  window_seqs <- BiocGenerics::unlist(BiocGenerics::lapply(
    grl_ranges,
    function(chr, res){
      seq(GenomicRanges::start(chr), GenomicRanges::end(chr), res)
    },
    res = resolution
  ))

  coverage_grl <- GenomicRanges::GRangesList(lapply(
    seq_along(gr_ranges),
    function(i, gr_ranges, window_seqs){

      seqname <- GenomicRanges::seqnames(gr_ranges[i])
      window <- window_seqs[[i]]

      GenomicRanges::GRanges(
        seqnames = rep(seqname, length(window)),
        ranges = IRanges::IRanges(
          start = window, width = rep(resolution, length(window))),
        strand = rep("*", length(window))
      )

    },
    gr_ranges = gr_ranges,
    window_seqs = window_seqs
  ))

  coverage_pos <- coverage_grl

  coverage_pos <- GenomicRanges::GRangesList(lapply(
    coverage_pos,
    function(x){
      GenomicRanges::strand(x) <- rep("+", length(x))
      x
    }
  ))

  coverage_neg <- coverage_grl

  coverage_neg <- GenomicRanges::GRangesList(lapply(
    coverage_pos,
    function(x){
      GenomicRanges::strand(x) <- rep("-", length(x))
      x
    }
  ))

  dplyr::bind_rows(lapply(
    seq_along(coverage_grl),
    function(i, gr){

      as.data.frame(coverage_grl[[i]], row.names = NULL) %>%
        dplyr::select(seqnames, start, end, width) %>%
        dplyr::mutate(
          readCountsPos = GenomicRanges::countOverlaps(coverage_pos[[i]], gr),
          readCountsNeg = GenomicRanges::countOverlaps(coverage_neg[[i]], gr)
        ) %>%
        dplyr::arrange(seqnames)

    },
    gr = gr
  ))

}

#' Combine a list of ShortRead objects
#'
#' @param split.seqs list of ShortRead objects
#'
#' @description Utility function that combines ShortRead objects.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

serialAppendS4 <- function(split.seqs){
  stopifnot(class(split.seqs) == "list")
  Reduce(f = ShortRead::append, x = split.seqs)
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
#' @author Christopher Nobles, Ph.D.
#' @export
#' @importFrom magrittr %>%

clusterKV <- function(key, val, return = "standard"){

  # Check inputs
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
#' @author Christopher Nobles, Ph.D.
#' @export

generateGenomicRegions <- function(ref, res, drop.alt.chr = TRUE){

  if( class(ref) == "BSgenome" ) ref <- GenomeInfoDb::seqinfo(ref)

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
#'   "C" = LETTERS[2:6]
#' )
#' vcollapse(df)
#' vcollapse(df, sep = "-")
#' vcollapse(df, sep = "-", fill = "z")
#'
#' @author Christopher Nobles, Ph.D.
#' @export

vcollapse <- function(d, sep = "", fill = "NA"){

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

#' Remove matching sequences from input vector and replace with character
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
#' @author Christopher Nobles, Ph.D.
#' @export

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
#' @author Christopher Nobles, Ph.D.
#' @export

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

#' Create Nuclease and Treatment Symbols
#' 
#' @usage combo_symbols(x)
#' 
#' @param x integer vector, commonly generated by `seq_len()` or `seq_along`.
#' 
#' @description Generates a vector of letters for the indices that are passed in
#' as "x".  Outputs an Alpha index ("A", "B", "C") that is appropriate for the 
#' length of indices queried. If great than 26 indicies are queried, the output
#' changes to double letter identifiers ("AA", "AB", "AC", ...). And so on.
#' 
#' @author Christopher Nobles, Ph.D.
#' @export
#' 
combo_symbols <- function(x){
  y <- LETTERS
  while (length(x) > length(y)) y <- sapply(y, function(z) paste0(z, LETTERS))
  y[x]
}


#' A Binary Ambiguous Nucleotide scoring Matrix (BAN Mat)
#'
#' @usage banmat()
#'
#' @description Meant for comparing ambiguous sequences against "A", "T", "G",
#' "C", and "N" containing sequences. Currently matches between ambiuous
#' nucleotides are considered mismatch. Constructed based on NUC4.4. Function
#' outputs a single square matrix.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

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

#' A Symmetric Scaled Ambiguous Nucleotide scoring Matrix (SSAN Mat)
#'
#' @usage ssanmat()
#'
#' @description Meant for comparing ambiguous sequences against other sequences
#' possibly containing ambiguous nucleotides. The scores are scaled between 0
#' and 1 for ambiguous matches based on overlaping proportions and the matrix is
#' symmetric. Constructed originally from NUC4.4. Function outputs a single
#' square matrix.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

ssanmat <- function(){
  matrix(c(
    1,0,0,0, 0   , 1   , 1   , 0   , 0   , 1   , 0   , 1   , 1   , 1   , 1,
    0,1,0,0, 0   , 1   , 0   , 1   , 1   , 0   , 1   , 0   , 1   , 1   , 1,
    0,0,1,0, 1   , 0   , 1   , 0   , 1   , 0   , 1   , 1   , 0   , 1   , 1,
    0,0,0,1, 1   , 0   , 0   , 1   , 0   , 1   , 1   , 1   , 1   , 0   , 1,
    0,0,1,1, 1   , 0   , 0.5 , 0.5 , 0.5 , 0.5 , 0.67, 0.67, 0.33, 0.33, 1,
    1,1,0,0, 0   , 1   , 0.5 , 0.5 , 0.5 , 0.5 , 0.33, 0.33, 0.67, 0.67, 1,
    1,0,1,0, 0.5 , 0.5 , 1   , 0   , 0.5 , 0.5 , 0.33, 0.67, 0.33, 0.67, 1,
    0,1,0,1, 0.5 , 0.5 , 0   , 1   , 0.5 , 0.5 , 0.67, 0.33, 0.67, 0.33, 1,
    0,1,1,0, 0.5 , 0.5 , 0.5 , 0.5 , 1   , 0   , 0.67, 0.33, 0.33, 0.67, 1,
    1,0,0,1, 0.5 , 0.5 , 0.5 , 0.5 , 0   , 1   , 0.33, 0.67, 0.67, 0.33, 1,
    0,1,1,1, 0.67, 0.33, 0.33, 0.67, 0.67, 0.33, 1   , 0.67, 0.67, 0.67, 1,
    1,0,1,1, 0.66, 0.33, 0.67, 0.33, 0.33, 0.67, 0.67, 1   , 0.67, 0.67, 1,
    1,1,0,1, 0.33, 0.67, 0.33, 0.67, 0.33, 0.67, 0.67, 0.67, 1   , 0.67, 1,
    1,1,1,0, 0.33, 0.67, 0.67, 0.33, 0.67, 0.33, 0.67, 0.67, 0.67, 1   , 1,
    1,1,1,1, 1   , 1   , 1   , 1   , 1   , 1   , 1   , 1   , 1   , 1   , 1),
    ncol = 15,
    nrow = 15,
    byrow = TRUE,
    dimnames = list(
      c("A", "T", "G", "C", "S", "W", "R", "Y",
        "K", "M", "B", "V", "H", "D", "N"),
      c("A", "T", "G", "C", "S", "W", "R", "Y",
        "K", "M", "B", "V", "H", "D", "N")))
}

#' An Asymmetric Scaled Ambiguous Nucleotide scoring Matrix (USAN Mat)
#'
#' @usage asanmat()
#'
#' @description Meant for comparing ambiguous sequences against other sequences
#' possibly containing ambiguous nucleotides. The scores are scaled between 0
#' and 1 for ambiguous matches based on overlaping proportions and the matrix is
#' asymmetric. This is meant for proper scoreing between different ambiguous
#' bases. For example, comparing a "S" (meaning a "C" or "G") to a "V" (meaning
#' a "C", "G", or "A") would be if either nucleotide was chosen (score = 1), but
#' comparing a "V" to an "S" would only be correct 2 out of three times
#' (score = 0.67).Constructed originally from NUC4.4. Function outputs a single
#' square matrix.
#'
#' @author Christopher Nobles, Ph.D.
#' @export

asanmat <- function(){
  matrix(c(
    # A T G C  S     W     R     Y     K     M     B     V     H     D     N
    1,0,0,0, 0   , 1   , 1   , 0   , 0   , 1   , 0   , 1   , 1   , 1   , 1, #A
    0,1,0,0, 0   , 1   , 0   , 1   , 1   , 0   , 1   , 0   , 1   , 1   , 1, #T
    0,0,1,0, 1   , 0   , 1   , 0   , 1   , 0   , 1   , 1   , 0   , 1   , 1, #G
    0,0,0,1, 1   , 0   , 0   , 1   , 0   , 1   , 1   , 1   , 1   , 0   , 1, #C
    0,0,1,1, 1   , 0   , 0.5 , 0.5 , 0.5 , 0.5 , 1   , 1   , 0.5 , 0.5 , 1, #S
    1,1,0,0, 0   , 1   , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 , 1   , 1   , 1, #W
    1,0,1,0, 0.5 , 0.5 , 1   , 0   , 0.5 , 0.5 , 0.5 , 1   , 0.5 , 1   , 1, #R
    0,1,0,1, 0.5 , 0.5 , 0   , 1   , 0.5 , 0.5 , 1   , 0.5 , 1   , 0.5 , 1, #Y
    0,1,1,0, 0.5 , 0.5 , 0.5 , 0.5 , 1   , 0   , 1   , 0.5 , 0.5 , 1   , 1, #K
    1,0,0,1, 0.5 , 0.5 , 0.5 , 0.5 , 0   , 1   , 0.5 , 1   , 1   , 0.5 , 1, #M
    0,1,1,1, 0.67, 0.33, 0.33, 0.67, 0.67, 0.33, 1   , 0.67, 0.67, 0.67, 1, #B
    1,0,1,1, 0.66, 0.33, 0.67, 0.33, 0.33, 0.67, 0.67, 1   , 0.67, 0.67, 1, #V
    1,1,0,1, 0.33, 0.67, 0.33, 0.67, 0.33, 0.67, 0.67, 0.67, 1   , 0.67, 1, #H
    1,1,1,0, 0.33, 0.67, 0.67, 0.33, 0.67, 0.33, 0.67, 0.67, 0.67, 1   , 1, #D
    1,1,1,1, 1   , 1   , 1   , 1   , 1   , 1   , 1   , 1   , 1   , 1   , 1),#N
    ncol = 15,
    nrow = 15,
    byrow = TRUE,
    dimnames = list(
      c("A", "T", "G", "C", "S", "W", "R", "Y",
        "K", "M", "B", "V", "H", "D", "N"),
      c("A", "T", "G", "C", "S", "W", "R", "Y",
        "K", "M", "B", "V", "H", "D", "N")))
}
