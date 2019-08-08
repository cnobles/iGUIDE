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
#' @author Christopher Nobles, Ph.D.
#' @export

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
#' @author Christopher Nobles, Ph.D.
#' @export

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
#' @author Christopher Nobles, Ph.D.
#' @export

makeSquare <- function(p, dims, fudge=1){

  dims <- heatmap_dims(p)
  p + ggplot2::theme(aspect.ratio = (dims$nrows / dims$ncols) * fudge)

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
#' @author Christopher Nobles, Ph.D.
#' @export

plotGenomicDensity <- function(grl, res = 1E7, grp.col = NULL, cutoff = 2,
                               abund.method = "norm.log.count",
                               drop.alt.chr = TRUE, clean = FALSE){

  if( class(grl) == "GRanges" ) grl <- GenomicRanges::GRangesList(grl)

  if( !is.null(grp.col) ){

    stopifnot(all(BiocGenerics::unlist(BiocGenerics::lapply(
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

  gen_densities <- BiocGenerics::lapply(
    grl,
    function(x){

      if( is.null(grp.col) ){
        x <- list(x)
      }else{
        x <- GenomicRanges::split(
          x, GenomicRanges::mcols(x)[,grp.col, drop = FALSE]
        )
      }

      dplyr::bind_rows(BiocGenerics::lapply(
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
        .id = "cond"
      )

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
      GenomicRanges::split(
        gen_den[,match(abund.method, names(gen_den))], 
      gen_den$type
    ), 
    max
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
#' @author Christopher Nobles, Ph.D.
#' @export
#' @importFrom magrittr %>%

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
