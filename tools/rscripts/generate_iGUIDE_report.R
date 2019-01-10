#' For anyone reviewing the code below, the following is a small style guide 
#' outlining the various formats for the code. 
#' 
#' Names with "_": objects, inlucding data.frames, GRanges, vectors, ...
#' Names in caMel format: functions or components of objects (i.e. columns 
#' within a data.frame).
#' Names with ".": arguments / options for functions

options(stringsAsFactors = FALSE, scipen = 99, width = 999)


# Set up and gather command line arguments ----
parser <- argparse::ArgumentParser(
  description = "Generate an iGUIDE report for input run(s)."
)

parser$add_argument(
  "input", nargs = "+", type = "character",
  help = "Post-processing output .rds file(s)."
)

parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", 
  help = "Output report file, extension not required."
)

parser$add_argument(
  "-c", "--config", nargs = "+", type = "character",
  help = "Run specific config file(s) in yaml format."
)

parser$add_argument(
  "-s", "--support", nargs = "+", type = "character",
  help = "Supplementary data input, csv or tsv format."
)

parser$add_argument(
  "-f", "--figures", action = "store_true",
  help = "Generate figures along with output report (pdf and png formats)."
)

parser$add_argument(
  "-d", "--data", action = "store_true",
  help = "Data to generate the report will be saved as an R image with output."
)

parser$add_argument(
  "-t", "--format", nargs = 1, type = "character", default = "html",
  help = "Output format for report. Either 'pdf' or 'html' (default)."
)

parser$add_argument(
  "-g", "--graphic", action = "store_true",
  help = "Includes an opening graphic on the report."
)

parser$add_argument(
  "--install_path", nargs = 1, type = "character", default = "IGUIDE_DIR",
  help = "iGUIDE install directory path, do not change for normal applications."
)


args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

if( length(args$input) != length(args$config) ){
  stop("Must supply one config file for each input.")
}

report_formats <- c("html" = "html_document", "pdf" = "pdf_document")

if( !args$format %in% names(report_formats) ){
  stop("Please input either 'html' or 'pdf' for format.\n",
       "Other formats not supported.")
}

output_format <- report_formats[args$format]

input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(seq_along(args), function(i){
    paste(args[[i]], collapse = ", ")
  })
)

input_table <- input_table[
  match(
    c("input :", "output :", "config :", "support :", 
      "figures :", "data :", "graphic :", "format :"),
    input_table$Variables),
]

cat("\niGUIDE Report Inputs:")

print(
  data.frame(input_table),
  right = FALSE, 
  row.names = FALSE
)


# Load dependancies ----
cat("\nLoading dependencies.")

add_packs <- c("magrittr", "knitr")

add_packs_loaded <- suppressMessages(
  sapply(add_packs, require, character.only = TRUE)
)

if( !all(add_packs_loaded) ){
  
  print(
    data.frame(
      "R-Packages" = names(add_packs_loaded), 
      "Loaded" = add_packs_loaded
    ), 
    right = FALSE,
    row.names = FALSE
  )
  
  stop("Check dependancies.")
  
}

code_dir <- dirname(sub(
  pattern = "--file=", 
  replacement = "", 
  x = grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
))

source(file.path(code_dir, "supporting_scripts/post_process_support.R"))

## Additional supporting functions ----
generate_genomic_regions <- function(ref, res, drop.alt.chr = TRUE){
  
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

genomic_density <- function(gr, res, cutoff = 2, adj = 1, drop.alt.chr = TRUE){
  
  stopifnot(class(gr) == "GRanges")
  
  if(
    is.na(sum(is.numeric(GenomeInfoDb::seqlengths(GenomicRanges::seqinfo(gr)))))
  ){
    stop("SeqInfo should be present for input GRanges.")
  }
  
  ref_regions <- generate_genomic_regions(
    GenomicRanges::seqinfo(gr), res, drop.alt.chr = drop.alt.chr
  )
  
  ref_regions$count <- GenomicRanges::countOverlaps(ref_regions, gr) + adj
  ref_regions$log.count <- log(ref_regions$count, base = 10)
  ref_regions$norm.log.count <- ref_regions$log.count/max(ref_regions$log.count)
  ref_regions <- ref_regions[ref_regions$count >= cutoff]
  ref_regions
  
}

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

plot_genomic_density <- function(grl, res, grp.col = NULL, cutoff = 2, 
                                 drop.alt.chr = TRUE, clean = FALSE){
  
  if( class(grl) == "GRanges" ) grl <- GenomicRanges::GRangesList(grl)
  
  if( !is.null(grp.col) ){
    
    stopifnot(all(unlist(lapply(
      grl, 
      function(x) grp.col %in% names(GenomicRanges::mcols(x))
    ))))
    
  }
  
  ref_len <- GenomeInfoDb::seqlengths(GenomicRanges::seqinfo(grl))
  
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
            y <- genomic_density(
              y, res = res, cutoff = cutoff, drop.alt.chr = drop.alt.chr
            )
            
            df <- as.data.frame(y, row.names = NULL) 
            
            cum_adj_pos <- ref_cum_len[
              match(df$seqnames, names(ref_cum_len)) - 1
            ]
            
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
  
  
  gen_den$score <- as.numeric(gen_den$type) + gen_den$norm.log.count
  
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
  
  y_breaks <- seq_along(grl)
  
  p <- ggplot(gen_den) 
  
  if( !clean ){
    
    p <- p +
      geom_hline(yintercept = y_breaks, color = "grey90") +
      geom_vline(xintercept = x_breaks, color = "grey90") +
      scale_x_continuous(
        breaks = x_lab_pos,
        labels = gsub("chr", "", names(x_lab_pos)))
    
  }else{
    
    p <- p + scale_x_continuous(labels = NULL)
    
  }
  
  p <- p + 
    geom_rect(
      aes(
        xmin = adj.start, xmax = adj.end, 
        ymin = as.integer(type), ymax = score, 
        fill = type)) +
    scale_y_continuous(
      limits = c(0, length(grl)+1), breaks = seq_along(grl), labels = NULL) +
    scale_fill_brewer(type = "qual", palette = "Set1", direction = -1) +
    labs(x = "Chromosome", fill = "Levels") + 
    coord_polar() +
    theme_bw() +
    theme(
      panel.background = element_rect(color = "white"),
      panel.border = element_rect(color = "white"),
      panel.grid = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks = element_blank())
  
  if( clean ){
    p + theme(legend.position = "none", axis.title = element_blank())
  }else{
    p
  }
  
}

div_seq <- function(seqs, ref, match.chr = "."){
  
  seqs <- as.character(seqs)
  ref <- as.character(ref)
  
  stopifnot( all(nchar(seqs) == nchar(ref)) )
  
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

seq_diverge_plot <- function(df, ref, nuc.col = NULL, padding = 4, 
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
  if( convert.seq ) seqs <- div_seq(seqs, ref)
  
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
    as.data.frame() %>%
    dplyr::mutate(pos.y = -(seq_len(n()))) %>%
    tidyr::gather(key = "var", value = "value", -pos.y) %>%
    #reshape2::melt(id.vars = "pos.y") %>%
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
    #melt(sup_df, id.vars = "pos.y") %>%
    dplyr::mutate(
      pos.x = nuc_len + (match(var, names(sup_df))) * padding,
      color = "#FFFFFF"
    ) %>%
    dplyr::select(pos.x, pos.y, value, color)
  
  plot_melt <- dplyr::bind_rows(nuc_melt, sup_melt)
  
  plot_colors <- structure(
    unique(plot_melt$color), 
    names = unique(plot_melt$color)
  )

  p <- ggplot(plot_melt, aes(x = pos.x, y = pos.y)) +
    geom_tile(aes(fill = color)) +
    geom_text(
      aes(label = value), size = text.size, 
      family = font.family, fontface = font.face) +
    scale_fill_manual(values = plot_colors) +
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "none")
  
  if( force.sq ){
    p <- p + theme(aspect.ratio = with(plot_melt, max(abs(pos.y))/max(pos.x)))
  }
  
  p
  
}


# Import metadata and consolidate into report objects ----
message("Importing experimental data and configurations.")

## Load config files
configs <- lapply(args$config, yaml::yaml.load_file)
names(configs) <- sapply(configs, "[[", "Run_Name")

## Load reference genome 
if( grepl(".fa", unique(sapply(configs, "[[", "Ref_Genome"))) ){
  
  if( !file.exists(config$Ref_Genome) ){
    stop("Specified reference genome file not found.")
  }
  
  ref_file_type <- ifelse(
    grepl(".fastq", unique(sapply(configs, "[[", "Ref_Genome"))), 
    "fastq", 
    "fasta"
  )
  
  ref_genome <- Biostrings::readDNAStringSet(
    filepath = unique(sapply(configs, "[[", "Ref_Genome")), 
    format = ref_file_type
  )
  
}else{
  
  ref_genome <- unique(sapply(configs, "[[", "Ref_Genome"))
  
  genome <- grep(
    pattern = ref_genome, 
    x = unique(BSgenome::installed.genomes()), 
    value = TRUE
  )
  
  if( length(genome) == 0 ){
    
    cat("\nInstalled genomes include:")
    print(unique(BSgenome::installed.genomes()))
    cat("\nSelected reference genome not in list.")
    stop("Error: Genome not available.")
    
  }else if( length(genome) > 1 ){
    
    cat("\nInstalled genomes include:")
    print(unique(BSgenome::installed.genomes()))
    cat(
      "\nPlease be more specific about reference genome.",
      "Multiple matches to input."
    )
    stop("Error: Multiple genomes requested.")
    
  }
  
  suppressMessages(library(genome, character.only = TRUE))
  
  ref_genome <- get(genome)
  
}

## Get versioning

root_dir <- Sys.getenv("IGUIDE_DIR")

soft_version <- as.character(read.delim(
  file = file.path(root_dir, ".version"), header = FALSE))

build_version <- list.files(file.path(root_dir, "etc")) %>%
  grep(pattern = "build", x = ., value = TRUE) %>%
  stringr::str_extract(pattern = "v[0-9]+\\.[0-9]+.[0-9]+")


## Load reference files
ref_genes <- loadRefFiles(
  configs[[1]]$refGenes, 
  type = "GRanges", 
  freeze = configs[[1]]$Ref_Genome
)

onco_genes <- loadRefFiles(
  configs[[1]]$oncoGeneList, 
  type = "gene.list", 
  freeze = configs[[1]]$Ref_Genome
)

special_genes <- loadRefFiles(
  configs[[1]]$specialGeneList, 
  type = "gene.list", 
  freeze = config[[1]]$Ref_Genome
)

umitag_option <- all(unlist(lapply(configs, "[[", "UMItags")))

upstream_dist <- unique(sapply(configs, function(x) x$upstreamDist))
downstream_dist <- unique(sapply(configs, function(x) x$downstreamDist))

if( length(upstream_dist) > 1 | length(downstream_dist) > 1 ){
  stop(
    "Inconsistant upstream or downstream distances between config files.\n",
    "Comparisons between groups with different run specific criteria\n", 
    "is not recommended.")
}

## Combine sampleInfo files

sample_info <- dplyr::bind_rows(lapply(
      sapply(configs, "[[", "Sample_Info"), 
    function(x){
      data.table::fread(x, data.table = FALSE)
    }), .id = "run_set")

sample_name_col <- unique(sapply(configs, "[[", "Sample_Name_Column"))

if( length(sample_name_col) != 1 ) stop("SampleInfo files not in same format.")

sample_info$specimen <- stringr::str_extract(
  string = sample_info[,sample_name_col], 
  pattern = "[\\w]+"
)

specimen_levels <- unique(sample_info$specimen)

sample_info$specimen <- factor(sample_info$specimen, levels = specimen_levels)

## Identify all gRNAs used from config files
gRNAs <- lapply(
  do.call(c, lapply(configs, "[[", "Guide_RNA_Sequences")), 
  toupper
)

gRNAs_grps <- stringr::str_extract(
  string = names(gRNAs), 
  pattern = "[\\w\\-\\_]+"
)

names(gRNAs) <- sub("[\\w\\-\\_]+.", "", names(gRNAs), perl = TRUE)
gRNAs <- split(gRNAs, gRNAs_grps)

pams <- lapply(
  do.call(c, lapply(configs, "[[", "PAM_Sequence")), 
  toupper
)

pams <- split(pams, stringr::str_extract(
  string = names(pams), 
  pattern = "[\\w\\-\\_]+"
))

gRNAs <- dplyr::bind_rows(
    mapply(
      function(seqs, pam){
        
        pam_mat <- matrix(unlist(
          lapply(
            pam, 
            function(pat) stringr::str_detect(seqs, paste0(pat,"$")) 
          )), 
          ncol = length(pam)
        )
        
        if( any(rowSums(pam_mat) > 1) ){ 
          stop("Multiple PAM sequences detected on a single guide RNA.") 
        }
        
        rownames(pam_mat) <- stringr::str_extract(
          string = names(seqs), 
          pattern = "[\\w\\-\\_\\.]+$"
        )
        
        colnames(pam_mat) <- unlist(pam)
        
        data.frame(
          row.names = rownames(pam_mat),
          "Guide" = rownames(pam_mat),
          "gRNA" = sapply(
            seq_along(seqs), 
            function(i){
              stringr::str_replace(
                string = seqs[[i]], 
                pattern = paste0(colnames(pam_mat)[pam_mat[i,]], "$"), 
                replacement = ""
              )
            }
          ),
          "PAM" = colnames(pam_mat)[
            sapply(
              seq_len(nrow(pam_mat)), 
              function(i) which(pam_mat[i,])
            )
          ]
        )
        
      }, 
      seqs = gRNAs, 
      pam = pams, 
      SIMPLIFY = FALSE
    ), 
    .id = "run.set"
  ) %>%
  dplyr::distinct(Guide, gRNA, PAM)

## Identify on-target edit sites from config files
on_targets <- unlist(lapply(configs, "[[", "On_Target_Sites"))
names(on_targets) <- stringr::str_replace(
  string = names(on_targets), pattern = stringr::fixed("."), replacement = ":"
)

names(on_targets) <- stringr::str_extract(
    string = names(on_targets), 
    pattern = "[\\w\\_\\-\\'\\.]+$"
  ) %>%
  stringr::str_extract(pattern = "[\\w\\_\\-\\.]+")

on_targets <- structure(
  unique(on_targets), 
  names = names(on_targets)[match(unique(on_targets), on_targets)]
)

## Treatment across runs
treatments <- lapply(configs, "[[", "Treatment")

if( any(grepl("sampleInfo:", treatments[[1]])) ){
  
  info_col <- match(
    stringr::str_extract(string = treatments[[1]], pattern = "[\\w]+$"), 
    names(sample_info)
  )
  
  if( length(info_col) != 1 ){
    stop("Cannot parse treatment data. Check config yaml and sampleInfo.")
  }
  
  treatment_df <- data.frame(
      sampleName = sample_info[,sample_name_col], 
      treatment = sample_info[,info_col]
    ) %>%
    dplyr::mutate(
      specimen = stringr::str_extract(string = sampleName, pattern = "[\\w]+"),
      specimen = factor(specimen, levels = specimen_levels)
    ) %>%
    dplyr::distinct(specimen, treatment) %>%
    dplyr::arrange(specimen) %>%
    dplyr::mutate(treatment = ifelse(is.na(treatment), "Mock", treatment)) %>%
    dplyr::select(specimen, treatment)
  
  treatment <- strsplit(treatment_df$treatment, ";")
  names(treatment) <- treatment_df$specimen
  
}else if( any(grepl("all", names(treatments[[1]]))) ){
  
  treatment_df <- data.frame(
      sampleName = sample_info[,sample_name_col], 
      treatment = unique(unlist(treatments))
    ) %>%
    dplyr::mutate(
      specimen = stringr::str_extract(string = sampleName, pattern = "[\\w]+"),
      specimen = factor(specimen, levels = specimen_levels)
    ) %>%
    dplyr::distinct(specimen, treatment) %>%
    dplyr::arrange(specimen) %>%
    dplyr::mutate(treatment = ifelse(is.na(treatment), "Mock", treatment)) %>%
    dplyr::select(specimen, treatment)
  
  treatment <- strsplit(treatment_df$treatment, ";")
  names(treatment) <- treatment_df$specimen
  
}else{
  
  stop(
    "Treatment information not accurately parsed from config(s).\n", 
    "Check config(s) formating.")
  
}

## Load in supporting information ----
if( length(args$support) > 0 ){
  
  supp_data <- data.table::fread(args$support, data.table = FALSE)
  specimen_levels <- supp_data$specimen[supp_data$specimen %in% specimen_levels]
  
  supp_data <- dplyr::filter(supp_data, specimen %in% specimen_levels) %>%
    dplyr::mutate(specimen = factor(specimen, levels = specimen_levels))
  
  treatment_df <- dplyr::filter(treatment_df, specimen %in% specimen_levels) %>%
    dplyr::mutate(
      specimen = factor(as.character(specimen), levels = specimen_levels)
    ) %>%
    dplyr::arrange(specimen)
  
  treatment <- treatment[names(treatment) %in% specimen_levels]
  sample_info <- dplyr::filter(sample_info, specimen %in% specimen_levels) %>%
    dplyr::mutate(
      specimen = factor(as.character(specimen), levels = specimen_levels)
    ) %>%
    dplyr::arrange(specimen)
  
}else{
  
  supp_data <- data.frame()
  
}


## Consolidate supplementary data ----
if( is.null(args$support) ){
  spec_overview <- treatment_df
}else{
  spec_overview <- supp_data
}

cond_overview <- spec_overview %>%
  dplyr::mutate(
    condition = vcollapse(
      d = dplyr::select(spec_overview, -specimen), 
      sep = " - ", 
      fill = "NA"
    ),
    condition = factor(condition, levels = c(unique(condition), "Mock"))
  ) %>%
  dplyr::select(specimen, condition)


## Read in experimental data and contatenate different sets ----
input_data <- lapply(args$input, readRDS)
names(input_data) <- names(configs)
data_type_names <- names(input_data[[1]])

input_data <- lapply(
  seq_along(input_data[[1]]), 
  function(i) dplyr::bind_rows(lapply(input_data, "[[", i), .id = "run.set")
)

names(input_data) <- data_type_names

input_data <- lapply(
  input_data, 
  function(x) x[x$specimen %in% spec_overview$specimen,] 
)


cat("\nStarting analysis.")

# Opening graphic ----
graphic_order <- c("algnmts", "pile_up_algns", "paired_algns", "matched_algns")
graphic_data <- input_data[graphic_order]

graphic_grl <- GenomicRanges::GRangesList(lapply(
  graphic_data, 
  GenomicRanges::makeGRangesFromDataFrame, 
  seqinfo = GenomicRanges::seqinfo(ref_genome)
))


# Specimen summary ----
# Summarize components and append to specimen table
tbl_algn_counts <- input_data$algnmts %>% 
  dplyr::mutate(specimen = factor(specimen, levels = specimen_levels)) %>%
  dplyr::group_by(specimen)


if( umitag_option ){
  
  tbl_algn_counts <- dplyr::summarise(
    tbl_algn_counts, 
    Reads = sum(count), UMItags = sum(umitag), Alignments = sum(contrib)
  )
  
}else{
  
  tbl_algn_counts <- dplyr::summarise(
    tbl_algn_counts, Reads = sum(count), Alignments = sum(contrib)
  )
  
}


spec_overview_join <- dplyr::left_join(
  spec_overview, tbl_algn_counts, by = "specimen"
) 


# On-target summary ----
# Algnmts
tbl_ot_algn <- input_data$algnmts %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(
    ot_algns = pNums(sum(contrib * as.integer(edit.site %in% on_targets))),
    ot_algns_pct = 100 * sum(contrib * as.integer(edit.site %in% on_targets)) /
      sum(contrib)
  ) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# Probable edited sites
tbl_ot_prob <- input_data$probable_algns %>% 
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(
    ot_prob = pNums(sum(contrib * as.integer(edit.site %in% on_targets))),
    ot_prob_pct = 100 * sum(contrib * as.integer(edit.site %in% on_targets)) /
      sum(contrib)
  ) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# Pile ups of read alignments
tbl_ot_pile <- input_data$pile_up_algns %>% 
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(
    ot_pile = pNums(sum(contrib * as.integer(edit.site %in% on_targets))),
    ot_pile_pct = 100 * sum(contrib * as.integer(edit.site %in% on_targets)) /
      sum(contrib)
  ) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# Paired or flanking algnments
tbl_ot_pair <- input_data$paired_regions %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen))),
    on.off.target = factor(
      on.off.target, levels = c("On-target", "Off-target")
    )
  ) %>%
  dplyr::group_by(specimen, on.off.target) %>%
  dplyr::summarise(cnt = sum(algns)) %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(
    ot_pair = pNums(sum(ifelse(on.off.target == "On-target", cnt, 0))),
    ot_pair_pct = 100 * sum(ifelse(on.off.target == "On-target", cnt, 0)) /
      sum(cnt)
  ) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# Guide RNA matched within 6 mismatches
tbl_ot_match <- input_data$matched_summary %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::group_by(specimen, on.off.target) %>%
  dplyr::summarise(cnt = sum(algns)) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(specimen) %>%
  dplyr::summarise(
    ot_match = pNums(sum(ifelse(on.off.target == "On-target", cnt, 0))),
    ot_match_pct = 100 * sum(ifelse(on.off.target == "On-target", cnt, 0)) /
      sum(cnt)
  ) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# Summary table
ot_tbl_summary <- dplyr::left_join(
    treatment_df, cond_overview, by = "specimen"
  ) %>%
  dplyr::mutate(
    condition = factor(
      ifelse(is.na(condition), "Mock", paste(condition)), 
      levels = levels(condition)
    )
  )

ot_tbl_summary <- Reduce(
    function(x,y) dplyr::left_join(x, y, by = "specimen"),
    list(
      tbl_ot_algn[,c(1,3)], tbl_ot_pile[,c(1,3)],
      tbl_ot_pair[,c(1,3)], tbl_ot_match[,c(1,3)]
    ),
    init = ot_tbl_summary
  ) %>%
  dplyr::mutate(specimen = factor(specimen, levels = specimen_levels)) %>%
  dplyr::arrange(specimen)


names(ot_tbl_summary) <- c(
  "Specimen", "Treatment", "Condition", "All\nAlign.", "Align.\nPileups", 
  "Flanking\nPairs", "gRNA\nMatched")


# On-target incorporation distribution ----
on_tar_dists <- input_data$matched_algns %>%
  dplyr::filter(on.off.target == "On-target") %>%
  dplyr::mutate(
    gRNA = stringr::str_extract(string = guideRNA.match, pattern = "[\\w]+"),
    pos = as.numeric(
      stringr::str_extract(string = edit.site, pattern = "[0-9]+$")
    ),
    edit.site.dist = ifelse(strand == "+", start - pos, end - pos),
    specimen = factor(specimen, levels = levels(cond_overview$specimen))
  ) %>%
  dplyr::left_join(cond_overview, by = "specimen") %>%
  dplyr::select(
    run.set, specimen, gRNA, condition, 
    edit.site, edit.site.dist, strand, contrib)

on_tar_dens <- lapply(
  split(on_tar_dists, on_tar_dists$condition), 
  function(x){
    
    if( nrow(x) >= 10 ){
      return(
        density(abs(x$edit.site.dist), from = 0, to = upstream_dist, bw = 1)
      )
    }else{
      return(NA)
    }
    
  }
)

on_tar_dists <- dplyr::group_by(
    on_tar_dists, 
    condition, gRNA, edit.site.dist, strand
  ) %>%
  dplyr::summarise(cnt = sum(contrib)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    strand.cnt = ifelse(
      strand == "+", log(cnt, base = 10), -log(cnt, base = 10))
  )

if( length(unique(cond_overview$condition)) == 1 ){
  on_tar_dists$condition <- " "
}

if( is.null(args$support) ){
  sites_included <- on_tar_dists %>% dplyr::group_by(gRNA)
}else{
  sites_included <- on_tar_dists %>% dplyr::group_by(condition, gRNA)
}

sites_included <- dplyr::summarise(
    sites_included,
    prop = 100 * sum(cnt[
      abs(edit.site.dist) <= upstream_dist & 
        abs(edit.site.dist) >= -downstream_dist
    ]) / sum(cnt),
    x_pos = upstream_dist,
    y_pos = 0.8 * min(strand.cnt[
      abs(edit.site.dist) <= upstream_dist & 
        abs(edit.site.dist) >= -downstream_dist
    ])
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(prop = paste0(pNums(prop, digits = 4), "%"))


# Off-target summary ----
# All alignments
tbl_ft_algn <- input_data$algnmts %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::filter(!edit.site %in% on_targets) %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(ft_algns = dplyr::n_distinct(clus.ori)) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# Probable edit sites
tbl_ft_prob <- input_data$probable_algns %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::filter(on.off.target == "Off-target") %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(ft_prob = dplyr::n_distinct(clus.ori)) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# Pile ups
tbl_ft_pile <- input_data$pile_up_algns %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::filter(on.off.target == "Off-target") %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(ft_pile = dplyr::n_distinct(clus.ori)) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# Paired or flanked loci
tbl_ft_pair <- input_data$paired_regions %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::filter(on.off.target == "Off-target") %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(ft_pair = n()) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# gRNA sequence matched
tbl_ft_match <- input_data$matched_summary %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = sort(unique(sample_info$specimen)))
  ) %>%
  dplyr::filter(on.off.target == "Off-target") %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise(ft_match = n()) %>%
  dplyr::ungroup() %>% 
  as.data.frame()

# Summary table
ft_tbl_summary <- dplyr::left_join(
    treatment_df, cond_overview, by = "specimen"
  ) %>%
  dplyr::mutate(
    condition = factor(
      ifelse(is.na(condition), "Mock", paste(condition)), 
      levels = levels(condition)
    )
  )

ft_tbl_summary <- Reduce(
    function(x,y) dplyr::left_join(x, y, by = "specimen"),
    list(tbl_ft_algn, tbl_ft_pile, tbl_ft_pair, tbl_ft_match),
    init = ft_tbl_summary
  ) %>%
  dplyr::mutate(specimen = factor(specimen, levels = specimen_levels)) %>%
  dplyr::arrange(specimen)

names(ft_tbl_summary) <- c(
  "Specimen", "Treatment", "Condition", "All\nAlign.", "Align.\nPileups", 
  "Flanking\nPairs", "gRNA\nMatched")


# Onco-gene enrichment analysis ----
rand_sites <- selectRandomSites(
  num = nrow(input_data$paired_regions) + nrow(input_data$matched_summary), 
  ref.genome = ref_genome, 
  drop.extra.seqs = TRUE, 
  rnd.seed = 714)

rand_sites$gene_id <- suppressMessages(assignGeneID(
  seqnames = seqnames(rand_sites), 
  positions = start(rand_sites), 
  reference = ref_genome, 
  ref.genes = ref_genes, 
  onco.genes = onco_genes, 
  special.genes = special_genes
))

rand_df <- data.frame(
  condition = "Random", 
  "total" = length(rand_sites), 
  "onco" = sum(stringr::str_detect(rand_sites$gene_id, "~")), 
  "special" = sum(stringr::str_detect(rand_sites$gene_id, "!"))
)

paired_list <- split(
  x = input_data$paired_regions, 
  f = cond_overview$condition[
    match(input_data$paired_regions$specimen, cond_overview$specimen)
  ]
)

paired_df <- dplyr::bind_rows(
  lapply(
    paired_list, 
    function(df){
    
    data.frame(
      "total" = nrow(df), 
      "onco" = sum(stringr::str_detect(df$gene_id, "~")), 
      "special" = sum(stringr::str_detect(df$gene_id, "!"))
    )
    
    }
  ), 
  .id = "condition"
)

matched_list <- split(
  x = input_data$matched_summary, 
  f = cond_overview$condition[
    match(input_data$matched_summary$specimen, cond_overview$specimen)
])

matched_df <- dplyr::bind_rows(
  lapply(
    matched_list, 
    function(df){
      
      data.frame(
        "total" = nrow(df), 
        "onco" = sum(stringr::str_detect(df$gene_id, "~")), 
        "special" = sum(stringr::str_detect(df$gene_id, "!"))
      )
      
    }
  ), 
  .id = "condition"
)

enrich_df <- dplyr::bind_rows(
    list(
      "Reference" = rand_df, 
      "Flanking Pairs" = paired_df, 
      "gRNA Matched" = matched_df), 
    .id = "origin"
  ) %>%
  dplyr::filter(total > 0)

enrich_df$onco.p.value <- p.adjust(
  sapply(
    seq_len(nrow(enrich_df)), 
    function(i){
      
      ref <- enrich_df[1, c("total", "onco")]
      query <- enrich_df[i, c("total", "onco")]
      ref$diff <- abs(diff(as.numeric(ref)))
      query$diff <- abs(diff(as.numeric(query)))
      
      fisher.test(as.matrix(rbind(
        ref[,c("diff", "onco")], query[,c("diff", "onco")]
      )))$p.value
      
    }
  ), 
  method = "BH"
)

enrich_df$special.p.value <- p.adjust(
  sapply(
    seq_len(nrow(enrich_df)), 
    function(i){
      
      ref <- enrich_df[1, c("total", "special")]
      query <- enrich_df[i, c("total", "special")]
      ref$diff <- abs(diff(as.numeric(ref)))
      query$diff <- abs(diff(as.numeric(query)))
      
      fisher.test(as.matrix(rbind(
        ref[,c("diff", "special")], query[,c("diff", "special")]
      )))$p.value
      
    }
  ), 
  method = "BH"
)

names(enrich_df) <- c(
  "Origin", "Condition", "Total Gene Count", "Onco Related Count", 
  "Special Gene Count", "Onco Enrich. p-value", "Special Enrich. p-value")


# Genomic Distribution of edited sites ----
genomic_grl <- GenomicRanges::GRangesList(lapply(
  graphic_data, 
  function(x){
    
    y <- makeGRangesFromDataFrame(x, seqinfo = seqinfo(ref_genome))
    mcols(y) <- cond_overview[
      match(x$specimen, cond_overview$specimen), "condition", drop = FALSE]
    y
    
  }
))

num_conds <- max(length(unique(cond_overview$condition)), 1)

names(genomic_grl) <- c(
  "All Align.", "Pileup Align.", "Flanking Pairs", "gRNA Matched")


# Off-target sequence analysis ----
ft_MESL <- input_data$matched_algns %>%
  dplyr::mutate(
    edit.site.dist = abs(ifelse(
      strand == "+", 
      start - as.numeric(stringr::str_extract(edit.site, "[0-9]+$")), 
      as.numeric(stringr::str_extract(edit.site, "[0-9]+$")) - end
    )),
    specimen = factor(specimen, levels = levels(cond_overview$specimen))
  ) %>%
  dplyr::left_join(cond_overview, by = "specimen") %>%
  dplyr::mutate(order = seq_len(n())) 

if( nrow(ft_MESL) > 0 ){
  
  ft_MESL <- ft_MESL %>%
    dplyr::group_by(order) %>%
    dplyr::mutate(
      ESL = predictESProb(
        x = edit.site.dist, density = on_tar_dens[[condition]]
      ),
      gene_id = input_data$matched_summary$gene_id[
        match(edit.site, input_data$matched_summary$edit.site)
      ]
    ) %>%
    dplyr::group_by(condition, edit.site, gene_id) %>%
    dplyr::summarise(MESL = 100 * max(c(0,ESL), na.rm = TRUE)) %>%
    dplyr::ungroup()
  
}else{
  
  ft_MESL <- ft_MESL %>%
    dplyr::mutate(
      MESL = vector(mode = "numeric"), gene_id = vector(mode = "character")) %>%
    dplyr::select(condition, edit.site, gene_id, MESL)
  
}

ft_seqs <- input_data$matched_summary %>%
  dplyr::select(
    specimen, aligned.sequence, guideRNA.match, edit.site,
    guideRNA.mismatch, on.off.target, algns, gene_id
  ) %>%
  dplyr::mutate(
    specimen = factor(specimen, levels = levels(cond_overview$specimen))
  ) %>%
  dplyr::left_join(cond_overview, by = "specimen") %>%
  dplyr::left_join(ft_MESL, by = c("condition", "edit.site", "gene_id"))

if( is.null(args$support) ){
  
  ft_seqs <- dplyr::group_by(
      ft_seqs, 
      guideRNA.match, edit.site, aligned.sequence, 
      guideRNA.mismatch, on.off.target, gene_id
    ) %>%
    dplyr::summarise(algns = sum(algns), MESL = max(MESL, na.rm = TRUE))
  
}else{
  
  ft_seqs <- dplyr::group_by(
      ft_seqs,
      condition, guideRNA.match, edit.site, aligned.sequence, 
      guideRNA.mismatch, on.off.target, gene_id
    ) %>%
    dplyr::summarise(algns = sum(algns), MESL = max(MESL, na.rm = TRUE))
  
}

ft_seqs <- dplyr::arrange(
    ft_seqs, desc(algns), desc(MESL), guideRNA.mismatch
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    on.off.target = stringr::str_extract(on.off.target, "[\\w]+")
  ) %>%
  dplyr::rename(
    "target" = on.off.target, 
    "mismatch" = guideRNA.mismatch, 
    "gRNA" = guideRNA.match,
    "aligns" = algns
  )

if( is.null(args$support) ){
  
  ft_seqs_list <- split(ft_seqs, ft_seqs$gRNA)
  
}else{
  
  ft_seqs_conds <- dplyr::arrange(ft_seqs, condition, gRNA) %$% 
    unique(paste0(condition, " - ", gRNA))
  
  ft_seqs_list <- split(
    x = ft_seqs, 
    f = factor(
      paste0(ft_seqs$condition, " - ", ft_seqs$gRNA), levels = ft_seqs_conds
    )
  )
  
}

cat("\nAnalysis complete. Starting report generation.")

# Data passed to Rmd for report generation ----
set_names <- ifelse(
  length(configs) ==  1, 
  names(configs), 
  paste0(
    paste(names(configs)[seq_len(length(configs)-1)], collapse = ", "),
    ", and ", 
    names(configs)[length(configs)]
  )
)

# Normalize file output path
write(c(), file = args$output)
args$output <- normalizePath(args$output)
unlink(args$output)

output_path <- unlist(strsplit(args$output, "/"))
output_dir <- paste(output_path[seq_len(length(output_path)-1)], collapse = "/")
output_file <- output_path[length(output_path)]

if( args$format == "html" & !stringr::str_detect(output_file, ".html$") ){
  output_file <- paste0(output_file, ".html")
}

if( args$format == "pdf" & !stringr::str_detect(output_file, ".pdf$") ){
  output_file <- paste0(output_file, ".pdf")
}

if( args$data ){

  if( args$format == "html" ){
    
    save.image(file = file.path(
      output_dir, stringr::str_replace(output_file, ".html$", ".RData")
    )) 
    
  }else if( args$format == "pdf" ){
    
    save.image(file = file.path(
      output_dir, stringr::str_replace(output_file, ".pdf$", ".RData")
    )) 
    
  }
  
}

if( args$format == "html" ){
  
  template_path <- normalizePath(
    file.path(code_dir, "iGUIDE_report_template.Rmd")
  )
  
  css_path <- normalizePath(file.path(code_dir, "iguide.css"))
  
  rmarkdown::render(
    input = template_path,
    output_format = output_format, 
    output_file = output_file,
    output_dir = output_dir,
    output_options = list("css" = css_path)
  )
  
}else{
  
  template_path <- normalizePath(
    file.path(code_dir, "iGUIDE_report_template.Rmd")
  )
  
  rmarkdown::render(
    input = template_path,
    output_format = output_format, 
    output_file = output_file,
    output_dir = output_dir
  )
  
}

q()
