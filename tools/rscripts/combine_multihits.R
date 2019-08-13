#' For those reviewing the code below, the following is a small style guide 
#' outlining the various formats for the code. 
#' 
#' Names with "_": objects, inlucding data.frames, GRanges, vectors, ...
#' Names in caMel format: functions or components of objects (i.e. columns 
#' within a data.frame).
#' Names with ".": arguments / options for functions

#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, scipen = 99, width = 120)

code_dir <- dirname(sub(
  pattern = "--file=", 
  replacement = "", 
  x = grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
))

# Set up and gather command line arguments
parser <- argparse::ArgumentParser(
  description = "Script for combining multihit objects together.",
  usage = paste(
    "Rscript combine_multihits.R -d <directory> -p <pattern>",
    "[-h/--help, -v/--version] [optional args]"
  )
)

parser$add_argument(
  "-d", "--dir", nargs = 1, type = "character", 
  help = paste(
    "Directory where to look for multihit files. Combine with 'pattern'",
    "to select specific files."
  )
)

parser$add_argument(
  "-p", "--pattern", nargs = 1, type = "character", default = ".",
  help = paste(
    "Pattern to identify files within the directory specified to combine.",
    "Regex patterns supported through R. Default: '.'"
  )
)

parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", required = TRUE,
  help = "Output file name. Output format only supports R-based rds format.
)

parser$add_argument(
  "--stat", nargs = 1, type = "character", default = FALSE, 
  help = "Stat output name. Stats output in long format."
)


args <- parser$parse_args(commandArgs(trailingOnly = TRUE))


# Check output file name
if( stringr::str_detect(args$output, ".rds$") ){

  stop(paste(
    "\n  Output file name must be in rds format.",
    "\n  Please change name to have the proper extension (*.rds).\n"
  ))

}

# Check stat file name if provided
if( args$stat != FALSE ){
  if( stringr::str_detect(args$output, ".rds$") ){

    stop(paste(
      "\n  Output file name must be in rds format.",
      "\n  Please change name to have the proper extension (*.rds).\n"
    ))

  }
}


# Print Inputs to terminal
input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(seq_along(args), function(i){
    paste(args[[i]], collapse = ", ")
  })
)

input_table <- input_table[
  match(
    c("dir :", "pattern :", "output :", "stat :"),
    input_table$Variables
  ),
]

cat("\nCombine Multihit Inputs:\n")
print(
  data.frame(input_table),
  right = FALSE, 
  row.names = FALSE
)


# Check for input files
input_files <- list.files(
  path = args$dir, pattern = ifelse(args$pattern == ".", NULL, args$pattern)
)

if( length(input_files) == 0 ){

  cat("\n  Warning:\n  No input files identified, writing empty output files.\n")

  saveRDS(
    object = list(
      "unclustered_multihits" = GenomicRanges::GRanges(),
      "clustered_multihit_positions" = GenomicRanges::GRangesList(),
      "clustered_multihit_lengths" = IRanges::RleList()
    ),
    file = args$output
  )

  if( args$stat != FALSE ){
    write.table(
      x = data.frame(), file = args$stat, 
      sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE
    )
  }

}else{

  cat(paste(
    "\n  A few multihit files to join together:\n ", 
    paste(head(file.path(args$dir, input_files)), collapse = "\n  ")
  ))

}


# Load supporting scripts
source(file.path(code_dir, "supporting_scripts", "printHead.R"))


## Set up stat object
if( args$stat != FALSE ){
  
  sampleName <- unlist(strsplit(args$output, "/"))
  
  sampleName <- unlist(
    strsplit(sampleName[length(sampleName)], ".", fixed = TRUE)
  )[1]
  
  stat <- data.frame(
    sampleName = vector("character"),
    metric = vector("character"),
    count = vector("character")
  )
  
}











# Group and characterize multihits 
# Multihits are reads that align to multiple locations in the reference 
# genome. There are bound to always be a certain proportion of reads aligning
# to repeated sequence due to the high level degree of repeated DNA elements
# within genomes. The final object generated, "multihitData", is a list of 
# three objects. "unclustered_multihits" is a GRanges object where every 
# alignment for every multihit read is present in rows. 
# "clustered_multihit_positions" returns all the possible integration site 
# positions for the multihit. Lastly, "clustered_multihit_lengths" contains the
# length of the templates mapping to the multihit clusters, used for
# abundance calculations.
if( !is.null(args$multihits) ){
  
  unclustered_multihits <- GenomicRanges::GRanges()
  clustered_multihit_positions <- GenomicRanges::GRangesList()
  clustered_multihit_lengths <- list()
  
  if( length(multihit_rpks) > 0 ){
    
    #' Only consider readPairKeys that aligned to multiple genomic loci
    multihit_templates <- coupled_loci[
      coupled_loci$lociPairKey %in% unlist(rpk_loci_key[multihit_rpks])
    ]
    
    multihit_templates <- multihit_templates[
      match(unlist(rpk_loci_key[multihit_rpks]), multihit_templates$lociPairKey)
    ]
    
    multihit_templates$readPairKey <- as.character(S4Vectors::Rle(
      values = multihit_rpks, lengths = lengths(rpk_loci_key[multihit_rpks])
    ))
    
    #' As the loci are expanded from the coupled_loci object, unique templates 
    #' and readPairKeys are present in the readPairKeys unlisted from the 
    #' paired_loci object.
    multihit_keys <- keys[keys$readPairKey %in% multihit_rpks,]
    
    multihit_keys$sampleName <- stringr::str_extract(
      string = as.character(multihit_keys$anchorSeqID), pattern = "^[\\w-]+"
    )
    
    multihit_keys$ID <- multihit_keys$readNames
    
    #' Medians are based on all the potential sites for a given read, which will
    #' be identical for all reads associated with a readPairKey.
    multihit_medians <- round(
      median(GenomicRanges::width(split(
        x = multihit_templates, 
        f = multihit_templates$readPairKey
      )))
    )
    
    multihit_keys$medians <- multihit_medians[multihit_keys$readPairKey]
    
    multihits_pos <- GenomicRanges::flank(
      x = multihit_templates, width = -1, start = TRUE
    )
    
    multihits_red <- GenomicRanges::reduce(
      x = multihits_pos, min.gapwidth = 5L, with.revmap = TRUE
    )  #! Should make min.gapwidth a option
    
    revmap <- multihits_red$revmap
    
    axil_nodes <- as.character(S4Vectors::Rle(
      values = multihit_templates$readPairKey[min(revmap)], 
      lengths = lengths(revmap)
    ))
    
    nodes <- multihit_templates$readPairKey[unlist(revmap)]
    edgelist <- unique(matrix( c(axil_nodes, nodes), ncol = 2 ))
    
    multihits_cluster_data <- igraph::clusters(
      igraph::graph.edgelist(el = edgelist, directed = FALSE)
    )
    
    clus_key <- data.frame(
      row.names = unique(as.character(t(edgelist))),
      "clusID" = multihits_cluster_data$membership
    )
    
    multihits_pos$clusID <- clus_key[multihits_pos$readPairKey, "clusID"]
    multihits_pos <- multihits_pos[order(multihits_pos$clusID)]
    clustered_multihit_index <- as.data.frame(mcols(multihits_pos))
    
    multihit_loci_rle <- S4Vectors::Rle(factor(
      x = clustered_multihit_index$lociPairKey, 
      levels = unique(clustered_multihit_index$lociPairKey)
    ))
    
    multihit_loci_intL <- split(
      multihit_loci_rle, clustered_multihit_index$clusID
    )
    
    clustered_multihit_positions <- GenomicRanges::granges(
      x = multihits_pos[
        match(
          x = unlist(S4Vectors::runValue(multihit_loci_intL)), 
          table = clustered_multihit_index$lociPairKey)
      ]
    )
    
    clustered_multihit_positions <- split(
      x = clustered_multihit_positions,
      f = S4Vectors::Rle(
        values = seq_along(multihit_loci_intL), 
        lengths = S4Vectors::width(S4Vectors::runValue(
          multihit_loci_intL
        )@partitioning)
      )
    )
    
    readPairKey_cluster_index <- unique(
      clustered_multihit_index[,c("readPairKey", "clusID")]
    )
    
    multihit_keys$clusID <- readPairKey_cluster_index$clusID[
      match(multihit_keys$readPairKey, readPairKey_cluster_index$readPairKey)
    ]
    
    multihit_keys <- multihit_keys[order(multihit_keys$medians),]
    
    clustered_multihit_lengths <- split(
      x = S4Vectors::Rle(multihit_keys$medians), 
      f = multihit_keys$clusID
    )
    
    #' Expand the multihit_templates object from readPairKey specific to read
    #' specific.
    multihit_keys <- multihit_keys[order(multihit_keys$readPairKey),]
    
    multihit_readPair_read_exp <- IRanges::IntegerList(
      split(x = seq_len(nrow(multihit_keys)), f = multihit_keys$readPairKey)
    )
    
    unclustered_multihits <- multihit_templates
    
    multihit_readPair_read_exp <- multihit_readPair_read_exp[
      as.character(unclustered_multihits$readPairKey)
    ]
    
    unclustered_multihits <- unclustered_multihits[S4Vectors::Rle(
      values = seq_along(unclustered_multihits),
      lengths = S4Vectors::width(multihit_readPair_read_exp@partitioning)
    )]
    
    names(unclustered_multihits) <- multihit_keys$ID[
      unlist(multihit_readPair_read_exp)
    ]
    
    unclustered_multihits$ID <- multihit_keys$ID[
      unlist(multihit_readPair_read_exp)
    ]
    
    unclustered_multihits$sampleName <- multihit_keys$sampleName[
      unlist(multihit_readPair_read_exp)
    ]
    
  }
  
  stopifnot(
    length(clustered_multihit_positions) == length(clustered_multihit_lengths)
  )
  
  multihitData <- list(
    unclustered_multihits, 
    clustered_multihit_positions, 
    clustered_multihit_lengths
  )
  
  names(multihitData) <- c(
    "unclustered_multihits", 
    "clustered_multihit_positions", 
    "clustered_multihit_lengths"
  )
  
  writeOutputFile(multihitData, file = args$multihits, format = "rds")
  
  printHead(
    data.frame(
      "multihit_reads" = length(unique(names(unclustered_multihits))),
      "multihit_alignments" = length(unique(unclustered_multihits)),
      "multihit_clusters" = length(clustered_multihit_positions),
      "multihit_lengths" = sum(lengths(clustered_multihit_lengths))
    ),
    title = "Multihit metrics", 
    caption = "Metrics highlighting the observation of multiple aligning reads."
  )
  
  if( args$stat != FALSE ){
    
    add_stat <- data.frame(
      sampleName = sampleName,
      metric = c("multihit.reads", "multihit.lengths", "multihit.clusters"),
      count = c(
        length(unique(names(unclustered_multihits))), 
        sum(lengths(clustered_multihit_lengths)), 
        length(clustered_multihit_positions))
    )
    
    stat <- rbind(stat, add_stat)
    
  }
  
}

if( args$stat != FALSE ){

  write.table(
    x = stat, file = args$stat, 
    sep = ",", row.names = FALSE, 
    col.names = FALSE, quote = FALSE
  )
  
}

if( !is.null(args$saveImage) ) save.image(args$saveImage)

q()
