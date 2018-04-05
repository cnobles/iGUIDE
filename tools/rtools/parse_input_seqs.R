seq_file <- ShortRead::readFastq(snakemake@input[[1]])
seq_file <- split(
  seq_file, ceiling(
    seq_along(seq_file)/(length(seq_file)/snakemake@params[[1]])))
input_path <- unlist(stringr::str_split(snakemake@input[[1]], "/"))
output_path <- paste(input_path[1:(length(input_path) - 1)], collapse = "/")
sample_name <- unlist(
  stringr::str_split(tail(input_path, n = 1), stringr::fixed(".")))
sample_name <- paste(sample_name[1:2], collapse = ".")
null <- lapply(seq_along(seq_file), function(i){
  sam_name <- paste0(
    sample_name, ".chunk", stringr::str_pad(i, 2, pad = "0"), ".fastq.gz")
  ShortRead::writeFastq(
    object = seq_file[[i]], 
    file = file.path(output_path, sam_name), 
    compress = TRUE)
})
