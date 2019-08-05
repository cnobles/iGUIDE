context("Flanking pairs clustering")

input_gr <- GenomicRanges::GRanges(
  seqnames = c("chr4"),
  ranges = IRanges::IRanges(
    start = c(1, 2, 3, 5, 8, 13, 21, 34, 55),
    width = 2
  ),
  strand = c("-", "+", "-", "+", "-", "+", "-", "+", "-")
)


expected_clustering <- c(1, 1, 1, 1, 2, 2, NA, NA, NA)

test_that(
  desc = "Cluster ranges that flank each other in an appropriate manner",
  code = {

    output <- identifyPairedAlgnmts(input_gr, maxgap = 5L, maxovlp = 3L)
    expect_equal(output, expected_clustering)

  }
)
