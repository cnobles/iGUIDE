context("Pileup clustering")

input_gr <- GenomicRanges::GRanges(
  seqnames = c("chr4"),
  ranges = IRanges::IRanges(
    start = c(1, 2, 3, 5, 8, 13, 21),
    width = 5
  ),
  strand = "+"
)

expected_clustering <- c(rep("chr4:+:1", 5), "chr4:+:13", "chr4:+:21")

test_that(
  desc = "Cluster identified ranges by pileup method",
  code = {

    output <- pileupCluster(input_gr, return = "simple")
    expect_equal(output, expected_clustering)

  }
)
