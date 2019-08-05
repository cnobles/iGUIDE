context("Group pileups")

input_gr <- GenomicRanges::GRanges(
  seqnames = c("chr4"),
  ranges = IRanges::IRanges(
    start = c(1, 2, 3, 5, 8, 13, 21),
    width = 5
  ),
  strand = "+"
)

expected_grouping <- c(1, 1, 1, 1, 1, 2, 3)

test_that(
  desc = "Group ranges that overlap by pileup",
  code = {

    output <- groupPileups(gr = input_gr, strand = "+", maxgap = 0L)
    expect_equal(output$clusID, expected_grouping)

  }
)

