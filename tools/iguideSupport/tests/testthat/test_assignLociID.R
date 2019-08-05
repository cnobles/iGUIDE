context("Assign loci ID")

input_gr <- GenomicRanges::GRanges(
  seqnames = c("chr4"),
  ranges = IRanges::IRanges(
    start = c(1, 2, 3, 5, 8, 13, 21, 34, 55),
    width = 2
  ),
  strand = c("-", "+", "-", "+", "-", "+", "-", "+", "-")
)

expected_lociIDs <- c(1L, 1L, 1L, 1L, 1L, 1L, 3L, 2L, 4L)

test_that(
  desc = "Assign locational IDs based on clustering algorithms",
  code = {

    output <- assignLociID(input_gr, pairgap = 5L)
    expect_equal(output, expected_lociIDs)

  }
)
