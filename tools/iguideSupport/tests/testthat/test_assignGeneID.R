context("Assign gene ID")

library(GenomicRanges)

input_gr <- GenomicRanges::GRanges(
  seqnames = c("chr4"),
  ranges = IRanges::IRanges(
    start = c(1, 2, 3, 5, 8, 13, 21, 34, 55),
    width = 2
  ),
  strand = c("-", "+", "-", "+", "-", "+", "-", "+", "-")
)

ref_range <- GenomicRanges::GRanges(
  seqnames = c("chr4"),
  ranges = IRanges::IRanges(
    start = c(2, 5, 13, 20, 32, 60),
    width = 3
  ),
  strand = c("-", "+", "-", "+", "-", "+"),
  annot_sym = LETTERS[1:6]
)

sim_onco <- LETTERS[3:6]
sim_special <- LETTERS[5:6]

expected_geneIDs <- c("A", "A*", "A*", "B*", "B", "C*~", "D*~", "E*~!", "F~!")

test_that(
  desc = "Assign Gene IDs based on reference lists",
  code = {

    output <- assignGeneID(
      seqnames = GenomicRanges::seqnames(input_gr),
      positions = GenomicRanges::start(input_gr),
      reference = input_gr,
      ref.genes = ref_range,
      onco.genes = sim_onco,
      special.genes = sim_special,
      annotations = TRUE
    )

    expect_equal(output, expected_geneIDs)

  }
)
