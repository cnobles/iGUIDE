context("Divergent sequences")

input_ref <- "ABCDEFG"
input_seqs <- c("GBCDEFA", "ACBDEFG", "ABCDFEG")
expected_output <- c("G.....A", ".CB....", "....FE.")

test_that(
  desc = "Indicate diverging sequences within a character string",
  code = {

    output <- divSeq(input_seqs, input_ref)
    expect_equal(output, expected_output)

  }
)
