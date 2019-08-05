context("Expand position ID strings")

input_posids <- c(
  "chr2:+:392-394", "chr2:*:392", "chr2:*:390-391", "chr2:+:400"
)

expected_posids <- list(
  c("chr2:+:392", "chr2:+:393", "chr2:+:394"),
  c("chr2:+:392", "chr2:-:392"),
  c("chr2:+:390", "chr2:+:391", "chr2:-:390", "chr2:-:391"),
  c("chr2:+:400")
)

test_that(
  desc = "Select a number of random sites across the reference genome",
  code = {

    output_posids <- lapply(input_posids, expandPosStr)
    expect_equal(output_posids, expected_posids)

  }
)
