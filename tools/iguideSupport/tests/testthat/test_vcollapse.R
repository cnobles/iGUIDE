context("Vector collapse")

input_df <- data.frame(
  "A" = letters[1:5],
  "B" = 3:7,
  "C" = LETTERS[2:6],
  stringsAsFactors = FALSE
)

expected_output <- c("a-3-B", "b-4-C", "c-5-D", "d-6-E", "e-7-F")

test_that(
  desc = "Collapse multiple columns into a character string vector",
  code = {

    output <- vcollapse(input_df, "-")
    expect_equal(output, expected_output)

  }
)
