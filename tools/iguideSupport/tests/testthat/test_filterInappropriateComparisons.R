context("Filter inappropriate comparisons")

gRNA_match <- c("B2M", "TRAC", "TRBC", "B2M")

specimen <- c("A", "A", "B", "C")

treatment <- list(
  "A" = c("B2M"),
  "B" = c("TRAC", "TRBC"),
  "C" = c("B2M", "TRAC", "TRBC")
)

expected_output <- c("B2M", "No_valid_match", "TRBC", "B2M")

test_that(
  desc = "Filter inappropriate comparisons made during analysis",
  code = {

    output <- filterInappropriateComparisons(gRNA_match, specimen, treatment)
    expect_equal(output, expected_output)

  }
)
