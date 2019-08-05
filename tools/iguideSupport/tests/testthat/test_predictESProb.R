context("Predict edit site probability")

input_vals <- c(-1, 1, 2, 3, 5, 8, 13, 21, 34, 55)
test_vals <- c(rep(0, 10), rep(1, 5), rep(2, 3), 3, 3, 4, 6, 9, 15, 27)
test_density <- density(test_vals)

expected_output <- c(0.96, 0.64, 0.44, 0.3, 0.17, 0.12, 0.08, 0.04, NA, NA)

test_that(
  desc = "Predict edit site probability given a distribution",
  code = {

    output <- predictESProb(input_vals, test_density)
    expect_equal(round(output, digits = 2), expected_output)

  }
)
