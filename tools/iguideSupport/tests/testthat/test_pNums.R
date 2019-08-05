context("Pretty numbers")

nums <- c(1, 10, 1000, 1.0, 1.001, 0.0001)
expected_nums <- c("1", "10", "1,000", "1", "1.001", "1e-04")

test_that(
  desc = "Pretty formating for numbers in reports",
  code = {

    tested_nums <- unlist(sapply(nums, pNums))
    expect_equal(tested_nums, expected_nums)

  }
)
