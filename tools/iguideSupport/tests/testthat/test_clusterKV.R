context("Key Value (KV) clustering")

input_keys <- c("A", "A", "A", "B", "B", "B", "C", "C", "C")
input_values <- c(1, 2, 3, 3, 4, 5, 6, 7, 8)

expected_clustering <- list("A" = 1, "B" = 1, "C" = 2)

test_that(
  desc = "Cluster keys based on value content",
  code = {

    cluskv_output <- clusterKV(input_keys, input_values, return = "simple")
    expect_equal(as.list(cluskv_output), expected_clustering)

  }
)
