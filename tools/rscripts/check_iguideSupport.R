# This script executes unit tests for the iguideSupport R-package

tto <- devtools::test(pkg = "tools/iguideSupport")

num_success <- sum(
  sapply(seq_along(tto), function(i){
    tto[[i]]$results[[1]]$message}
  ) == "success"
)

num_failed <- length(tto) - num_success

q(save = "no", status = num_failed)
