# This script executes unit tests for the iguideSupport R-package

tto <- devtools::test(pkg = "tools/iguideSupport")

tto_df <- as.data.frame(tto)

num_failed <- sum(tto_df$failed)

q(save = "no", status = num_failed)
