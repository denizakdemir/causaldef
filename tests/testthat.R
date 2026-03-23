# Launch tests for causaldef package
# Preload glue so older waldo/testthat stacks on old R register comparison
# methods without tripping delayed S3 registration bugs.
if (requireNamespace("glue", quietly = TRUE)) {
  library(glue)
}
library(testthat)
library(causaldef)

test_check("causaldef")
