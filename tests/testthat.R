library(testthat)
library(gscounts)

if (requireNamespace("gsDesign", quietly = TRUE)) {
  if (requireNamespace("mvtnorm", quietly = TRUE)) {
    test_check("gscounts")
  }
} 

