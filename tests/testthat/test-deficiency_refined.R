
test_that("TV distance deficiency behaves correctly", {
  # Case 1: Randomized Trial (Perfect Overlap)
  # Deficiency should be close to 0 (allow for sampling noise)
  set.seed(42)
  n <- 500
  W <- rnorm(n)
  A <- rbinom(n, 1, 0.5) # Randomized
  Y <- A + W + rnorm(n)
  df <- data.frame(W = W, A = A, Y = Y)
  spec <- causal_spec(df, "A", "Y", "W")
  
  res <- estimate_deficiency(spec, methods = "unadjusted", n_boot = 0, verbose = FALSE)
  expect_lt(res$estimates[["unadjusted"]], 0.1) 
  
  # Case 2: Strong Confounding (Poor Overlap)
  # Deficiency should be high for unadjusted
  set.seed(42)
  W <- rnorm(n)
  # Propensity is sigmoid(3*W). High W -> A=1, Low W -> A=0
  A <- rbinom(n, 1, plogis(3 * W)) 
  Y <- A + W + rnorm(n)
  df_bad <- data.frame(W = W, A = A, Y = Y)
  spec_bad <- suppressWarnings(causal_spec(df_bad, "A", "Y", "W"))
  
  res_bad <- estimate_deficiency(spec_bad, methods = c("unadjusted", "iptw"), n_boot = 0, verbose = FALSE)
  
  # Unadjusted should be high (confounded)
  expect_gt(res_bad$estimates[["unadjusted"]], 0.2)
  
  # IPTW should fix the overlap (lower deficiency)
  # Note: With extreme weights, IPTW might be unstable, but deficiency (TV of weighted PS)
  # should functionally be lower than the raw mismatch.
  expect_lt(res_bad$estimates[["iptw"]], res_bad$estimates[["unadjusted"]])
})

test_that("TV calculation handles only-intercept models gracefully", {
  set.seed(42)
  n <- 100
  A <- rbinom(n, 1, 0.5)
  Y <- rnorm(n)
  df <- data.frame(A = A, Y = Y) 
  # No covariates!
  spec <- causal_spec(df, "A", "Y", NULL)
  
  res <- estimate_deficiency(spec, methods = "unadjusted", n_boot = 0, verbose = FALSE)
  
  # With no covariates, propensity is constant.
  # TV between two point masses?
  # If p(x) and q(x) are identical (constant), TV is 0.
  expect_lt(res$estimates[["unadjusted"]], 0.05)
})
