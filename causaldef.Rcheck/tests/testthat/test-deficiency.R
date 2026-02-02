# =============================================================================
# Tests for estimate_deficiency()
# =============================================================================

test_that("estimate_deficiency returns correct class", {
  spec <- .make_test_spec(n = 200)
  result <- estimate_deficiency(spec, methods = "iptw", n_boot = 10)
  
  expect_s3_class(result, "deficiency")
  expect_true("estimates" %in% names(result))
  expect_true("se" %in% names(result))
  expect_true("method" %in% names(result))
})

test_that("backdoor kernel achieves δ ≈ 0 when identified", {
  df <- .simulate_backdoor(n = 500, seed = 123)
  spec <- causal_spec(df, "A", "Y", "W")
  result <- estimate_deficiency(spec, methods = "aipw", n_boot = 50)
  
  # Deficiency should be small when backdoor is satisfied
  expect_lt(result$estimates["aipw"], 0.15)
})

test_that("deficiency increases with confounding (Theorem 4.1)", {
  # Test the theoretical formula from Theorem 4.1 directly
  # rather than the SMD proxy which has no signal without covariates
  gammas <- c(0.5, 1.0, 2.0, 4.0)
  alpha <- 1.0
  sigma_A <- 1
  sigma_Y <- 1
  
  # Use the internal deficiency formula
  deltas <- vapply(gammas, function(gamma) {
    causaldef:::.deficiency_gaussian(alpha, gamma, sigma_A, sigma_Y, C = 0.25)
  }, numeric(1))
  
  # Deficiency should monotonically increase with gamma
  expect_true(all(diff(deltas) > 0))
})

test_that("multiple methods can be compared", {
  spec <- .make_test_spec(n = 200)
  result <- estimate_deficiency(
    spec, 
    methods = c("unadjusted", "iptw", "aipw"),
    n_boot = 20
  )
  
  expect_length(result$estimates, 3)
  expect_true(all(c("unadjusted", "iptw", "aipw") %in% names(result$estimates)))
})

test_that("bootstrap confidence intervals are computed", {
  spec <- .make_test_spec(n = 200)
  result <- estimate_deficiency(spec, methods = "iptw", n_boot = 50)
  
  expect_true(!is.null(result$ci))
  expect_equal(nrow(result$ci), 1)
  expect_equal(ncol(result$ci), 2)
  expect_true(result$ci[1, 1] <= result$estimates["iptw"])
  expect_true(result$ci[1, 2] >= result$estimates["iptw"])
})

test_that("error on invalid method", {
  spec <- .make_test_spec()
  expect_error(
    estimate_deficiency(spec, methods = "invalid_method")
  )
})

test_that("deficiency estimates are bounded [0, 1]", {
  spec <- .make_test_spec(n = 200)
  result <- estimate_deficiency(spec, methods = c("unadjusted", "iptw"), n_boot = 0)
  
  expect_true(all(result$estimates >= 0))
  expect_true(all(result$estimates <= 1 + 1e-6))
})

test_that("print.deficiency works", {
  spec <- .make_test_spec(n = 100)
  result <- estimate_deficiency(spec, methods = "iptw", n_boot = 10)
  
  expect_output(print(result), "Deficiency")
})
