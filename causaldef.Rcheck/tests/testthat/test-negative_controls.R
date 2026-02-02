# =============================================================================
# Tests for nc_diagnostic()
# =============================================================================

test_that("nc_diagnostic returns correct class", {
  df <- .simulate_with_nc(200, seed = 1)
  spec <- causal_spec(df, "A", "Y", "W", negative_control = "Y_nc")
  result <- nc_diagnostic(spec, method = "iptw", n_boot = 20)
  
  expect_s3_class(result, "nc_diagnostic")
  expect_true("delta_nc" %in% names(result))
  expect_true("delta_bound" %in% names(result))
  expect_true("falsified" %in% names(result))
  expect_true("p_value" %in% names(result))
})

test_that("nc_diagnostic requires negative control", {
  df <- .simulate_backdoor(100, seed = 1)
  spec <- causal_spec(df, "A", "Y", "W")  # No negative control
  
  expect_error(
    nc_diagnostic(spec, method = "iptw"),
    "negative control"
  )
})

test_that("nc_diagnostic detects confounding violation", {
  # Strong confounding should be detectable
  df <- .simulate_with_nc(500, alpha = 2, gamma_y = 2, gamma_nc = 2, seed = 42)
  spec <- causal_spec(df, "A", "Y", "W", negative_control = "Y_nc")
  
  # Without proper adjustment, should show high NC deficiency
  result_unadj <- nc_diagnostic(spec, method = "unadjusted", n_boot = 50)
  
  # With adjustment, should show lower NC deficiency
  result_adj <- nc_diagnostic(spec, method = "iptw", n_boot = 50)
  
  # Adjusted should have lower NC deficiency
  expect_lt(result_adj$delta_nc, result_unadj$delta_nc)
})

test_that("nc_diagnostic respects kappa parameter", {
  df <- .simulate_with_nc(200, seed = 1)
  spec <- causal_spec(df, "A", "Y", "W", negative_control = "Y_nc")
  
  result_k1 <- nc_diagnostic(spec, method = "iptw", kappa = 1.0, n_boot = 20)
  result_k2 <- nc_diagnostic(spec, method = "iptw", kappa = 2.0, n_boot = 20)
  
  # delta_bound = kappa * delta_nc
  expect_equal(result_k2$delta_bound, 2 * result_k1$delta_bound, tolerance = 0.01)
})

test_that("nc_diagnostic falsification flag works", {
  df <- .simulate_with_nc(200, seed = 1)
  spec <- causal_spec(df, "A", "Y", "W", negative_control = "Y_nc")
  
  result <- nc_diagnostic(spec, method = "iptw", alpha = 0.05, n_boot = 50)
  
  expect_type(result$falsified, "logical")
})

test_that("print.nc_diagnostic works", {
  df <- .simulate_with_nc(100, seed = 1)
  spec <- causal_spec(df, "A", "Y", "W", negative_control = "Y_nc")
  result <- nc_diagnostic(spec, method = "iptw", n_boot = 10)
  
  expect_output(print(result), "Negative Control")
})
