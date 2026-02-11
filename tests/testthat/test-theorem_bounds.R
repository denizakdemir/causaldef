# =============================================================================
# Tests for theorem-aligned utilities (Results 1, 3, 4, 6, 8)
# =============================================================================

test_that("rkhs_rate_bound decreases with n", {
  b1 <- rkhs_rate_bound(n = 200, beta = 1, d_w = 1, eta = 0.2, xi = 0.05)
  b2 <- rkhs_rate_bound(n = 800, beta = 1, d_w = 1, eta = 0.2, xi = 0.05)
  expect_true(is.finite(b1) && is.finite(b2))
  expect_gt(b1, b2)
})

test_that("wasserstein_deficiency_gaussian matches closed form", {
  delta <- wasserstein_deficiency_gaussian(alpha = 1, gamma = 2, sigma_A = 3, a = 1.5)
  expect_equal(delta, 0.3, tolerance = 1e-12)
})

test_that("sharp_lower_bound returns matching TV bounds", {
  out <- sharp_lower_bound(alpha = 1.2, gamma = 0.7, sigma_A = 1.0, sigma_Y = 1.0, metric = "tv")
  expect_equal(out$lower, out$upper, tolerance = 1e-12)
  expect_true(out$lower >= 0 && out$lower <= 1)
})

test_that("partial_id_set uses correct half-widths", {
  mean_set <- partial_id_set(estimate = 0.5, delta = 0.1, estimand = "mean", outcome_range = c(0, 1))
  expect_equal(mean_set$half_width, 0.2, tolerance = 1e-12)
  expect_equal(mean_set$lower, 0.3, tolerance = 1e-12)
  expect_equal(mean_set$upper, 0.7, tolerance = 1e-12)

  ate_set <- partial_id_set(estimate = 0.0, delta = 0.1, estimand = "ate", outcome_range = c(0, 1))
  expect_equal(ate_set$half_width, 0.4, tolerance = 1e-12)
})

test_that("policy_regret_bound_vc adds a complexity penalty", {
  b <- policy_regret_bound_vc(deficiency = 0.1, vc_dim = 2, n = 1000, xi = 0.05,
                              utility_range = c(0, 1), obs_regret = 0.2, C = 2)
  expect_s3_class(b, "policy_bound")
  expect_true(!is.null(b$complexity_penalty))
  expect_true(b$complexity_penalty > 0)
  expect_equal(b$transfer_penalty, 0.1, tolerance = 1e-12)
  expect_equal(b$minimax_floor, 0.05, tolerance = 1e-12)
  expect_equal(b$regret_bound, 0.2 + b$transfer_penalty + b$complexity_penalty, tolerance = 1e-12)
})

test_that("overlap_diagnostic returns propensity summaries and ESS", {
  set.seed(1)
  n <- 300
  W <- rnorm(n)
  A <- rbinom(n, 1, plogis(1.2 * W))
  Y <- A + W + rnorm(n)
  df <- data.frame(W = W, A = A, Y = Y)
  spec <- causal_spec(df, "A", "Y", "W")

  od <- overlap_diagnostic(spec, trim = 0.05)
  expect_s3_class(od, "overlap_diagnostic")
  expect_true(is.numeric(od$ess_iptw) && od$ess_iptw > 0)
  expect_true(length(od$ps_summary) >= 5)
})

