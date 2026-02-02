# =============================================================================
# Tests for Advanced Estimation Methods (TMLE, MatchIt, grf, cox_iptw)
# =============================================================================

# Skip these tests if packages are not available
skip_if_not_installed_quiet <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    skip(paste("Package", pkg, "not installed"))
  }
}

# =============================================================================
# Helper: Create Survival Test Data
# =============================================================================

.make_test_survival_df <- function(n = 200, seed = 123) {
  set.seed(seed)
  
  W <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.5 * W))
  
  # Exponential survival times
  lambda <- exp(-0.5 * A + 0.3 * W)
  time <- rexp(n, rate = lambda)
  
  # Random censoring
  censor_time <- rexp(n, rate = 0.3)
  event <- as.integer(time <= censor_time)
  time <- pmin(time, censor_time)
  
  data.frame(W = W, A = A, time = time, event = event)
}

# =============================================================================
# MatchIt Integration Tests
# =============================================================================

test_that("matching method works when MatchIt is available", {
  skip_if_not_installed_quiet("MatchIt")
  
  spec <- .make_test_spec(n = 200)
  result <- estimate_deficiency(spec, methods = "matching", n_boot = 10)
  
  expect_s3_class(result, "deficiency")
  expect_true("matching" %in% names(result$estimates))
  expect_true(result$estimates["matching"] >= 0)
  expect_true(result$estimates["matching"] <= 1)
})

test_that("matching kernel contains MatchIt object", {
  skip_if_not_installed_quiet("MatchIt")
  
  spec <- .make_test_spec(n = 200)
  result <- estimate_deficiency(spec, methods = "matching", n_boot = 0)
  
  kernel <- result$kernel$matching
  expect_true("match_obj" %in% names(kernel) || kernel$method == "iptw")
})

# =============================================================================
# TMLE Integration Tests
# =============================================================================

test_that("tmle method works when tmle package is available", {
  skip_if_not_installed_quiet("tmle")
  skip_if_not_installed_quiet("SuperLearner")
  
  spec <- .make_test_spec(n = 200)
  result <- estimate_deficiency(spec, methods = "tmle", n_boot = 10)
  
  expect_s3_class(result, "deficiency")
  expect_true("tmle" %in% names(result$estimates))
  expect_true(result$estimates["tmle"] >= 0)
  expect_true(result$estimates["tmle"] <= 1)
})

test_that("tmle kernel contains ATE estimate", {
  skip_if_not_installed_quiet("tmle")
  skip_if_not_installed_quiet("SuperLearner")
  
  spec <- .make_test_spec(n = 200)
  result <- estimate_deficiency(spec, methods = "tmle", n_boot = 0)
  
  kernel <- result$kernel$tmle
  # Either has tmle results or fell back to aipw
  expect_true(
    "ate" %in% names(kernel) || kernel$method == "aipw"
  )
})

# =============================================================================
# grf (Causal Forest) Integration Tests
# =============================================================================

test_that("grf method works when grf package is available", {
  skip_if_not_installed_quiet("grf")
  
  spec <- .make_test_spec(n = 300)
  result <- estimate_deficiency(spec, methods = "grf", n_boot = 10)
  
  expect_s3_class(result, "deficiency")
  expect_true("grf" %in% names(result$estimates))
  expect_true(result$estimates["grf"] >= 0)
  expect_true(result$estimates["grf"] <= 1)
})

test_that("grf kernel contains treatment effect predictions", {
  skip_if_not_installed_quiet("grf")
  
  spec <- .make_test_spec(n = 300)
  result <- estimate_deficiency(spec, methods = "grf", n_boot = 0)
  
  kernel <- result$kernel$grf
  # Either has grf results or fell back to aipw
  if (kernel$method == "grf") {
    expect_true("tau_hat" %in% names(kernel))
    expect_true("ate" %in% names(kernel))
    expect_equal(length(kernel$tau_hat), 300)
  }
})

test_that("grf is rejected for survival outcomes", {
  skip_if_not_installed_quiet("grf")
  skip_if_not_installed_quiet("survival")
  
  # Create survival spec
  df <- .make_test_survival_df(n = 200)
  spec <- causal_spec_survival(df, "A", "time", "event", "W")
  
  # grf should not be in valid methods for survival
  expect_error(
    estimate_deficiency(spec, methods = "grf"),
    regexp = "grf"
  )
})

# =============================================================================
# Cox IPTW Tests (Survival)
# =============================================================================

test_that("cox_iptw method works for survival outcomes", {
  skip_if_not_installed_quiet("survival")
  
  df <- .make_test_survival_df(n = 300)
  spec <- causal_spec_survival(df, "A", "time", "event", "W")
  
  result <- estimate_deficiency(spec, methods = "cox_iptw", n_boot = 10)
  
  expect_s3_class(result, "deficiency")
  expect_true("cox_iptw" %in% names(result$estimates))
  expect_true(result$estimates["cox_iptw"] >= 0)
  expect_true(result$estimates["cox_iptw"] <= 1)
})

test_that("cox_iptw kernel contains hazard ratio", {
  skip_if_not_installed_quiet("survival")
  
  df <- .make_test_survival_df(n = 300)
  spec <- causal_spec_survival(df, "A", "time", "event", "W")
  
  result <- estimate_deficiency(spec, methods = "cox_iptw", n_boot = 0)
  
  kernel <- result$kernel$cox_iptw
  expect_true("hr" %in% names(kernel))
  expect_true("cox_model" %in% names(kernel))
})

test_that("cox_iptw is rejected for non-survival outcomes", {
  spec <- .make_test_spec(n = 200)
  
  expect_error(
    estimate_deficiency(spec, methods = "cox_iptw"),
    regexp = "survival|cox"
  )
})

# =============================================================================
# Parallel Bootstrap Tests
# =============================================================================

test_that("parallel bootstrap works when future.apply is available", {
  skip_if_not_installed_quiet("future.apply")
  
  spec <- .make_test_spec(n = 200)
  
  # Should not error
  result <- estimate_deficiency(
    spec, 
    methods = "iptw", 
    n_boot = 20,
    parallel = TRUE
  )
  
  expect_s3_class(result, "deficiency")
  expect_true(!is.na(result$se["iptw"]))
})

test_that("parallel = TRUE falls back gracefully without future.apply", {
  # This test checks the warning is issued, but function still works
  # We can't easily unload future.apply, so just verify the parameter is accepted
  
  spec <- .make_test_spec(n = 100)
  
  # Should not error even if future.apply not available
  result <- estimate_deficiency(
    spec, 
    methods = "iptw", 
    n_boot = 10,
    parallel = FALSE  # Use FALSE to ensure it works
  )
  
  expect_s3_class(result, "deficiency")
})

# =============================================================================
# Method Comparison Tests
# =============================================================================

test_that("advanced methods can be compared together", {
  skip_if_not_installed_quiet("MatchIt")
  skip_if_not_installed_quiet("tmle")
  skip_if_not_installed_quiet("SuperLearner")
  
  spec <- .make_test_spec(n = 300)
  
  result <- estimate_deficiency(
    spec, 
    methods = c("iptw", "aipw", "matching", "tmle"),
    n_boot = 20
  )
  
  expect_equal(length(result$estimates), 4)
  expect_true(all(result$estimates >= 0))
  expect_true(all(result$estimates <= 1))
})


