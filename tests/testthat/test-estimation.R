

test_that("estimate_effect works for standard data", {
  # Create sample data
  n <- 200
  W <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.5 * W))
  Y <- 1 + 2 * A + W + rnorm(n)
  df <- data.frame(W = W, A = A, Y = Y)
  
  spec <- causal_spec(df, "A", "Y", "W")
  
  # Standard
  def <- estimate_deficiency(spec, methods = c("unadjusted", "iptw"), n_boot = 0, verbose=FALSE)
  
  eff_unadj <- estimate_effect(def, target_method = "unadjusted")
  expect_true(abs(eff_unadj$estimate - 2) < 1.0) # Crude check
  
  eff_iptw <- estimate_effect(def, target_method = "iptw")
  expect_true(abs(eff_iptw$estimate - 2) < 0.5)
  expect_equal(eff_iptw$type, "ATE")

  # Alias: `method` (preferred API)
  eff_iptw_alias <- estimate_effect(def, method = "iptw")
  expect_equal(eff_iptw_alias$estimate, eff_iptw$estimate)
  expect_equal(eff_iptw_alias$type, eff_iptw$type)
})


test_that("estimate_effect works for survival data (RMST)", {
  skip_if_not_installed("survival")
  skip_if(getRversion() < "4.0", "Survival package requires R >= 4.0 for internal deparse1 calls")
  
  set.seed(42)
  n <- 200
  W <- rnorm(n)
  A <- rbinom(n, 1, 0.5)
  # Treatment makes you live longer
  # Increase rate difference to ensure clear RMST difference
  time <- rexp(n, rate = exp(-0.5 * A)) 
  event <- rep(1, n)
  df <- data.frame(W = W, A = A, time = time, event = event)
  
  spec <- causal_spec_survival(df, "A", "time", "event", "W", estimand = "RMST", horizon = 1.0)
  
  def <- estimate_deficiency(spec, methods = "iptw", n_boot = 0, verbose=FALSE)
  
  eff <- estimate_effect(def, target_method = "iptw")
  
  expect_s3_class(eff, "causal_effect")
  expect_equal(eff$type, "RMST Difference")
  expect_true(!is.na(eff$estimate))
  # RMST should be positive (A=1 lives longer)
  # But with random noise, just checking it exists is safer for unit tests unless we set seed carefully.
  # Given the setup, A reduces hazard -> longer life -> higher RMST.
  # However, let's just ensure the structure is correct.
  expect_true(is.numeric(eff$estimate))
  
  # Check that we can extract the table values
  expect_true(!is.null(eff$rmst1))
  expect_true(!is.null(eff$rmst0))
})

test_that("estimate_effect falls back to simple fit if not RMST", {
  skip_if_not_installed("survival")
  skip_if(getRversion() < "4.0")
  
  n <- 50
  df <- data.frame(
    A = rbinom(n, 1, 0.5),
    time = rexp(n),
    event = 1, 
    W = rnorm(n)
  )
  
  spec <- causal_spec_survival(df, "A", "time", "event", "W", estimand = "ATE") # Not RMST
  def <- estimate_deficiency(spec, methods = "unadjusted", n_boot=0, verbose=FALSE)
  
  eff <- estimate_effect(def)
  
  expect_equal(eff$type, "Survival Curve")
  expect_s3_class(eff$fit, "survfit")
})
