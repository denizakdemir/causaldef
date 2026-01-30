# =============================================================================
# Tests for causal_spec()
# =============================================================================

test_that("causal_spec creates valid object", {
  df <- .simulate_backdoor(100, seed = 1)
  spec <- causal_spec(df, "A", "Y", "W")
  
  expect_s3_class(spec, "causal_spec")
  expect_equal(spec$treatment, "A")
  expect_equal(spec$outcome, "Y")
  expect_equal(spec$covariates, "W")
  expect_equal(spec$n, 100)
})

test_that("causal_spec validates inputs", {
  df <- .simulate_backdoor(100, seed = 1)
  
  # Missing treatment variable

  expect_error(
    causal_spec(df, "MISSING", "Y", "W")
  )
  
  # Missing outcome variable
  expect_error(
    causal_spec(df, "A", "MISSING", "W")
  )
  
  # Empty data
  expect_error(
    causal_spec(data.frame(), "A", "Y")
  )
})

test_that("causal_spec handles different estimands", {
  df <- .simulate_backdoor(100, seed = 1)
  
  spec_ate <- causal_spec(df, "A", "Y", "W", estimand = "ATE")
  expect_equal(spec_ate$estimand, "ATE")
  
  spec_att <- causal_spec(df, "A", "Y", "W", estimand = "ATT")
  expect_equal(spec_att$estimand, "ATT")
  
  spec_atc <- causal_spec(df, "A", "Y", "W", estimand = "ATC")
  expect_equal(spec_atc$estimand, "ATC")
})

test_that("causal_spec infers treatment type", {
  # Binary treatment
  df_binary <- .simulate_backdoor(100, seed = 1)
  spec_binary <- causal_spec(df_binary, "A", "Y", "W")
  expect_equal(spec_binary$treatment_type, "binary")
  
  # Continuous treatment
  df_cont <- data.frame(
    W = rnorm(100),
    A = rnorm(100),
    Y = rnorm(100)
  )
  spec_cont <- causal_spec(df_cont, "A", "Y", "W")
  expect_equal(spec_cont$treatment_type, "continuous")
})

test_that("causal_spec handles negative controls", {
  df <- .simulate_with_nc(100, seed = 1)
  spec <- causal_spec(df, "A", "Y", "W", negative_control = "Y_nc")
  
  expect_equal(spec$negative_control, "Y_nc")
})

test_that("causal_spec handles missing values", {
  df <- .simulate_backdoor(100, seed = 1)
  df$Y[1:5] <- NA  # Introduce missing values
  
  expect_warning(
    spec <- causal_spec(df, "A", "Y", "W"),
    "dropped"
  )
  expect_equal(spec$n, 95)
})

test_that("causal_spec handles outcome types", {
  df <- .simulate_backdoor(100, seed = 1)
  
  spec_cont <- causal_spec(df, "A", "Y", "W", outcome_type = "continuous")
  expect_equal(spec_cont$outcome_type, "continuous")
  
  # Binary outcome
  df$Y_bin <- rbinom(100, 1, 0.5)
  spec_bin <- causal_spec(df, "A", "Y_bin", "W", outcome_type = "binary")
  expect_equal(spec_bin$outcome_type, "binary")
})

test_that("print.causal_spec works", {
  df <- .simulate_backdoor(100, seed = 1)
  spec <- causal_spec(df, "A", "Y", "W")
  
  expect_output(print(spec), "Causal Specification")
  expect_output(print(spec), "Treatment")
  expect_output(print(spec), "Outcome")
})
