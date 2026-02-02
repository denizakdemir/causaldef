# =============================================================================
# Tests for Front-Door, Transport, Competing Risks, and IV
# =============================================================================

# =============================================================================
# Front-Door Tests
# =============================================================================

test_that("frontdoor_effect returns correct class", {
  # Simulate front-door data
  set.seed(42)
  n <- 200
  U <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.5 * U))
  M <- 0.5 + 1.2 * A + rnorm(n, sd = 0.5)
  Y <- 1 + 0.8 * M + 0.5 * U + rnorm(n)
  
  df <- data.frame(A = A, M = M, Y = Y)
  spec <- causal_spec(df, "A", "Y", covariates = NULL)
  
  result <- frontdoor_effect(spec, mediator = "M", n_boot = 20)
  
  expect_s3_class(result, "frontdoor_effect")
  expect_true("estimate" %in% names(result))
  expect_true("deficiency" %in% names(result))
  expect_true(is.numeric(result$estimate))
})

test_that("frontdoor_effect detects effect through mediator", {
  set.seed(42)
  n <- 300
  A <- rbinom(n, 1, 0.5)  # Randomized treatment
  M <- 0.5 + 2.0 * A + rnorm(n, sd = 0.3)  # Strong mediation
  Y <- 1 + 1.5 * M + rnorm(n)  # Y only through M
  
  df <- data.frame(A = A, M = M, Y = Y)
  spec <- causal_spec(df, "A", "Y", covariates = NULL)
  
  result <- frontdoor_effect(spec, mediator = "M", n_boot = 0)
  
  # Effect should be approximately 2.0 * 1.5 = 3.0
  expect_true(abs(result$estimate - 3.0) < 1.0)
})

test_that("print.frontdoor_effect works", {
  set.seed(42)
  n <- 100
  A <- rbinom(n, 1, 0.5)
  M <- A + rnorm(n)
  Y <- M + rnorm(n)
  
  df <- data.frame(A = A, M = M, Y = Y)
  spec <- causal_spec(df, "A", "Y", covariates = NULL)
  result <- frontdoor_effect(spec, "M", n_boot = 10)
  
  expect_output(print(result), "Effect estimate")
})

# =============================================================================
# Transport Deficiency Tests
# =============================================================================

test_that("transport_deficiency returns correct class", {
  set.seed(42)
  
  # Source population
  n_s <- 200
  age_s <- rnorm(n_s, 50, 10)
  A_s <- rbinom(n_s, 1, 0.5)
  Y_s <- 10 + 2 * A_s - 0.1 * age_s + rnorm(n_s)
  source_df <- data.frame(age = age_s, A = A_s, Y = Y_s)
  
  # Target population
  n_t <- 100
  age_t <- rnorm(n_t, 60, 8)
  target_df <- data.frame(age = age_t)
  
  source_spec <- causal_spec(source_df, "A", "Y", "age")
  
  result <- transport_deficiency(
    source_spec,
    target_data = target_df,
    transport_vars = "age",
    n_boot = 20
  )
  
  expect_s3_class(result, "transport_deficiency")
  expect_true("delta_transport" %in% names(result))
  expect_true("ate_source" %in% names(result))
  expect_true("ate_target" %in% names(result))
})

test_that("transport_deficiency increases with covariate shift", {
  set.seed(42)
  
  n_s <- 200
  age_s <- rnorm(n_s, 50, 10)
  A_s <- rbinom(n_s, 1, 0.5)
  Y_s <- 10 + 2 * A_s - 0.1 * age_s + rnorm(n_s)
  source_df <- data.frame(age = age_s, A = A_s, Y = Y_s)
  source_spec <- causal_spec(source_df, "A", "Y", "age")
  
  # Small shift
  target_small <- data.frame(age = rnorm(100, 52, 10))
  result_small <- transport_deficiency(source_spec, target_small, n_boot = 0)
  

  # Large shift
  target_large <- data.frame(age = rnorm(100, 75, 5))
  result_large <- transport_deficiency(source_spec, target_large, n_boot = 0)
  
  # Large shift should have higher deficiency
  expect_true(result_large$delta_transport > result_small$delta_transport)
})

test_that("print.transport_deficiency works", {
  set.seed(42)
  
  source_df <- data.frame(
    age = rnorm(100, 50, 10),
    A = rbinom(100, 1, 0.5),
    Y = rnorm(100)
  )
  target_df <- data.frame(age = rnorm(50, 60, 8))
  
  spec <- causal_spec(source_df, "A", "Y", "age")
  result <- transport_deficiency(spec, target_df, n_boot = 0)
  
  expect_output(print(result), "delta_transport")
})

# =============================================================================
# Competing Risks Tests
# =============================================================================

test_that("causal_spec_competing creates valid object", {
  skip_if_not_installed("survival")
  
  set.seed(42)
  n <- 200
  W <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.3 * W))
  
  time_1 <- rexp(n, rate = exp(-0.5 * A))
  time_2 <- rexp(n, rate = exp(0.3 * A))
  time_c <- runif(n, 0, 3)
  
  obs_time <- pmin(time_1, time_2, time_c)
  event <- ifelse(obs_time == time_1, 1,
                  ifelse(obs_time == time_2, 2, 0))
  
  df <- data.frame(W = W, A = A, time = obs_time, event = event)
  
  spec <- causal_spec_competing(df, "A", "time", "event", "W", event_of_interest = 1)
  
  expect_s3_class(spec, "causal_spec_competing")
  expect_s3_class(spec, "causal_spec_survival")
  expect_equal(spec$event_of_interest, 1)
  expect_equal(spec$n_events, 2)
})

test_that("estimate_deficiency_competing returns deficiency object", {
  skip_if_not_installed("survival")
  
  set.seed(42)
  n <- 200
  W <- rnorm(n)
  A <- rbinom(n, 1, 0.5)
  
  time_1 <- rexp(n, rate = exp(-0.3 * A + 0.2 * W))
  time_2 <- rexp(n, rate = exp(0.2 * A - 0.1 * W))
  time_c <- runif(n, 0, 3)
  
  obs_time <- pmin(time_1, time_2, time_c)
  event <- ifelse(obs_time == time_1, 1,
                  ifelse(obs_time == time_2, 2, 0))
  
  df <- data.frame(W = W, A = A, time = obs_time, event = event)
  spec <- causal_spec_competing(df, "A", "time", "event", "W")
  
  result <- estimate_deficiency_competing(spec, method = "cshr", n_boot = 10)
  
  expect_s3_class(result, "deficiency")
  expect_true("estimates" %in% names(result))
  expect_true(result$estimates > 0 && result$estimates <= 1)
})

# =============================================================================
# Instrumental Variable Tests
# =============================================================================

test_that("iv_effect returns correct class", {
  set.seed(42)
  n <- 300
  U <- rnorm(n)
  Z <- rbinom(n, 1, 0.5)
  A <- as.numeric(0.3 + 0.4 * Z + 0.3 * U + rnorm(n, sd = 0.3) > 0.5)
  Y <- 1 + 2 * A + 0.8 * U + rnorm(n)
  
  df <- data.frame(Z = Z, A = A, Y = Y)
  spec <- causal_spec(df, "A", "Y", instrument = "Z")
  
  result <- iv_effect(spec, n_boot = 20)
  
  expect_s3_class(result, "iv_effect")
  expect_true("estimate" %in% names(result))
  expect_true("f_stat" %in% names(result))
  expect_true("weak_iv" %in% names(result))
  expect_true("deficiency" %in% names(result))
})

test_that("iv_effect recovers causal effect with strong instrument", {
  set.seed(42)
  n <- 1000
  U <- rnorm(n)
  Z <- rbinom(n, 1, 0.5)
  # Strong first stage
  A <- as.numeric(0.5 * Z + 0.1 * U + rnorm(n, sd = 0.2) > 0.25)
  # True effect = 3.0
  Y <- 1 + 3.0 * A + 0.5 * U + rnorm(n)
  
  df <- data.frame(Z = Z, A = A, Y = Y)
  spec <- causal_spec(df, "A", "Y", instrument = "Z")
  
  result <- iv_effect(spec, n_boot = 0)
  
  # Should be close to 3.0 with strong instrument
  expect_true(abs(result$estimate - 3.0) < 1.5)
  expect_false(result$weak_iv)
  expect_true(result$f_stat > 10)
})

test_that("iv_effect detects weak instruments", {
  set.seed(42)
  n <- 200
  U <- rnorm(n)
  Z <- rbinom(n, 1, 0.5)
  # Weak first stage (small coefficient on Z)
  A <- as.numeric(0.05 * Z + 0.8 * U + rnorm(n, sd = 0.5) > 0.5)
  Y <- 1 + 2 * A + 0.8 * U + rnorm(n)
  
  df <- data.frame(Z = Z, A = A, Y = Y)
  spec <- causal_spec(df, "A", "Y", instrument = "Z")
  
  result <- iv_effect(spec, n_boot = 0)
  
  # Should detect weak instrument
  expect_true(result$weak_iv || result$f_stat < 15)
  expect_true(result$deficiency > 0.3)
})

test_that("test_instrument provides diagnostic output", {
  set.seed(42)
  n <- 300
  Z <- rbinom(n, 1, 0.5)
  A <- 0.5 * Z + rnorm(n, sd = 0.3)
  Y <- 2 * A + rnorm(n)
  
  df <- data.frame(Z = Z, A = A, Y = Y)
  spec <- causal_spec(df, "A", "Y", instrument = "Z")
  
  result <- test_instrument(spec)
  
  expect_true("relevance" %in% names(result))
  expect_true("reduced_form" %in% names(result))
})

test_that("print.iv_effect works", {
  set.seed(42)
  n <- 200
  Z <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, 0.3 + 0.4 * Z)
  Y <- 1 + 2 * A + rnorm(n)
  
  df <- data.frame(Z = Z, A = A, Y = Y)
  spec <- causal_spec(df, "A", "Y", instrument = "Z")
  result <- iv_effect(spec, n_boot = 10)
  
  expect_output(print(result), "F-statistic")
})
