test_that("causal_spec_survival coerces factor event to 0/1", {
  skip_if_not_installed("survival")
  skip_if(getRversion() < "4.0", "Survival package requires R >= 4.0 for internal deparse1 calls")
  
  set.seed(1)
  n <- 50
  df <- data.frame(
    W = rnorm(n),
    A = rbinom(n, 1, 0.5),
    time = rexp(n),
    event_status = factor(sample(c("Death", "Censored"), n, replace = TRUE))
  )
  
  spec <- causal_spec_survival(
    data = df,
    treatment = "A",
    time = "time",
    event = "event_status",
    covariates = "W",
    estimand = "RMST",
    horizon = 1
  )
  
  expect_true(all(spec$data$event_status %in% c(0, 1)))
  expect_equal(spec$n_events, sum(spec$data$event_status == 1))
})

test_that("causal_spec_survival warns on multi-state factor events", {
  skip_if_not_installed("survival")
  skip_if(getRversion() < "4.0", "Survival package requires R >= 4.0 for internal deparse1 calls")
  
  set.seed(1)
  n <- 50
  df <- data.frame(
    W = rnorm(n),
    A = rbinom(n, 1, 0.5),
    time = rexp(n),
    event_status = factor(sample(c("Death", "Relapse", "Censored"), n, replace = TRUE))
  )
  
  expect_warning(
    causal_spec_survival(
      data = df,
      treatment = "A",
      time = "time",
      event = "event_status",
      covariates = "W",
      estimand = "RMST",
      horizon = 1
    ),
    "multiple non-censor levels"
  )
})

