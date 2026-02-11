# =============================================================================
# Test Helpers: Simulation Functions
# =============================================================================

#' Simulate backdoor-identified data
#' @noRd
.simulate_backdoor <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  W <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.5 * W))
  Y <- 1 + 2 * A + W + rnorm(n)
  data.frame(W = W, A = A, Y = Y)
}

#' Simulate confounded data (linear Gaussian confounding setup)
#' @noRd
.simulate_confounded <- function(n, alpha = 1, gamma = 1, beta = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  U <- rnorm(n)
  A <- alpha * U + rnorm(n)
  Y <- beta * A + gamma * U + rnorm(n)
  data.frame(A = A, Y = Y, U = U)
}

#' Simulate data with negative control
#' @noRd
.simulate_with_nc <- function(n, alpha = 1, gamma_y = 1, gamma_nc = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  U <- rnorm(n)
  W <- U + rnorm(n, sd = 0.5)
  # Binary treatment for compatibility with glm logistic
  A <- rbinom(n, 1, plogis(alpha * U + 0.5 * W))
  Y <- 1 * A + gamma_y * U + rnorm(n)
  Y_nc <- gamma_nc * U + rnorm(n)  # No causal effect from A
  data.frame(W = W, A = A, Y = Y, Y_nc = Y_nc)
}

#' Simulate survival data
#' @noRd
.simulate_survival <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  W <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.5 * W))
  lambda <- exp(-1 + 0.5 * A + 0.3 * W)
  time <- rexp(n, rate = lambda)
  censor_time <- rexp(n, rate = 0.1)
  observed_time <- pmin(time, censor_time)
  event <- as.integer(time <= censor_time)
  data.frame(W = W, A = A, time = observed_time, event = event)
}

#' Create minimal test spec
#' @noRd
.make_test_spec <- function(n = 100, seed = 1) {
  df <- .simulate_backdoor(n, seed)
  causal_spec(df, "A", "Y", "W")
}

#' Create minimal test survival spec
#' @noRd
.make_survival_spec <- function(n = 100, seed = 1) {
  df <- .simulate_survival(n, seed)
  causal_spec_survival(df, "A", "time", "event", "W")
}
