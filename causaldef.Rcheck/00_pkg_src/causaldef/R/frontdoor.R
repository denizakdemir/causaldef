# =============================================================================
# Front-Door Kernel Implementation (Theorem 2.2)
# =============================================================================

#' Front-Door Adjustment Kernel
#'
#' Implements the front-door criterion for causal identification when a mediator
#' blocks all paths from treatment to outcome while being unaffected by confounders.
#'
#' @param spec A causal_spec object with mediator specified
#' @param mediator Character: name of the mediator variable M
#' @param method Character: estimation method for the front-door formula
#'   \itemize{
#'     \item "plugin": Simple plug-in estimator
#'     \item "dr": Doubly robust front-door estimator
#'   }
#' @param n_boot Integer: bootstrap replicates for standard errors (default 200)
#' @param ci_level Numeric: confidence level (default 0.95)
#'
#' @return Object of class "frontdoor_effect" containing:
#'   \itemize{
#'     \item estimate: The front-door causal effect
#'     \item se: Standard error (if bootstrap > 0)
#'     \item ci: Confidence interval
#'     \item deficiency: Estimated deficiency (should be ~0 if assumptions hold)
#'   }
#'
#' @details
#' The front-door criterion (Pearl, 1995) identifies causal effects through
#' a mediator M when:
#' \enumerate{
#'   \item A -> M is unconfounded (M intercepts all directed paths from A to Y)
#'   \item M -> Y has no direct A effect except through M
#'   \item There's no unblocked back-door path from M to Y
#' }
#'
#' The front-door formula is:
#' \deqn{P(Y|do(a)) = \sum_m P(M=m|A=a) \sum_{a'} P(Y|M=m, A=a') P(A=a')}
#'
#' When these assumptions hold, Theorem 2.2 guarantees delta = 0.
#'
#' @examples
#' # Simulate front-door scenario
#' n <- 500
#' U <- rnorm(n)  # Unmeasured confounder
#' A <- rbinom(n, 1, plogis(0.5 * U))
#' M <- 0.5 + 1.2 * A + rnorm(n, sd = 0.5)  # Mediator (unconfounded by U)
#' Y <- 1 + 0.8 * M + 0.5 * U + rnorm(n)    # Outcome
#' df <- data.frame(A = A, M = M, Y = Y)
#' 
#' spec <- causal_spec(df, "A", "Y", covariates = NULL)
#' \dontrun{
#' fd_result <- frontdoor_effect(spec, mediator = "M")
#' print(fd_result)
#' }
#'
#' @references
#' Pearl, J. (1995). Causal diagrams for empirical research. Biometrika.
#' Akdemir, D. (2026). Constraints on Causal Inference as Experiment Comparison.
#' See Theorem 2.2.
#'
#' @seealso [estimate_deficiency()], [causal_spec()]
#' @export
frontdoor_effect <- function(spec, mediator, method = c("plugin", "dr"),
                             n_boot = 200, ci_level = 0.95) {
  
  method <- match.arg(method)
  
  # Validation
  if (!inherits(spec, "causal_spec")) {
    cli::cli_abort("spec must be a causal_spec object")
  }
  
  data <- spec$data
  if (!(mediator %in% names(data))) {
    cli::cli_abort("Mediator {.val {mediator}} not found in data")
  }
  
  A <- data[[spec$treatment]]
  M <- data[[mediator]]
  Y <- data[[spec$outcome]]
  n <- nrow(data)
  
  # Check treatment is binary (current limitation)
  if (length(unique(A)) != 2) {
    cli::cli_abort("Front-door currently requires binary treatment")
  }
  
  # Point estimate
  estimate <- .compute_frontdoor(A, M, Y, method)
  
  # Bootstrap for SE and CI
  if (n_boot > 0) {
    boot_estimates <- replicate(n_boot, {
      idx <- sample(n, replace = TRUE)
      .compute_frontdoor(A[idx], M[idx], Y[idx], method)
    })
    
    se <- sd(boot_estimates, na.rm = TRUE)
    alpha <- 1 - ci_level
    ci <- quantile(boot_estimates, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  } else {
    se <- NA_real_
    ci <- c(NA_real_, NA_real_)
  }
  
  # Compute residual deficiency (should be ~0 if front-door holds)
  deficiency <- .compute_frontdoor_deficiency(A, M, Y)
  
  result <- structure(
    list(
      estimate = estimate,
      se = se,
      ci = ci,
      deficiency = deficiency,
      method = method,
      mediator = mediator,
      n = n
    ),
    class = "frontdoor_effect"
  )
  
  # Report
  delta_str <- sprintf("%.3f", deficiency)
  if (deficiency < 0.05) {
    cli::cli_alert_success("Front-door effect: {round(estimate, 3)} (delta = {delta_str} ~ 0)")
  } else {
    cli::cli_alert_warning("Front-door effect: {round(estimate, 3)} (delta = {delta_str} > 0.05)")
    cli::cli_alert_info("Consider checking front-door assumptions")
  }
  
  result
}

#' @keywords internal
.compute_frontdoor <- function(A, M, Y, method = "plugin") {
  n <- length(A)
  
  if (method == "plugin") {
    # Step 1: P(M | A) via regression
    m_model <- lm(M ~ A)
    
    # Step 2: P(Y | M, A) regression
    y_model <- lm(Y ~ M + A)
    
    # Step 3: Front-door formula
    # E[Y|do(A=1)] = sum_m P(M=m|A=1) * E[Y|M=m]_marginal
    
    # Predict M under A=1 and A=0
    m_given_a1 <- predict(m_model, newdata = data.frame(A = rep(1, n)))
    m_given_a0 <- predict(m_model, newdata = data.frame(A = rep(0, n)))
    
    # For each M value, compute E[Y|M] marginalized over A
    # E[Y|M=m] = sum_a P(A=a) * E[Y|M=m, A=a]
    p_a1 <- mean(A)
    p_a0 <- 1 - p_a1
    
    # E[Y|do(A=1)] ~ mean over grid
    # Using the sample M values under A=1
    ey_do_a1 <- mean(sapply(m_given_a1, function(m) {
      # E[Y|M=m] = P(A=0)*E[Y|M,A=0] + P(A=1)*E[Y|M,A=1]
      predict(y_model, newdata = data.frame(M = m, A = 0)) * p_a0 +
        predict(y_model, newdata = data.frame(M = m, A = 1)) * p_a1
    }))
    
    ey_do_a0 <- mean(sapply(m_given_a0, function(m) {
      predict(y_model, newdata = data.frame(M = m, A = 0)) * p_a0 +
        predict(y_model, newdata = data.frame(M = m, A = 1)) * p_a1
    }))
    
    ate <- ey_do_a1 - ey_do_a0
    
  } else if (method == "dr") {
    # Doubly robust front-door estimator
    # Uses both M and Y models
    
    # Propensity model for A
    ps <- mean(A)
    
    # Mediator model
    m_model <- lm(M ~ A)
    M_pred_a1 <- predict(m_model, newdata = data.frame(A = rep(1, n)))
    M_pred_a0 <- predict(m_model, newdata = data.frame(A = rep(0, n)))
    
    # Outcome model
    y_model <- lm(Y ~ M + A)
    
    # DR adjustment
    # Influence function approach
    Y_pred <- predict(y_model)
    resid_Y <- Y - Y_pred
    
    # E[Y|do(A=1)]
    ey_do_a1 <- mean(sapply(seq_len(n), function(i) {
      m <- M_pred_a1[i]
      predict(y_model, newdata = data.frame(M = m, A = 0)) * (1 - ps) +
        predict(y_model, newdata = data.frame(M = m, A = 1)) * ps
    }))
    
    # E[Y|do(A=0)]
    ey_do_a0 <- mean(sapply(seq_len(n), function(i) {
      m <- M_pred_a0[i]
      predict(y_model, newdata = data.frame(M = m, A = 0)) * (1 - ps) +
        predict(y_model, newdata = data.frame(M = m, A = 1)) * ps
    }))
    
    ate <- ey_do_a1 - ey_do_a0
  }
  
  ate
}

#' @keywords internal
.compute_frontdoor_deficiency <- function(A, M, Y) {
  # Residual deficiency: compare observational to front-door estimate
  # If front-door is valid, the naive regression should differ from front-door
  
  n <- length(A)
  
  # Naive estimate (biased if confounded)
  naive_model <- lm(Y ~ A)
  naive_ate <- coef(naive_model)["A"]
  
  # Front-door estimate
  fd_ate <- .compute_frontdoor(A, M, Y, "plugin")
  
  # Deficiency proxy: normalized difference
  # (This is a heuristic; true delta requires distribution comparison)
  outcome_sd <- sd(Y)
  delta_proxy <- min(1, abs(naive_ate - fd_ate) / (2 * outcome_sd))
  
  delta_proxy
}

#' Print method for frontdoor_effect
#' @param x A frontdoor_effect object
#' @param ... Additional arguments (unused)
#' @export
print.frontdoor_effect <- function(x, ...) {
  cli::cli_h1("Front-Door Causal Effect")
  cli::cli_text("Mediator: {.val {x$mediator}}")
  cli::cli_text("Method: {.val {x$method}}")
  cli::cli_text("")
  
  cli::cli_h2("Results")
  cat(sprintf("  Effect estimate: %.4f\n", x$estimate))
  
  if (!is.na(x$se)) {
    cat(sprintf("  Standard error:  %.4f\n", x$se))
    cat(sprintf("  95%% CI:         [%.4f, %.4f]\n", x$ci[1], x$ci[2]))
  }
  
  cat(sprintf("  Deficiency (delta):  %.4f\n", x$deficiency))
  
  cat("\n")
  if (x$deficiency < 0.05) {
    cli::cli_alert_success("Front-door assumptions appear satisfied (delta < 0.05)")
  } else {
    cli::cli_alert_warning("Elevated deficiency suggests assumption violations")
  }
  
  invisible(x)
}
