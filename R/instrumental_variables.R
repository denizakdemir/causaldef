# =============================================================================
# Instrumental Variable Estimation with Strength Diagnostics
# =============================================================================

#' Instrumental Variable Effect Estimation
#'
#' Estimates causal effects using instrumental variables (IV) together with
#' first-stage strength diagnostics.
#'
#' @param spec A causal_spec object with instrument specified
#' @param instrument Character: name of the instrumental variable (or use spec$instrument)
#' @param method Character: IV estimation method
#'   \itemize{
#'     \item "2sls": Two-stage least squares
#'     \item "wald": Wald estimator (ratio of reduced forms)
#'     \item "liml": Reserved for a future limited-information maximum likelihood
#'     implementation. The current release errors if requested.
#'   }
#' @param n_boot Integer: bootstrap replicates for CI
#' @param ci_level Numeric: confidence level
#' @param weak_iv_threshold Numeric: F-statistic threshold for weak IV (default 10)
#'
#' @return Object of class "iv_effect" containing:
#'   \itemize{
#'     \item estimate: The IV causal effect (LATE)
#'     \item se: Standard error
#'     \item ci: Confidence interval
#'     \item f_stat: First-stage F-statistic
#'     \item weak_iv: Logical, whether instrument is weak
#'     \item deficiency_proxy: Heuristic proxy derived from first-stage strength
#'     \item deficiency: Backward-compatible alias of `deficiency_proxy`
#'   }
#'
#' @details
#' Instrumental variables identify the Local Average Treatment Effect (LATE)
#' for compliers when:
#' \enumerate{
#'   \item Relevance: Z affects A (testable via F-stat)
#'   \item Exclusion: Z affects Y only through A (untestable)
#'   \item Independence: Z is independent of unmeasured confounders
#' }
#'
#' Weak instruments (small first-stage F-statistics) lead to unstable IV inference.
#' The current implementation reports a `deficiency_proxy` derived from instrument
#' strength; it should be interpreted as a screening diagnostic rather than an exact
#' deficiency.
#'
#' @examples
#' # Simulate IV setting
#' set.seed(42)
#' n <- 1000
#' U <- rnorm(n)  # Unmeasured confounder
#' Z <- rbinom(n, 1, 0.5)  # Instrument (randomized encouragement)
#' A <- 0.3 + 0.4 * Z + 0.3 * U + rnorm(n, sd = 0.3)  # Treatment (continuous)
#' A <- as.numeric(A > 0.5)  # Dichotomize
#' Y <- 1 + 2 * A + 0.8 * U + rnorm(n)  # Outcome
#'
#' df <- data.frame(Z = Z, A = A, Y = Y)
#' spec <- causal_spec(df, "A", "Y", instrument = "Z")
#' 
#' \donttest{
#' iv_result <- iv_effect(spec)
#' print(iv_result)
#' }
#'
#' @references
#' Angrist, J. D., Imbens, G. W., & Rubin, D. B. (1996). Identification of 
#' causal effects using instrumental variables. JASA.
#'
#' @seealso [causal_spec()], [estimate_deficiency()]
#' @export
iv_effect <- function(spec, instrument = NULL, 
                      method = c("2sls", "wald", "liml"),
                      n_boot = 200, ci_level = 0.95,
                      weak_iv_threshold = 10) {
  
  method <- match.arg(method)

  if (identical(method, "liml")) {
    cli::cli_abort(c(
      "Method {.val liml} is not available in this release.",
      "x" = "The previous approximation was not a valid LIML estimator.",
      "i" = "Use {.val 2sls} or {.val wald}."
    ))
  }
  
  # Validation
  if (!inherits(spec, "causal_spec")) {
    cli::cli_abort("spec must be a causal_spec object")
  }
  
  data <- spec$data
  
  # Get instrument
  if (is.null(instrument)) {
    if (!is.null(spec$instrument)) {
      instrument <- spec$instrument
    } else {
      cli::cli_abort("No instrument specified. Provide via {.arg instrument} or in {.fn causal_spec}.")
    }
  }
  
  if (!(instrument %in% names(data))) {
    cli::cli_abort("Instrument {.val {instrument}} not found in data")
  }
  
  A <- data[[spec$treatment]]
  Y <- data[[spec$outcome]]
  Z <- data[[instrument]]
  n <- nrow(data)
  
  # Get covariates for 2SLS
  W <- if (!is.null(spec$covariates)) {
    as.matrix(data[, spec$covariates, drop = FALSE])
  } else {
    NULL
  }
  
  # First stage: A ~ Z + W
  if (!is.null(W)) {
    first_stage <- lm(A ~ Z + W)
  } else {
    first_stage <- lm(A ~ Z)
  }
  
  # F-statistic for instrument strength
  first_stage_summary <- summary(first_stage)
  if ("Z" %in% rownames(first_stage_summary$coefficients)) {
    z_coef <- first_stage_summary$coefficients["Z", ]
    t_stat <- z_coef["t value"]
    f_stat <- t_stat^2
  } else {
    # Z might be renamed if factor
    z_idx <- grep("^Z", rownames(first_stage_summary$coefficients))
    if (length(z_idx) > 0) {
      t_stat <- first_stage_summary$coefficients[z_idx[1], "t value"]
      f_stat <- t_stat^2
    } else {
      f_stat <- NA
    }
  }
  
  weak_iv <- is.na(f_stat) || f_stat < weak_iv_threshold
  
  estimate <- .compute_iv_estimate(A = A, Y = Y, Z = Z, W = W, method = method)
  
  # Bootstrap for SE and CI
  if (n_boot > 0) {
    boot_fn <- function(i) {
      idx <- sample(n, replace = TRUE)
      boot_A <- A[idx]
      boot_Y <- Y[idx]
      boot_Z <- Z[idx]
      boot_W <- if (!is.null(W)) W[idx, , drop = FALSE] else NULL

      .compute_iv_estimate(
        A = boot_A,
        Y = boot_Y,
        Z = boot_Z,
        W = boot_W,
        method = method
      )
    }
    
    boot_estimates <- vapply(seq_len(n_boot), boot_fn, numeric(1))
    boot_estimates <- boot_estimates[is.finite(boot_estimates)]
    
    se <- sd(boot_estimates, na.rm = TRUE)
    alpha <- 1 - ci_level
    ci <- quantile(boot_estimates, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  } else {
    se <- NA_real_
    ci <- c(NA_real_, NA_real_)
  }
  
  # Heuristic proxy from instrument strength
  if (!is.na(f_stat) && f_stat > 0) {
    delta_iv <- min(1, sqrt(weak_iv_threshold / f_stat))
  } else {
    delta_iv <- 1.0
  }
  
  result <- structure(
    list(
      estimate = estimate,
      se = se,
      ci = ci,
      f_stat = f_stat,
      weak_iv = weak_iv,
      deficiency_proxy = delta_iv,
      deficiency = delta_iv,
      method = method,
      instrument = instrument,
      n = n,
      first_stage = first_stage
    ),
    class = "iv_effect"
  )
  
  # Report
  f_str <- if (!is.na(f_stat)) sprintf("%.1f", f_stat) else "NA"
  delta_str <- sprintf("%.3f", delta_iv)
  
  if (weak_iv) {
    cli::cli_alert_warning("Weak instrument (F = {f_str} < {weak_iv_threshold})")
    cli::cli_alert_warning("IV estimate: {round(estimate, 3)} (proxy = {delta_str})")
    cli::cli_alert_info("Consider finding a stronger instrument")
  } else {
    cli::cli_alert_success("Strong instrument (F = {f_str})")
    cli::cli_alert_success("IV estimate (LATE): {round(estimate, 3)} (proxy = {delta_str})")
  }
  
  result
}

#' @keywords internal
.compute_iv_estimate <- function(A, Y, Z, W = NULL, method = c("2sls", "wald")) {
  method <- match.arg(method)

  if (method == "wald") {
    cov_yz <- cov(Y, Z)
    cov_az <- cov(A, Z)

    if (abs(cov_az) < 1e-10) {
      cli::cli_abort("Instrument has no variation with treatment (Cov(A,Z) ~ 0)")
    }

    return(unname(cov_yz / cov_az))
  }

  if (!is.null(W)) {
    first_stage <- lm(A ~ Z + W)
    A_hat <- fitted(first_stage)
    second_stage <- lm(Y ~ A_hat + W)
  } else {
    first_stage <- lm(A ~ Z)
    A_hat <- fitted(first_stage)
    second_stage <- lm(Y ~ A_hat)
  }

  unname(coef(second_stage)["A_hat"])
}

#' Print method for iv_effect
#' @param x An iv_effect object
#' @param ... Additional arguments (unused)
#'
#' @return Invisibly returns `x`, unchanged. The printed output reports the
#'   instrument and estimation method, first-stage F-statistic, weak-IV flag,
#'   IV effect estimate, optional bootstrap uncertainty, and the heuristic
#'   strength proxy stored in `x$deficiency_proxy`.
#' @export
print.iv_effect <- function(x, ...) {
  cli::cli_h1("Instrumental Variable Effect")
  cli::cli_text("Instrument: {.val {x$instrument}}")
  cli::cli_text("Method: {.val {x$method}}")
  cli::cli_text("")
  
  cli::cli_h2("First Stage")
  f_str <- if (!is.na(x$f_stat)) sprintf("%.2f", x$f_stat) else "NA"
  cat(sprintf("  F-statistic: %s\n", f_str))
  cat(sprintf("  Weak IV: %s\n", if (x$weak_iv) "YES (F < 10)" else "No"))
  
  cli::cli_h2("IV Estimate (LATE)")
  cat(sprintf("  Effect:      %.4f\n", x$estimate))
  
  if (!is.na(x$se)) {
    cat(sprintf("  Std. Error:  %.4f\n", x$se))
    cat(sprintf("  95%% CI:     [%.4f, %.4f]\n", x$ci[1], x$ci[2]))
  }
  
  cat(sprintf("  Strength proxy: %.4f\n", x$deficiency_proxy))
  
  cat("\n")
  if (x$weak_iv) {
    cli::cli_alert_danger("WARNING: Weak instrument may bias estimates toward OLS")
    cli::cli_alert_info("Consider: (1) stronger instrument, (2) more observations, (3) alternative robust IV methods outside the current package release")
  } else if (x$deficiency_proxy < 0.2) {
    cli::cli_alert_success("First-stage strength looks good under the current proxy.")
  } else {
    cli::cli_alert_warning("Moderate instrument weakness under the current proxy.")
  }
  
  invisible(x)
}

#' Summary method for iv_effect
#' @param object An iv_effect object
#' @param ... Additional arguments (unused)
#'
#' @return Invisibly returns `object`, unchanged. The printed summary includes
#'   the instrument, estimation method, sample size, the coefficient table from
#'   the first-stage regression stored in `object$first_stage`, the IV effect
#'   estimate with optional bootstrap uncertainty, and the weak-instrument
#'   diagnostics.
#' @export
summary.iv_effect <- function(object, ...) {
  cat("\n=== Instrumental Variable Analysis Summary ===\n\n")
  
  cat("Instrument:", object$instrument, "\n")
  cat("Method:", object$method, "\n")
  cat("Sample size:", object$n, "\n\n")
  
  cat("First Stage:\n")
  print(summary(object$first_stage)$coefficients)
  
  cat("\nIV Estimate:\n")
  cat("  LATE =", round(object$estimate, 4), "\n")
  if (!is.na(object$se)) {
    cat("  SE   =", round(object$se, 4), "\n")
    cat("  95% CI: [", round(object$ci[1], 4), ",", round(object$ci[2], 4), "]\n")
  }
  
  cat("\nDiagnostics:\n")
  cat("  F-statistic:", round(object$f_stat, 2), "\n")
  cat("  Weak IV:", if (object$weak_iv) "Yes" else "No", "\n")
  cat("  Strength proxy:", round(object$deficiency_proxy, 4), "\n")
  
  invisible(object)
}

#' Test Instrument Validity
#'
#' Performs diagnostic tests for instrumental variable assumptions
#'
#' @param spec A causal_spec with instrument
#' @param instrument Instrument name
#' @param overid_test Logical: perform overidentification test if multiple instruments
#'
#' @return A named list of IV screening diagnostics with components such as:
#'   \itemize{
#'     \item `relevance`: first-stage F-statistic, its p-value, and a logical
#'     pass/fail flag relative to the weak-IV threshold of 10
#'     \item `balance`: when covariates are present, per-covariate balance
#'     p-values and an overall pass/fail flag for instrument-covariate balance
#'     \item `reduced_form`: coefficient, standard error, p-value, and
#'     significance flag for the reduced-form association of the instrument with
#'     the outcome
#'   }
#'   These outputs are diagnostics for testable parts of the IV design. They do
#'   not by themselves prove the exclusion restriction or full instrument
#'   validity.
#' @export
test_instrument <- function(spec, instrument = NULL, overid_test = TRUE) {
  
  data <- spec$data
  
  if (is.null(instrument)) {
    instrument <- spec$instrument
  }
  
  if (is.null(instrument)) {
    cli::cli_abort("No instrument specified")
  }
  
  A <- data[[spec$treatment]]
  Y <- data[[spec$outcome]]
  Z <- data[[instrument]]
  n <- nrow(data)
  
  results <- list()
  
  # Test 1: Relevance (F-test)
  first_stage <- lm(A ~ Z)
  f_test <- summary(first_stage)$fstatistic
  if (!is.null(f_test)) {
    f_stat <- f_test[1]
    f_pval <- pf(f_stat, f_test[2], f_test[3], lower.tail = FALSE)
    results$relevance <- list(
      f_stat = f_stat,
      p_value = f_pval,
      pass = f_stat > 10
    )
  }
  
  # Test 2: Balance check (if covariates available)
  if (!is.null(spec$covariates)) {
    balance_pvals <- vapply(spec$covariates, function(cov) {
      x <- data[[cov]]
      if (is.numeric(x)) {
        t.test(x ~ Z)$p.value
      } else {
        chisq.test(table(x, Z))$p.value
      }
    }, numeric(1))
    
    results$balance <- list(
      p_values = balance_pvals,
      pass = all(balance_pvals > 0.05)
    )
  }
  
  # Test 3: Reduced form (Z -> Y)
  reduced_form <- lm(Y ~ Z)
  rf_summary <- summary(reduced_form)
  z_effect <- rf_summary$coefficients["Z", ]
  
  results$reduced_form <- list(
    estimate = z_effect["Estimate"],
    se = z_effect["Std. Error"],
    p_value = z_effect["Pr(>|t|)"],
    significant = z_effect["Pr(>|t|)"] < 0.05
  )
  
  # Report
  cli::cli_h2("Instrument Validity Tests")
  
  cli::cli_text("Relevance (F-test):")
  if (!is.null(results$relevance)) {
    if (results$relevance$pass) {
      cli::cli_alert_success("  F = {round(results$relevance$f_stat, 2)} > 10 [OK]")
    } else {
      cli::cli_alert_danger("  F = {round(results$relevance$f_stat, 2)} < 10 [X]")
    }
  }
  
  if (!is.null(results$balance)) {
    cli::cli_text("Covariate Balance:")
    if (results$balance$pass) {
      cli::cli_alert_success("  All covariates balanced (p > 0.05) [OK]")
    } else {
      cli::cli_alert_warning("  Some covariates imbalanced")
    }
  }
  
  cli::cli_text("Reduced Form (Z -> Y):")
  if (results$reduced_form$significant) {
    cli::cli_alert_success("  Significant (p = {round(results$reduced_form$p_value, 4)}) [OK]")
  } else {
    cli::cli_alert_warning("  Not significant (p = {round(results$reduced_form$p_value, 4)})")
  }
  
  results
}
