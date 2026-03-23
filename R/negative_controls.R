# =============================================================================
# nc_diagnostic() - Negative Control Diagnostic
# =============================================================================

#' Negative Control Diagnostic
#'
#' Screens for residual association with a negative control outcome and computes
#' the corresponding \eqn{\kappa}-based sensitivity bound.
#'
#' The diagnostic uses the bound: \eqn{\delta(\hat{K}) \le \kappa \cdot \delta_{NC}(\hat{K})}, where \eqn{\delta_{NC}}
#' is the negative control deficiency.
#'
#' @param spec A causal_spec with negative_control specified
#' @param method Character: adjustment method to test
#' @param kappa Numeric: alignment constant in the theoretical bound. If \code{NULL},
#'   defaults to 1 (conservative).
#' @param kappa_range Numeric vector: range of kappa values for sensitivity analysis.
#'   If provided, returns bounds for each kappa and adds class "nc_diagnostic_sensitivity".
#' @param alpha Numeric: significance level for the screening test
#' @param n_boot Integer: number of permutation draws used for the screening test
#'
#' @return Object of class "nc_diagnostic" containing:
#'   \itemize{
#'     \item delta_nc: Observable residual-association proxy based on a weighted
#'     treatment/negative-control association
#'     \item delta_bound: Upper bound on true deficiency (kappa * delta_nc)
#'     \item falsified: Logical indicating if assumptions are falsified
#'     \item p_value: Permutation p-value for the observed residual association
#'     \item kappa: Alignment constant used
#'     \item screening: List with the observed association statistic, permutation
#'     p-value, and null standard deviation
#'     \item sensitivity: (if kappa_range provided) data.frame with bounds for each kappa
#'   }
#'
#' @details
#' A negative control outcome Y' is a variable that:
#' \enumerate{
#'   \item Shares confounders with Y (affected by U)
#'   \item Is NOT affected by treatment A
#' }
#' 
#' If adjustment correctly removes confounding, the residual association
#' between A and Y' should be zero. Non-zero association indicates
#' confounding remains.
#'
#' The current implementation separates two tasks:
#' \enumerate{
#'   \item a \strong{screening test} based on a weighted correlation statistic and a
#'   permutation null distribution, and
#'   \item a \strong{sensitivity bound} of the form \eqn{\delta \le \kappa \delta_{NC}}
#'   based on the observed residual-association proxy `delta_nc`.
#' }
#'
#' This keeps the observable screening step distinct from the theorem's
#' user-supplied alignment parameter `kappa`.
#' 
#' When kappa is unknown, use kappa_range to perform sensitivity analysis
#' across plausible values. The resulting plot shows how the deficiency
#' bound varies with the sensitivity parameter.
#'
#' @examples
#' # Create data with negative control
#' n <- 200
#' U <- rnorm(n)
#' W <- U + rnorm(n, sd = 0.5)
#' A <- rbinom(n, 1, plogis(0.5 * W))  # Binary treatment
#' Y <- 1 + 2 * A + U + rnorm(n)
#' Y_nc <- U + rnorm(n)  # Shares U but no effect from A
#' df <- data.frame(W = W, A = A, Y = Y, Y_nc = Y_nc)
#'
#' spec <- causal_spec(df, "A", "Y", "W", negative_control = "Y_nc")
#' \donttest{
#' # Standard diagnostic
#' nc_result <- nc_diagnostic(spec, method = "iptw", n_boot = 50)
#' print(nc_result)
#' 
#' # Sensitivity analysis
#' nc_sens <- nc_diagnostic(spec, method = "iptw", 
#'                          kappa_range = seq(0.5, 2.0, by = 0.25),
#'                          n_boot = 50)
#' plot(nc_sens)
#' }
#'
#' @references
#' Akdemir, D. (2026). Constraints on Causal Inference as Experiment Comparison.
#' DOI: 10.5281/zenodo.18367347. See `thm:nc_bound` (Negative Control Sensitivity Bound).
#' 
#' Lipsitch, M., Tchetgen, E., & Cohen, T. (2010). Negative controls: A tool
#' for detecting confounding and bias. Epidemiology, 21(3), 383-388.
#'
#' @seealso [causal_spec()], [estimate_deficiency()]
#' @export
nc_diagnostic <- function(spec, method = "iptw", kappa = NULL,
                          kappa_range = NULL,
                          alpha = 0.05, n_boot = 200) {
  
  validate_causal_spec(spec)
  
  # Check negative control is specified
  if (is.null(spec$negative_control)) {
    .msg_error("No negative control specified. Add negative_control to causal_spec()")
  }
  
  checkmate::assert_number(alpha, lower = 0.001, upper = 0.5)
  checkmate::assert_integerish(n_boot, lower = 10)
  
  # Use kappa if provided, otherwise default to 1
  if (is.null(kappa)) {
    kappa <- 1.0  # Conservative default
    .msg_info("Using kappa = 1 (conservative). Consider domain-specific estimation or sensitivity analysis via kappa_range.")
  }
  
  data <- spec$data
  
  # Build kernel using the specified method
  tr_var <- data[[spec$treatment]]
  treatment_value <- if (is.factor(tr_var)) {
    levels(tr_var)[2]
  } else if (is.character(tr_var)) {
    lvls <- sort(unique(tr_var))
    if (length(lvls) != 2) {
      cli::cli_abort("Treatment variable {.val {spec$treatment}} must have exactly 2 levels for nc_diagnostic().")
    }
    lvls[2]
  } else {
    1
  }
  kernel <- .build_kernel(spec, method, treatment_value = treatment_value)
  
  # Screening test for residual A-Y' association
  nc_result <- .screen_negative_control(spec, kernel, n_boot)
  
  # Apply negative-control bound: delta <= kappa * delta_nc
  delta_bound <- kappa * nc_result$delta_nc
  
  # Falsification decision
  falsified <- nc_result$p_value < alpha
  
  result <- new_nc_diagnostic(
    delta_nc = nc_result$delta_nc,
    delta_bound = delta_bound,
    falsified = falsified,
    p_value = nc_result$p_value,
    kappa = kappa,
    kernel = kernel,
    screening = list(
      statistic = nc_result$statistic,
      p_value = nc_result$p_value,
      null_sd = nc_result$null_sd,
      method = "weighted_permutation_correlation"
    ),
    test_statistic = nc_result$statistic,
    se = nc_result$null_sd,
    p_value_method = "permutation"
  )
  
  # Sensitivity analysis if kappa_range provided
  if (!is.null(kappa_range)) {
    checkmate::assert_numeric(kappa_range, lower = 0.01, min.len = 2)
    sensitivity <- data.frame(
      kappa = kappa_range,
      delta_bound = kappa_range * nc_result$delta_nc
    )
    result$sensitivity <- sensitivity
    result$kappa_range <- kappa_range
    class(result) <- c("nc_diagnostic_sensitivity", class(result))
    .msg_info(sprintf("Sensitivity analysis: kappa in [%.2f, %.2f], delta_bound in [%.4f, %.4f]",
                      min(kappa_range), max(kappa_range),
                      min(sensitivity$delta_bound), max(sensitivity$delta_bound)))
  }
  
  if (falsified) {
    .msg_danger(paste("Causal assumptions FALSIFIED (p =", format.pval(nc_result$p_value), ")"))
  } else {
    .msg_success(paste("No evidence against causal assumptions (p =", format.pval(nc_result$p_value), ")"))
  }
  
  result
}

# =============================================================================
# Internal: Negative Control Screening Test
# =============================================================================

#' @keywords internal
.screen_negative_control <- function(spec, kernel, n_boot) {
  
  data <- spec$data
  A <- data[[spec$treatment]]
  if (is.factor(A) || is.character(A)) {
    A <- as.numeric(as.factor(A)) - 1
  } else if (is.logical(A)) {
    A <- as.integer(A)
  } else if (!is.numeric(A)) {
    A <- as.numeric(A)
  }
  Y_nc <- data[[spec$negative_control]]
  if (!is.numeric(Y_nc)) Y_nc <- as.numeric(as.factor(Y_nc))
  n <- nrow(data)
  
  # Get weights from kernel
  weights <- if (!is.null(kernel$weights)) {
    kernel$weights
  } else {
    rep(1, n)
  }
  
  # Test statistic: weighted correlation between A and Y_nc
  compute_nc_association <- function(A, Y_nc, weights) {
    w_sum <- sum(weights)
    if (!is.finite(w_sum) || w_sum <= 0) return(0)
    A_wtd <- sum(weights * A) / w_sum
    Y_nc_wtd <- sum(weights * Y_nc) / w_sum
    
    cov_wtd <- sum(weights * (A - A_wtd) * (Y_nc - Y_nc_wtd)) / w_sum
    var_A <- sum(weights * (A - A_wtd)^2) / w_sum
    var_Y <- sum(weights * (Y_nc - Y_nc_wtd)^2) / w_sum
    
    if (var_A == 0 || var_Y == 0) return(0)
    
    cov_wtd / sqrt(var_A * var_Y)
  }
  
  # Point estimate
  stat_obs <- compute_nc_association(A, Y_nc, weights)
  
  # Permutation null distribution for residual association
  perm_stats <- replicate(n_boot, {
    perm_idx <- sample.int(n)
    compute_nc_association(A, Y_nc[perm_idx], weights)
  })
  
  null_sd <- sd(perm_stats)
  p_value <- (1 + sum(abs(perm_stats) >= abs(stat_obs))) / (n_boot + 1)
  
  # Convert correlation to deficiency-like measure [0, 1]
  delta_nc <- min(1, abs(stat_obs))
  
  list(
    delta_nc = delta_nc,
    statistic = stat_obs,
    null_sd = null_sd,
    p_value = p_value
  )
}
