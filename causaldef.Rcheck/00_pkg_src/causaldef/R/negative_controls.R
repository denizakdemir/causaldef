# =============================================================================
# nc_diagnostic() - Negative Control Diagnostic
# =============================================================================

#' Negative Control Diagnostic
#'
#' Tests causal assumptions using negative control outcomes. If the
#' adjustment strategy correctly removes confounding, it should also
#' remove the spurious association with a negative control outcome
#' (one that shares confounders but is unaffected by treatment).
#'
#' The diagnostic uses the bound: \eqn{\delta(\hat{K}) \le \kappa \cdot \delta_{NC}(\hat{K})}, where \eqn{\delta_{NC}}
#' is the negative control deficiency.
#'
#' @param spec A causal_spec with negative_control specified
#' @param method Character: adjustment method to test
#' @param kappa Numeric: alignment constant (NULL = use 1.0 as conservative default)
#' @param alpha Numeric: significance level for falsification test
#' @param n_boot Integer: bootstrap iterations
#'
#' @return Object of class "nc_diagnostic" containing:
#'   \itemize{
#'     \item delta_nc: Negative control deficiency (observable)
#'     \item delta_bound: Upper bound on true deficiency (kappa * delta_nc)
#'     \item falsified: Logical indicating if assumptions are falsified
#'     \item p_value: P-value for falsification test
#'     \item kappa: Alignment constant used
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
#' nc_result <- nc_diagnostic(spec, method = "iptw", n_boot = 50)
#' print(nc_result)
#' }
#'
#' @references
#' Akdemir, D. (2026). Constraints on Causal Inference as Experiment Comparison.
#' DOI: 10.5281/zenodo.18367347. See Theorem 5.2 for theoretical details.
#' 
#' Lipsitch, M., Tchetgen, E., & Cohen, T. (2010). Negative controls: A tool
#' for detecting confounding and bias. Epidemiology, 21(3), 383-388.
#'
#' @seealso [causal_spec()], [estimate_deficiency()]
#' @export
nc_diagnostic <- function(spec, method = "iptw", kappa = NULL,
                          alpha = 0.05, n_boot = 200) {
  
  validate_causal_spec(spec)
  
  # Check negative control is specified
  if (is.null(spec$negative_control)) {
    .msg_error("No negative control specified. Add negative_control to causal_spec()")
  }
  
  checkmate::assert_number(alpha, lower = 0.001, upper = 0.5)
  checkmate::assert_integerish(n_boot, lower = 10)
  
  # Estimate kappa if not provided
  if (is.null(kappa)) {
    kappa <- 1.0  # Conservative default
    .msg_info("Using kappa = 1 (conservative). Consider domain-specific estimation.")
  }
  
  data <- spec$data
  
  # Build kernel using the specified method
  kernel <- .build_kernel(spec, method, treatment_value = 1)
  
  # Compute NC deficiency (test for residual A-Y' association)
  nc_result <- .test_negative_control(spec, kernel, n_boot)
  
  # Apply Theorem 5.2 bound: delta <= kappa * delta_nc
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
    test_statistic = nc_result$statistic,
    se = nc_result$se
  )
  
  if (falsified) {
    .msg_danger(paste("Causal assumptions FALSIFIED (p =", format.pval(nc_result$p_value), ")"))
  } else {
    .msg_success(paste("No evidence against causal assumptions (p =", format.pval(nc_result$p_value), ")"))
  }
  
  result
}

# =============================================================================
# Internal: Negative Control Test
# =============================================================================

#' @keywords internal
.test_negative_control <- function(spec, kernel, n_boot) {
  
  data <- spec$data
  A <- data[[spec$treatment]]
  if (!is.numeric(A)) A <- as.numeric(A)
  Y_nc <- data[[spec$negative_control]]
  if (!is.numeric(Y_nc)) Y_nc <- as.numeric(Y_nc)
  n <- nrow(data)
  
  # Get weights from kernel
  weights <- if (!is.null(kernel$weights)) {
    kernel$weights
  } else {
    rep(1, n)
  }
  
  # Test statistic: weighted correlation between A and Y_nc
  # Under null (correct adjustment), this should be ~0
  compute_nc_association <- function(A, Y_nc, weights) {
    # Weighted correlation as test statistic
    w_sum <- sum(weights)
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
  
  # Bootstrap for inference
  boot_stats <- replicate(n_boot, {
    boot_idx <- sample(n, replace = TRUE)
    compute_nc_association(A[boot_idx], Y_nc[boot_idx], weights[boot_idx])
  })
  
  se <- sd(boot_stats)
  
  # Two-sided p-value (null: association = 0)
  z_stat <- abs(stat_obs) / se
  p_value <- 2 * (1 - pnorm(z_stat))
  
  # Convert correlation to deficiency-like measure [0, 1]
  # Using Fisher's z or just absolute correlation
  delta_nc <- min(1, abs(stat_obs))
  
  list(
    delta_nc = delta_nc,
    statistic = stat_obs,
    se = se,
    p_value = p_value
  )
}
