# =============================================================================
# Print Methods for causaldef Classes
# =============================================================================

#' @export
print.causal_spec <- function(x, ...) {
  cat("\n-- Causal Specification ", paste(rep("-", 50), collapse = ""), "\n\n", sep = "")
  cat("* Treatment:", x$treatment, "(", x$treatment_type, ")\n")
  cat("* Outcome:", x$outcome, "(", x$outcome_type, ")\n")
  cat("* Covariates:", paste(x$covariates, collapse = ", "), "\n")
  cat("* Sample size:", x$n, "\n")
  cat("* Estimand:", x$estimand, "\n")
  
  if (!is.null(x$negative_control)) {
    cat("* Negative control:", x$negative_control, "\n")
  }
  if (!is.null(x$instrument)) {
    cat("* Instrument:", x$instrument, "\n")
  }
  cat("\n")
  invisible(x)
}

#' @export
print.deficiency <- function(x, ...) {
  cat("\n-- Le Cam Deficiency Estimates ", paste(rep("-", 40), collapse = ""), "\n\n", sep = "")
  
  # Create results table
  results_df <- data.frame(
    Method = names(x$estimates),
    Delta = round(x$estimates, 4),
    SE = if (all(is.na(x$se))) "-" else round(x$se, 4),
    stringsAsFactors = FALSE
  )
  
  if (!all(is.na(x$ci))) {
    results_df$CI <- paste0("[", round(x$ci[, 1], 4), ", ", round(x$ci[, 2], 4), "]")
  }
  
  results_df$Quality <- sapply(x$estimates, function(d) {
    if (d < 0.05) return("Excellent (Green)")
    if (d < 0.10) return("Caution (Yellow)")
    return("Insufficient (Red)")
  })

  print(results_df, row.names = FALSE)
  
  # Best method
  best_idx <- which.min(x$estimates)
  cat("\nBest method:", names(x$estimates)[best_idx], 
      "(delta =", round(x$estimates[best_idx], 4), ")\n")
  
  invisible(x)
}

#' @export
print.nc_diagnostic <- function(x, ...) {
  cat("\n-- Negative Control Diagnostic ", paste(rep("-", 40), collapse = ""), "\n\n", sep = "")
  cat("* delta_NC (observable):", round(x$delta_nc, 4), "\n")
  cat("* delta bound (Theorem 5.2):", round(x$delta_bound, 4), "(kappa =", x$kappa, ")\n")
  cat("* p-value:", format.pval(x$p_value), "\n\n")
  
  if (x$falsified) {
    cat("RESULT: REJECTED. Evidence of residual confounding.\n")
  } else {
    cat("RESULT: NOT REJECTED. We failed to catch an error, but that doesn't mean an error isn't there.\n")
    cat("NOTE: Your effect estimate must exceed the Noise Floor (delta_bound) to be meaningful.\n")
  }
  cat("\n")
  invisible(x)
}

#' @export
print.policy_bound <- function(x, ...) {
  cat("\n-- Policy Regret Bound (Theorem 3.2) ", paste(rep("-", 35), collapse = ""), "\n\n", sep = "")
  cat("* Deficiency delta:", round(x$delta, 4), "\n")
  cat("* Utility range: [", x$utility_range[1], ", ", x$utility_range[2], "]\n", sep = "")
  cat("* Safety floor:", round(x$safety_floor, 4), "(minimum regret given delta)\n")
  
  if (!is.null(x$obs_regret)) {
    cat("\n* Observed regret:", round(x$obs_regret, 4), "\n")
    cat("* Interventional bound:", round(x$regret_bound, 4), "\n")
  }
  
  pct <- round(100 * x$safety_floor / x$M, 1)
  cat("\nInterpretation: Worst-case regret is", pct, "% of utility range due to confounding\n\n")
  
  invisible(x)
}

#' @export
print.confounding_frontier <- function(x, ...) {
  cat("\n-- Confounding Frontier (Theorem 4.1) ", paste(rep("-", 35), collapse = ""), "\n\n", sep = "")
  cat("* Grid size:", x$params$grid_size, "x", x$params$grid_size, "\n")
  cat("* alpha range: [", x$params$alpha_range[1], ", ", x$params$alpha_range[2], "]\n", sep = "")
  cat("* gamma range: [", x$params$gamma_range[1], ", ", x$params$gamma_range[2], "]\n", sep = "")
  cat("* Model:", x$model, "\n\n")
  
  # Summary statistics
  cat("-- Deficiency Summary --\n\n")
  cat("* Min delta:", round(min(x$grid$delta), 4), "\n")
  cat("* Max delta:", round(max(x$grid$delta), 4), "\n")
  cat("* Mean delta:", round(mean(x$grid$delta), 4), "\n\n")
  
  # Identification zone
  n_identified <- sum(x$grid$delta < 0.01)
  pct_identified <- round(100 * n_identified / nrow(x$grid), 1)
  cat(pct_identified, "% of grid has near-zero deficiency (delta < 0.01)\n\n")
  
  invisible(x)
}
