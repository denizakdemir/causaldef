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
  title <- if (!is.null(x$metric) && identical(x$metric, "ps_tv")) {
    "Deficiency Proxy Estimates (PS-TV)"
  } else {
    "Deficiency Proxy Estimates"
  }
  cat("\n-- ", title, " ", paste(rep("-", max(0, 40 - nchar(title))), collapse = ""), "\n\n", sep = "")
  
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

  if (!is.null(x$metric) && identical(x$metric, "ps_tv")) {
    cat("Note: delta is a propensity-score TV proxy (overlap/balance diagnostic).\n")
  }

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
  cat("* delta bound (NC bound):", round(x$delta_bound, 4), "(kappa =", x$kappa, ")\n")
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
  cat("\n-- Policy Regret Bounds ", paste(rep("-", 49), collapse = ""), "\n\n", sep = "")
  cat("* Deficiency delta:", round(x$delta, 4), "\n")
  if (!is.null(x$delta_mode)) {
    mode_label <- if (identical(x$delta_mode, "upper")) "upper CI" else "point"
    cat("* Delta mode:", mode_label, "\n")
  }
  if (!is.null(x$delta_method) && nzchar(x$delta_method)) {
    cat("* Delta method:", x$delta_method, "\n")
  }
  cat("* Utility range: [", x$utility_range[1], ", ", x$utility_range[2], "]\n", sep = "")
  if (!is.null(x$transfer_penalty)) {
    cat("* Transfer penalty:", round(x$transfer_penalty, 4), "(additive regret upper bound)\n")
  } else {
    cat("* Transfer penalty:", round(x$safety_floor, 4), "(additive regret upper bound)\n")
  }
  if (!is.null(x$minimax_floor)) {
    cat("* Minimax floor:", round(x$minimax_floor, 4), "(worst-case lower bound)\n")
  }
  if (!is.null(x$complexity_penalty)) {
    cat("* Complexity penalty:", round(x$complexity_penalty, 4), "(VC / finite-sample term)\n")
  }

  if (!is.null(x$obs_regret)) {
    cat("\n* Observed regret:", round(x$obs_regret, 4), "\n")
    cat("* Interventional bound:", round(x$regret_bound, 4), "\n")
  }
  
  penalty <- if (!is.null(x$transfer_penalty)) x$transfer_penalty else x$safety_floor
  pct <- round(100 * penalty / x$M, 1)
  cat("\nInterpretation: Transfer penalty is", pct, "% of utility range given delta\n\n")
  
  invisible(x)
}

#' @export
print.partial_id_set <- function(x, ...) {
  cat("\n-- Partial Identification Set ", paste(rep("-", 46), collapse = ""), "\n\n", sep = "")
  cat("* Estimand:", x$estimand, "\n")
  cat("* Point estimate:", round(x$estimate, 6), "\n")
  cat("* Delta:", round(x$delta, 6), "\n")
  cat("* Half-width:", round(x$half_width, 6), "\n")
  cat("* Interval: [", round(x$lower, 6), ", ", round(x$upper, 6), "]\n\n", sep = "")
  invisible(x)
}

#' @export
print.overlap_diagnostic <- function(x, ...) {
  cat("\n-- Overlap Diagnostic ", paste(rep("-", 52), collapse = ""), "\n\n", sep = "")
  cat("* n:", x$n, "\n")
  cat("* trim:", x$trim, "\n")
  cat("* extreme ps count:", x$extreme_n, "\n")
  cat("* kept after trim:", x$kept_n, "\n")
  cat("* ESS (IPTW):", round(x$ess_iptw, 2), "\n")
  cat("\nPropensity quantiles:\n")
  print(round(x$ps_summary, 4))
  cat("\n")
  invisible(x)
}

#' @export
print.confounding_frontier <- function(x, ...) {
  cat("\n-- Confounding Frontier (Confounding Lower Bound) ", paste(rep("-", 26), collapse = ""), "\n\n", sep = "")
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

#' @export
print.data_audit_report <- function(x, ...) {
  cat("\n")
  cat("==============================================================================\n")
  cat("                         Data Integrity Report\n")
  cat("==============================================================================\n\n")
  
  cat("Treatment:", x$treatment, "| Outcome:", x$outcome, "\n")
  cat("Variables audited:", x$summary_stats$n_vars_audited, "\n")
  cat("Issues found:", x$summary_stats$n_issues, "\n\n")
  
  if (nrow(x$issues) > 0) {
    # Filter to show only significant issues (not "Safe")
    significant <- x$issues[!x$issues$issue_type %in% c("Safe", "Valid Negative Control"), ]
    
    if (nrow(significant) > 0) {
      cat("-- Issues Detected ", paste(rep("-", 50), collapse = ""), "\n\n", sep = "")
      
      # Print table
      display <- data.frame(
        Variable = significant$variable,
        Type = significant$issue_type,
        `p-value` = format.pval(significant$p_value, digits = 3),
        Recommendation = significant$recommendation,
        check.names = FALSE
      )
      print(display, row.names = FALSE, right = FALSE)
      cat("\n")
    }
  }
  
  cat("-- Recommendations ", paste(rep("-", 50), collapse = ""), "\n\n", sep = "")
  for (rec in x$recommendations) {
    cat("*", rec, "\n")
  }
  cat("\n")
  
  invisible(x)
}
