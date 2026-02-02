# =============================================================================
# Transport Deficiency: Measuring Distribution Shift
# =============================================================================

#' Transport Deficiency Between Source and Target Populations
#'
#' Computes the deficiency for transporting causal effects from a source 
#' population to a target population. This quantifies how much information 
#' is lost when extrapolating treatment effects across populations.
#'
#' @param source_spec A causal_spec object for the source (training) population
#' @param target_data Data frame: the target population (may have fewer variables)
#' @param transport_vars Character vector: variables that may shift between populations
#' @param method Character: method for transport estimation
#'   \itemize{
#'     \item "iptw": Inverse probability of trial/source weighting
#'     \item "calibration": Calibration weighting
#'     \item "sace": Sample average causal effect
#'   }
#' @param n_boot Integer: bootstrap replicates (default 100)
#' @param ci_level Numeric: confidence level (default 0.95)
#'
#' @return Object of class "transport_deficiency" containing:
#'   \itemize{
#'     \item delta_transport: The transport deficiency
#'     \item covariate_shift: Per-variable shift diagnostics
#'     \item weights: Transport weights for source observations
#'     \item ate_source: ATE estimate in source population
#'     \item ate_target: Transported ATE estimate for target
#'   }
#'
#' @details
#' Transport deficiency measures how much causal information is lost when 
#' moving from source S to target T:
#'
#' \deqn{\delta_{transport} = \delta(\mathcal{E}_S, \mathcal{E}_T)}
#'
#' Key assumptions for valid transport:
#' \enumerate{
#'   \item Treatment effects are identifiable in source
#'   \item Effect homogeneity or sufficient transport variables
#'   \item Positivity: overlap between source and target covariate distributions
#' }
#'
#' @examples
#' # Source population (RCT)
#' set.seed(42)
#' n_source <- 500
#' age_s <- rnorm(n_source, 50, 10)
#' A_s <- rbinom(n_source, 1, 0.5)  # Randomized
#' Y_s <- 10 + 2 * A_s - 0.1 * age_s + rnorm(n_source)
#' source_df <- data.frame(age = age_s, A = A_s, Y = Y_s, S = 1)
#'
#' # Target population (different age distribution)
#' n_target <- 300
#' age_t <- rnorm(n_target, 65, 8)  # Older population
#' target_df <- data.frame(age = age_t, S = 0)
#'
#' source_spec <- causal_spec(source_df, "A", "Y", "age")
#' \dontrun{
#' transport <- transport_deficiency(
#'   source_spec,
#'   target_data = target_df,
#'   transport_vars = "age"
#' )
#' print(transport)
#' }
#'
#' @seealso [estimate_deficiency()], [policy_regret_bound()]
#' @export
transport_deficiency <- function(source_spec, target_data, 
                                  transport_vars = NULL,
                                  method = c("iptw", "calibration", "sace"),
                                  n_boot = 100, ci_level = 0.95) {
  
  method <- match.arg(method)
  
  # Validation
  if (!inherits(source_spec, "causal_spec")) {
    cli::cli_abort("source_spec must be a causal_spec object")
  }
  if (!is.data.frame(target_data)) {
    cli::cli_abort("target_data must be a data frame")
  }
  
  source_data <- source_spec$data
  n_source <- nrow(source_data)
  n_target <- nrow(target_data)
  
  # Infer transport variables if not specified
  if (is.null(transport_vars)) {
    transport_vars <- intersect(names(source_data), names(target_data))
    transport_vars <- setdiff(transport_vars, c(source_spec$treatment, source_spec$outcome))
    cli::cli_alert_info("Using transport variables: {.val {transport_vars}}")
  }
  
  # Check all transport vars exist in both datasets
  missing_source <- setdiff(transport_vars, names(source_data))
  missing_target <- setdiff(transport_vars, names(target_data))
  
  if (length(missing_source) > 0) {
    cli::cli_abort("Transport variables missing from source: {.val {missing_source}}")
  }
  if (length(missing_target) > 0) {
    cli::cli_abort("Transport variables missing from target: {.val {missing_target}}")
  }
  
  # Compute covariate shift diagnostics
  cov_shift <- .compute_covariate_shift(
    source_data[, transport_vars, drop = FALSE],
    target_data[, transport_vars, drop = FALSE]
  )
  
  # Compute transport weights
  weights_result <- .compute_transport_weights(
    source_data, target_data, transport_vars, method
  )
  
  # Estimate ATEs
  ate_source <- .estimate_source_ate(source_spec)
  ate_target <- .estimate_transported_ate(source_spec, weights_result$weights)
  
  # Compute transport deficiency
  delta_transport <- .compute_transport_delta(
    cov_shift, 
    weights_result$extreme_weights,
    n_source, n_target
  )
  
  # Bootstrap for uncertainty
  if (n_boot > 0) {
    boot_deltas <- replicate(n_boot, {
      # Resample source
      idx_s <- sample(n_source, replace = TRUE)
      boot_source <- source_data[idx_s, , drop = FALSE]
      
      # Resample target
      idx_t <- sample(n_target, replace = TRUE)
      boot_target <- target_data[idx_t, , drop = FALSE]
      
      boot_shift <- .compute_covariate_shift(
        boot_source[, transport_vars, drop = FALSE],
        boot_target[, transport_vars, drop = FALSE]
      )
      
      mean(boot_shift$shift_metric, na.rm = TRUE)
    })
    
    se <- sd(boot_deltas, na.rm = TRUE)
    alpha <- 1 - ci_level
    ci <- quantile(boot_deltas, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  } else {
    se <- NA_real_
    ci <- c(NA_real_, NA_real_)
  }
  
  result <- structure(
    list(
      delta_transport = delta_transport,
      se = se,
      ci = ci,
      covariate_shift = cov_shift,
      weights = weights_result$weights,
      effective_sample_size = weights_result$ess,
      ate_source = ate_source,
      ate_target = ate_target,
      method = method,
      transport_vars = transport_vars,
      n_source = n_source,
      n_target = n_target
    ),
    class = "transport_deficiency"
  )
  
  # Report
  delta_str <- sprintf("%.3f", delta_transport)
  ess_ratio <- sprintf("%.1f%%", 100 * weights_result$ess / n_source)
  
  if (delta_transport < 0.1) {
    cli::cli_alert_success("Transport deficiency: {delta_str} (ESS: {ess_ratio})")
    cli::cli_alert_info("Effect transport appears reliable")
  } else if (delta_transport < 0.25) {
    cli::cli_alert_warning("Transport deficiency: {delta_str} (ESS: {ess_ratio})")
    cli::cli_alert_info("Moderate distribution shift; interpret with caution")
  } else {
    cli::cli_alert_danger("Transport deficiency: {delta_str} (ESS: {ess_ratio})")
    cli::cli_alert_info("Severe distribution shift; transport may be unreliable")
  }
  
  result
}

#' @keywords internal
.compute_covariate_shift <- function(source_covs, target_covs) {
  var_names <- names(source_covs)
  
  shift_metrics <- vapply(var_names, function(v) {
    x_s <- source_covs[[v]]
    x_t <- target_covs[[v]]
    
    if (is.numeric(x_s) && is.numeric(x_t)) {
      # Standardized mean difference
      pooled_sd <- sqrt((var(x_s) + var(x_t)) / 2)
      if (pooled_sd > 0) {
        smd <- abs(mean(x_t) - mean(x_s)) / pooled_sd
      } else {
        smd <- 0
      }
      smd
    } else {
      # For categorical: total variation
      tab_s <- prop.table(table(x_s))
      tab_t <- prop.table(table(x_t))
      all_levels <- union(names(tab_s), names(tab_t))
      
      p_s <- setNames(rep(0, length(all_levels)), all_levels)
      p_t <- setNames(rep(0, length(all_levels)), all_levels)
      p_s[names(tab_s)] <- tab_s
      p_t[names(tab_t)] <- tab_t
      
      0.5 * sum(abs(p_s - p_t))
    }
  }, numeric(1))
  
  data.frame(
    variable = var_names,
    shift_metric = shift_metrics,
    severity = cut(shift_metrics, 
                   breaks = c(-Inf, 0.1, 0.25, 0.5, Inf),
                   labels = c("low", "moderate", "high", "severe"))
  )
}

#' @keywords internal
.compute_transport_weights <- function(source_data, target_data, 
                                        transport_vars, method) {
  n_source <- nrow(source_data)
  n_target <- nrow(target_data)
  
  # Combine data with source indicator
  source_subset <- source_data[, transport_vars, drop = FALSE]
  target_subset <- target_data[, transport_vars, drop = FALSE]
  
  source_subset$S <- 1
  target_subset$S <- 0
  
  combined <- rbind(source_subset, target_subset)
  
  if (method == "iptw") {
    # Probability of being in target vs source
    ps_formula <- as.formula(paste("S ~", paste(transport_vars, collapse = " + ")))
    
    ps_model <- tryCatch({
      glm(ps_formula, data = combined, family = binomial())
    }, error = function(e) NULL)
    
    if (!is.null(ps_model)) {
      ps_source <- predict(ps_model, newdata = source_subset, type = "response")
      # Weight source to look like target: (1-p)/p
      ps_bounded <- pmax(0.01, pmin(0.99, ps_source))
      weights <- (1 - ps_bounded) / ps_bounded
      weights <- weights / sum(weights) * n_source
    } else {
      weights <- rep(1, n_source)
    }
    
  } else if (method == "calibration") {
    # Calibration weighting (entropy balancing approximation)
    # Target: match means of transport_vars
    target_means <- colMeans(target_subset[, transport_vars, drop = FALSE])
    source_mat <- as.matrix(source_subset[, transport_vars, drop = FALSE])
    
    # Simple iterative calibration
    weights <- rep(1, n_source)
    for (iter in 1:20) {
      for (j in seq_along(transport_vars)) {
        current_mean <- weighted.mean(source_mat[, j], weights)
        if (abs(current_mean) > 1e-10) {
          weights <- weights * (target_means[j] / current_mean)
        }
      }
      weights <- weights / sum(weights) * n_source
    }
    
  } else if (method == "sace") {
    # Sample average: no reweighting (assumes random sample)
    weights <- rep(1, n_source)
  }
  
  # Compute effective sample size
  ess <- sum(weights)^2 / sum(weights^2)
  
  # Check for extreme weights
  extreme_weights <- sum(weights > 10 * mean(weights)) / n_source
  
  list(
    weights = weights,
    ess = ess,
    extreme_weights = extreme_weights
  )
}

#' @keywords internal
.estimate_source_ate <- function(spec) {
  # Simple difference in means or regression
  data <- spec$data
  A <- data[[spec$treatment]]
  Y <- data[[spec$outcome]]
  
  if (!is.null(spec$covariates)) {
    formula_str <- paste(spec$outcome, "~", spec$treatment, "+",
                         paste(spec$covariates, collapse = " + "))
    model <- lm(as.formula(formula_str), data = data)
    coef(model)[spec$treatment]
  } else {
    mean(Y[A == 1]) - mean(Y[A == 0])
  }
}

#' @keywords internal
.estimate_transported_ate <- function(spec, weights) {
  data <- spec$data
  A <- data[[spec$treatment]]
  Y <- data[[spec$outcome]]
  
  # Weighted difference in means
  w_treated <- weights[A == 1]
  w_control <- weights[A == 0]
  
  weighted_mean_treated <- sum(Y[A == 1] * w_treated) / sum(w_treated)
  weighted_mean_control <- sum(Y[A == 0] * w_control) / sum(w_control)
  
  weighted_mean_treated - weighted_mean_control
}

#' @keywords internal
.compute_transport_delta <- function(cov_shift, extreme_weights, n_source, n_target) {
  # Combine covariate shift and weight instability
  mean_shift <- mean(cov_shift$shift_metric, na.rm = TRUE)
  
  # Penalty for sample size imbalance
  size_penalty <- abs(log(n_source / n_target)) / 10
  
  # Penalty for extreme weights
  weight_penalty <- extreme_weights * 0.5
  
  # Combined deficiency (heuristic)
  delta <- min(1, mean_shift * 0.5 + size_penalty + weight_penalty)
  
  delta
}

#' Print method for transport_deficiency
#' @param x A transport_deficiency object
#' @param ... Additional arguments (unused)
#' @export
print.transport_deficiency <- function(x, ...) {
  cli::cli_h1("Transport Deficiency Analysis")
  cli::cli_text("Method: {.val {x$method}}")
  cli::cli_text("Source n: {.val {x$n_source}} | Target n: {.val {x$n_target}}")
  cli::cli_text("")
  
  cli::cli_h2("Transport Deficiency")
  cat(sprintf("  delta_transport:           %.4f\n", x$delta_transport))
  if (!is.na(x$se)) {
    cat(sprintf("  Standard error:        %.4f\n", x$se))
    cat(sprintf("  95%% CI:               [%.4f, %.4f]\n", x$ci[1], x$ci[2]))
  }
  cat(sprintf("  Effective sample size: %.1f (%.1f%% of source)\n", 
              x$effective_sample_size,
              100 * x$effective_sample_size / x$n_source))
  
  cli::cli_h2("Covariate Shift")
  print(x$covariate_shift)
  
  cli::cli_h2("Effect Estimates")
  cat(sprintf("  ATE (source):     %.4f\n", x$ate_source))
  cat(sprintf("  ATE (transport):  %.4f\n", x$ate_target))
  
  invisible(x)
}

#' Plot method for transport_deficiency
#' @param x A transport_deficiency object
#' @param type Character: type of plot ("shift" or "weights")
#' @param ... Additional arguments passed to plot
#' @export
plot.transport_deficiency <- function(x, type = c("shift", "weights"), ...) {
  type <- match.arg(type)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggplot2} is required for plotting")
  }
  
  if (type == "shift") {
    df <- x$covariate_shift
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = reorder(variable, shift_metric), 
                                           y = shift_metric, fill = severity)) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "Covariate Shift: Source vs Target",
        x = "Variable",
        y = "Shift Metric (SMD / TV)",
        fill = "Severity"
      ) +
      ggplot2::scale_fill_manual(values = c(
        "low" = "forestgreen",
        "moderate" = "orange",
        "high" = "coral",
        "severe" = "darkred"
      )) +
      ggplot2::geom_hline(yintercept = c(0.1, 0.25), linetype = "dashed", alpha = 0.5) +
      ggplot2::theme_minimal()
    
  } else if (type == "weights") {
    df <- data.frame(weight = x$weights)
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = weight)) +
      ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "white") +
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
      ggplot2::labs(
        title = "Transport Weight Distribution",
        subtitle = sprintf("ESS: %.1f / %d (%.1f%%)", 
                          x$effective_sample_size, x$n_source,
                          100 * x$effective_sample_size / x$n_source),
        x = "Weight",
        y = "Count"
      ) +
      ggplot2::theme_minimal()
  }
  
  p
}
