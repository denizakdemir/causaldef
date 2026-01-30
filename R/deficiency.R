# =============================================================================
# estimate_deficiency() - Main Deficiency Estimation
# =============================================================================

#' Estimate Le Cam Deficiency
#'
#' Computes the deficiency \eqn{\delta(E_{obs}, E_{do})} between observational and
#' interventional experiments under various adjustment strategies.
#'
#' @param spec A causal_spec object
#' @param methods Character vector: adjustment methods to compare
#'   \itemize{
#'     \item "unadjusted": No adjustment (baseline)
#'     \item "iptw": Inverse probability weighting
#'     \item "aipw": Augmented IPW (doubly robust)
#'   }
#' @param treatment_value Numeric: intervention value for do(A=a)
#' @param n_boot Integer: bootstrap replicates for CI (0 = none)
#' @param ci_level Numeric: confidence level (default 0.95)
#' @param verbose Logical: show progress
#'
#' @return Object of class "deficiency" containing:
#'   \itemize{
#'     \item estimates: Named vector of \eqn{\delta} estimates per method
#'     \item se: Standard errors (if n_boot > 0)
#'     \item ci: Confidence intervals
#'     \item method: Methods used
#'     \item kernel: Fitted kernels for each method
#'   }
#'
#' @details
#' The Le Cam deficiency quantifies the information gap:
#' \deqn{\delta(\mathcal{E}_{obs}, \mathcal{E}_{do}) =
#'       \inf_K \sup_\theta ||KP^{obs}_\theta - P^{do}_\theta||_{TV}}
#'
#' When \eqn{\delta = 0}, perfect causal identification is possible.
#' When \eqn{\delta > 0}, some information loss is unavoidable.
#'
#' @examples
#' # Create sample data
#' n <- 200
#' W <- rnorm(n)
#' A <- rbinom(n, 1, plogis(0.5 * W))
#' Y <- 1 + 2 * A + W + rnorm(n)
#' df <- data.frame(W = W, A = A, Y = Y)
#'
#' spec <- causal_spec(df, "A", "Y", "W")
#' results <- estimate_deficiency(spec, methods = c("unadjusted", "iptw"), n_boot = 50)
#' print(results)
#'
#' @seealso [causal_spec()], [nc_diagnostic()], [policy_regret_bound()]
#' @export
estimate_deficiency <- function(spec, methods = c("iptw", "aipw"),
                                treatment_value = NULL, n_boot = 200,
                                ci_level = 0.95, verbose = TRUE) {
  
  # Validation
  validate_causal_spec(spec)
  valid_methods <- c("unadjusted", "iptw", "aipw", "matching")
  checkmate::assert_subset(methods, valid_methods)
  checkmate::assert_integerish(n_boot, lower = 0)
  checkmate::assert_number(ci_level, lower = 0.5, upper = 0.999)
  
  # Enforce binary treatment restriction (current implementation limitation)
  if (spec$treatment_type != "binary") {
    cli::cli_abort(c(
      "Unsupported treatment type: {.val {spec$treatment_type}}",
      "x" = "The current version of {.pkg causaldef} only supports binary treatments for deficiency estimation.",
      "i" = "The underlying estimators (IPTW/AIPW) currently assume a logistic propensity model.",
      "i" = "Please dichotomize your treatment variable (e.g., High vs. Low) or subset your data to two groups."
    ))
  }
  
  # Infer treatment value if not provided
  if (is.null(treatment_value)) {
    tr_var <- spec$data[[spec$treatment]]
    if (is.factor(tr_var)) {
      treatment_value <- levels(tr_var)[2] # Default to second level (usually treated)
      .msg_info(paste("Inferred treatment value:", treatment_value))
    } else {
      treatment_value <- 1
    }
  }

  # Estimate deficiency for each method
  results <- lapply(methods, function(m) {
    if (verbose) .msg_info(paste("Estimating deficiency:", m))
    .estimate_deficiency_single(spec, m, treatment_value, n_boot, ci_level)
  })
  names(results) <- methods
  
  # Combine into single deficiency object
  estimates <- vapply(results, function(r) r$estimate, numeric(1))
  names(estimates) <- methods
  
  se <- if (n_boot > 0) {
    vapply(results, function(r) r$se, numeric(1))
  } else {
    rep(NA_real_, length(methods))
  }
  names(se) <- methods
  
  ci <- if (n_boot > 0) {
    do.call(rbind, lapply(results, function(r) r$ci))
  } else {
    matrix(NA_real_, nrow = length(methods), ncol = 2)
  }
  rownames(ci) <- methods
  colnames(ci) <- c("lower", "upper")
  
  kernels <- lapply(results, function(r) r$kernel)
  names(kernels) <- methods
  
  new_deficiency(
    estimates = estimates,
    se = se,
    ci = ci,
    method = methods,
    kernel = kernels,
    spec = spec
  )
}

# =============================================================================
# Internal: Single Method Deficiency Estimation
# =============================================================================

#' @keywords internal
.estimate_deficiency_single <- function(spec, method, treatment_value, 
                                        n_boot, ci_level) {
  
  data <- spec$data
  A <- data[[spec$treatment]]
  
  if (inherits(spec, "causal_spec_survival")) {
    Y <- data[[spec$time]]
  } else {
    Y <- data[[spec$outcome]]
  }
  W <- if (!is.null(spec$covariates)) {
    as.matrix(data[, spec$covariates, drop = FALSE])
  } else {
    NULL
  }
  
  # Build kernel and compute point estimate
  kernel <- .build_kernel(spec, method, treatment_value)
  estimate <- .compute_deficiency(spec, kernel, method)
  
  # Bootstrap if requested
  if (n_boot > 0) {
    boot_estimates <- replicate(n_boot, {
      boot_idx <- sample(nrow(data), replace = TRUE)
      boot_data <- data[boot_idx, , drop = FALSE]
      boot_spec <- spec
      boot_spec$data <- boot_data
      boot_spec$n <- nrow(boot_data)
      boot_kernel <- .build_kernel(boot_spec, method, treatment_value)
      .compute_deficiency(boot_spec, boot_kernel, method)
    })
    
    se <- sd(boot_estimates)
    alpha <- 1 - ci_level
    ci <- quantile(boot_estimates, probs = c(alpha/2, 1 - alpha/2))
  } else {
    se <- NA_real_
    ci <- c(NA_real_, NA_real_)
  }
  
  list(
    estimate = estimate,
    se = se,
    ci = ci,
    kernel = kernel
  )
}

# =============================================================================
# Internal: Kernel Building
# =============================================================================

#' @keywords internal
.build_kernel <- function(spec, method, treatment_value) {
  
  data <- spec$data
  A <- data[[spec$treatment]]
  W <- if (!is.null(spec$covariates)) {
    as.matrix(data[, spec$covariates, drop = FALSE])
  } else {
    NULL
  }
  
  switch(method,
    "unadjusted" = list(
      method = "unadjusted",
      treatment_value = treatment_value,
      weights = rep(1, nrow(data))
    ),
    "iptw" = {
      if (is.null(W)) {
        .msg_error("IPTW requires covariates in causal_spec")
      }
      # Fit propensity score model
      ps_formula <- as.formula(paste(spec$treatment, "~", 
                                     paste(spec$covariates, collapse = " + ")))
      ps_model <- glm(ps_formula, data = data, family = binomial())
      ps <- predict(ps_model, type = "response")
      
      # IPW weights
      weights <- ifelse(A == treatment_value, 
                        1 / ps, 
                        1 / (1 - ps))
      weights <- weights / sum(weights) * length(weights)  # Normalize
      
      list(
        method = "iptw",
        treatment_value = treatment_value,
        ps_model = ps_model,
        ps = ps,
        weights = weights
      )
    },
    "aipw" = {
      if (is.null(W)) {
        .msg_error("AIPW requires covariates in causal_spec")
      }
      
      # Propensity score
      ps_formula <- as.formula(paste(spec$treatment, "~", 
                                     paste(spec$covariates, collapse = " + ")))
      ps_model <- glm(ps_formula, data = data, family = binomial())
      ps <- predict(ps_model, type = "response")
      
      # Outcome model
      if (inherits(spec, "causal_spec_survival")) {
        y_name <- spec$time
      } else {
        y_name <- spec$outcome
      }
      
      out_formula <- as.formula(paste(y_name, "~ (",
                                      spec$treatment, ") + (",
                                      paste(spec$covariates, collapse = " + "), ")"))
      out_model <- lm(out_formula, data = data)
      
      # Predicted outcomes under treatment
      data_treated <- data
      data_treated[[spec$treatment]] <- treatment_value
      mu_a <- predict(out_model, newdata = data_treated)
      
      list(
        method = "aipw",
        treatment_value = treatment_value,
        ps_model = ps_model,
        out_model = out_model,
        ps = ps,
        mu_a = mu_a
      )
    },
    "matching" = {
      # Simple placeholder - would use MatchIt in full implementation
      list(
        method = "matching",
        weights = rep(1, nrow(data))
      )
    }
  )
}

# =============================================================================
# Internal: Deficiency Computation
# =============================================================================

#' @keywords internal
.compute_deficiency <- function(spec, kernel, method) {
  
  data <- spec$data
  A <- data[[spec$treatment]]
  n <- nrow(data)
  
  # For deficiency estimation, we need a common ground to compare distributions.
  # The Propensity Score (PS) is the standard summary score for confounding.
  # We estimate the TV distance between the effective distributions implied by the method
  # and the target interventional distribution (approximated by the full sample marginal).
  
  # 1. Get Propensity Scores (PS)
  # If not computed in kernel (e.g. unadjusted), we must estimate them for the metric
  if (!is.null(kernel$ps)) {
    ps <- kernel$ps
  } else {
    # Fit nuisance PS model if not present
    if (!is.null(spec$covariates)) {
      ps_formula <- as.formula(paste(spec$treatment, "~", 
                                     paste(spec$covariates, collapse = " + ")))
      # Use tryCatch for robust automatic calculations
      ps_model <- tryCatch({
        glm(ps_formula, data = data, family = binomial())
      }, error = function(e) return(NULL))
      
      if (!is.null(ps_model)) {
        ps <- predict(ps_model, type = "response")
      } else {
        return(1.0) # Maximum deficiency if model fails
      }
    } else {
      # No covariates -> randomization assumption or uninformative
      ps <- rep(mean(A == kernel$treatment_value), n)
    }
  }
  
  # 2. Define Weights based on Method
  if (method == "unadjusted") {
    # Unadjusted compares Treated vs Control directly (or Treated vs Population)
    weights <- ifelse(A == kernel$treatment_value, 1, 0)
    # Normalize to sum to n
    weights <- weights / sum(weights) * n
    
  } else if (method %in% c("iptw", "aipw")) {
    # Use the weights defined by the method, OR derive them from PS
    if (!is.null(kernel$weights)) {
      weights <- kernel$weights
    } else {
      # AIPW uses same weights as IPTW for the deficiency metric (propensity overlap)
      # weights = 1/e for treated (if target), 1/(1-e) for control
      # This estimates the full population ATE.
      # If kernel$treatment_value corresponds to A=1, and we want ATE:
      # W = Z/e + (1-Z)/(1-e)
      # But wait, we want to construct the "Pseudo-Population" that mimics the TARGET.
      # For ATE, target is P(X).
      # The weighted sample should mimic P(X).
      # Standard IPW weights: w = 1/e(X) for A=1, w = 1/(1-e(X)) for A=0.
      
      # Ensure ps is bounded to avoid explosion
      ps_bounded <- pmax(0.001, pmin(0.999, ps))
      
      # We assume treatment_value is the '1' in binary for now.
      # If treatment variable matches treatment_value
      match_val <- A == kernel$treatment_value
      
      weights <- ifelse(match_val, 1/ps_bounded, 1/(1-ps_bounded))
    }
    
    # Normalize
    weights <- weights / sum(weights) * n
  } else {
    weights <- rep(1, n)
  }

  
  # 3. Compute TV Distance using Histogram approximation on Propensity Score
  # This provides a 1-dimensional bounded estimate of total variation.
  # TV(P, Q) = 0.5 * integral |dP - dQ|
  # Target Q is the distinct marginal distribution P(W) (weights = 1 for all)
  # Source P is the weighted distribution
  
  # Define bins (e.g., 20 bins or adaptive)
  n_bins <- min(50, floor(n / 10))
  if (n_bins < 5) n_bins <- 5
  
  # Bin PS scores
  breaks <- seq(0, 1, length.out = n_bins + 1)
  
  # H_obs: Weighted counts in each bin
  h_obs <- .weighted_hist(ps, weights, breaks)
  
  # H_target: Unweighted counts (Target distribution is the full population for ATE)
  # Note: For ATT, target would be different. Assuming ATE here.
  h_target <- .weighted_hist(ps, rep(1, n), breaks)
  
  # Normalize to densities (sum = 1)
  d_obs <- h_obs / sum(h_obs)
  d_target <- h_target / sum(h_target)
  
  # Calculate TV
  tv_dist <- 0.5 * sum(abs(d_obs - d_target))
  
  # 4. Refine for AIPW
  # AIPW should theoretically strictly improve on IPTW if the outcome model helps.
  # However, strictly speaking, deficiency is a property of the identification strategy.
  # If we trust the outcome model, deficiency drops. 
  # We use the residual confounding heuristic:
  # If outcome model is perfect, TV -> 0.
  # We can't know it's perfect. We take the min of the weighted TV estimates 
  # and the outcome-model residual balance (not implemented for simplicity).
  # For now, we report the TV distance of the weights, which acts as a "safety floor".
  # We do NOT arbitrarily discount AIPW. AIPW and IPTW share the same propensity metric
  # but AIPW is more robust. We enable the bound to reflect the propensity fit.
  
  # Ensure bounded [0, 1]
  delta <- max(0, min(1, tv_dist))
  return(delta)
}

# =============================================================================
# Internal: Weighted Histogram Helper
# =============================================================================

#' @keywords internal
.weighted_hist <- function(x, w, breaks) {
  # Assign x to bins
  bins <- cut(x, breaks, include.lowest = TRUE, labels = FALSE)
  
  # Sum weights per bin
  # We use tapply or aggregate.
  # Ensure all bins 1..n_bins are present
  bin_sums <- rep(0, length(breaks) - 1)
  
  if (length(x) > 0) {
    sums <- tapply(w, bins, sum)
    bin_sums[as.numeric(names(sums))] <- sums
  }
  
  bin_sums
}

# Note: The original SMD computation functions (.compute_smd, .compute_weighted_smd)
# have been removed as they provided heuristic proxies rather than theoretical bounds.

