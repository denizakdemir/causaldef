# =============================================================================
# Effect Estimation using Causal Kernels
# =============================================================================

#' Estimate Causal Effects from Deficiency Objects
#'
#' Uses the fitted causal kernels (from `estimate_deficiency`) to estimate
#' the causal effect (ATE, RMST, etc.).
#'
#' @param object A `deficiency` object returned by `estimate_deficiency()`
#' @param ... Additional arguments passed to specific estimators
#'
#' @return A list containing:
#'   \itemize{
#'     \item estimand: The estimated effect (numeric)
#'     \item type: Type of estimand (e.g., "ATE", "RMST Diff")
#'     \item method: Method used (e.g., "iptw")
#'     \item contrast: Contrast levels
#'     \item ci: Confidence interval (if bootstrap supported - future work)
#'     \item curves: (Survival only) Adjusted survival curves
#'   }
#' 
#' @examples
#' \dontrun{
#' spec <- causal_spec(df, "A", "Y", "W")
#' def <- estimate_deficiency(spec, methods = "iptw")
#' effect <- estimate_effect(def)
#' print(effect)
#' }
#' 
#' @export
estimate_effect <- function(object, ...) {
  UseMethod("estimate_effect")
}

#' @rdname estimate_effect
#' @param method Character: which adjustment method to use. If `NULL`, defaults to
#'   the method with the lowest estimated deficiency.
#' @param target_method Deprecated alias for `method`.
#' @param contrast Vector: c(treated, control) values for contrast. 
#'   Defaults to c(treatment_value, control_value).
#' @export
estimate_effect.deficiency <- function(object, target_method = NULL,
                                       method = NULL,
                                       contrast = NULL, ...) {

  # Backward/forward compatible method selection:
  # - Prefer `method` (public API, consistent with other package functions)
  # - Support `target_method` for older code
  if (!is.null(method)) {
    if (!is.null(target_method)) {
      cli::cli_abort("Specify only one of {.arg method} or {.arg target_method}.")
    }
    target_method <- method
  }
  
  # Select method
  if (is.null(target_method)) {
    # Default to method with lowest deficiency
    target_method <- names(which.min(object$estimates))
    if (length(target_method) == 0) target_method <- object$method[1]
    .msg_info(paste("Using method with lowest deficiency:", target_method))
  }
  
  if (!target_method %in% names(object$kernel)) {
    .msg_error(paste("Method", target_method, "not found in deficiency object."))
  }
  
  kernel <- object$kernel[[target_method]]
  spec <- object$spec
  data <- spec$data
  
  # Determine contrast levels
  if (is.null(contrast)) {
    if (is.factor(data[[spec$treatment]])) {
      lvls <- levels(data[[spec$treatment]])
      treated <- kernel$treatment_value
      control <- setdiff(lvls, treated)[1] # Pick first non-treated
      if (is.na(control)) control <- lvls[1] # Fallback
    } else {
      treated <- 1
      control <- 0
    }
    contrast <- c(treated, control)
  }
  
  # Dispatch based on spec type
  if (inherits(spec, "causal_spec_survival")) {
    .estimate_effect_survival(spec, kernel, target_method, contrast, ...)
  } else {
    .estimate_effect_standard(spec, kernel, target_method, contrast, ...)
  }
}

.estimate_effect_standard <- function(spec, kernel, method, contrast, ...) {
  data <- spec$data
  treated_val <- contrast[1]
  control_val <- contrast[2]
  
  if (method == "iptw") {
    # Calculate IPW estimates for Treated and Control
    # Note: Kernel usually stores weights for the "Target" treatment value (treated_val)
    # But IPTW weights are typically 1/PS for Treated and 1/(1-PS) for Control
    
    weights <- kernel$weights # Assuming these are the proper IPW weights for the whole sample
    
    # Estimate E[Y(1)] and E[Y(0)]
    Y <- data[[spec$outcome]]
    A <- data[[spec$treatment]]
    
    # We need to ensure we identify the groups correctly
    is_treated <- A == treated_val
    is_control <- A == control_val
    
    # Weighted means
    # Start with robust checks
    if (sum(is_treated) == 0 || sum(is_control) == 0) {
      .msg_error("Treatment or control group empty in effect estimation")
    }
    
    mu1 <- weighted.mean(Y[is_treated], weights[is_treated], na.rm = TRUE)
    mu0 <- weighted.mean(Y[is_control], weights[is_control], na.rm = TRUE)
    
    ate <- mu1 - mu0
    
    structure(list(
      estimate = ate,
      mu1 = mu1,
      mu0 = mu0,
      type = "ATE",
      method = method,
      contrast = contrast
    ), class = "causal_effect")
    
  } else if (method == "aipw") {
    # AIPW provides direct estimates
    # kernel$mu_a is typically E[Y|A=a, W] for the *target* treatment value
    # To get ATE, we need E[Y(1)] and E[Y(0)].
    # The current deficiency implementation for AIPW calculates mu_a for *treatment_value*
    
    # Limitation: The stored kernel might only have info for ONE treatment level if not fully stored
    # But for AIPW we likely fitted a model `out_model`
    
    if (is.null(kernel$out_model)) {
      .msg_error("AIPW outcome model not found in kernel")
    }
    
    # Predict for both
    data_1 <- data
    data_1[[spec$treatment]] <- treated_val
    mu1_ind <- predict(kernel$out_model, newdata = data_1)
    
    data_0 <- data
    data_0[[spec$treatment]] <- control_val
    mu0_ind <- predict(kernel$out_model, newdata = data_0)
    
    # DR correction
    # E[Y(a)]_dr = mean( (I(A=a)/P(A=a|W)) * (Y - mu(a,W)) + mu(a,W) )
    
    A <- data[[spec$treatment]]
    Y <- data[[spec$outcome]]
    ps <- kernel$ps
    
    # Weights for propensity adjustment
    w1 <- 1/ps
    w0 <- 1/(1-ps) # Assuming binary PS
    
	    is_treated <- A == treated_val
	    is_control <- A == control_val
	
	    dr1 <- mean( (as.numeric(is_treated)/ps) * (Y - mu1_ind) + mu1_ind )
	    dr0 <- mean( (as.numeric(is_control)/(1-ps)) * (Y - mu0_ind) + mu0_ind )
    
    ate <- dr1 - dr0
    
    structure(list(
      estimate = ate,
      mu1 = dr1,
      mu0 = dr0,
      type = "ATE (Doubly Robust)",
      method = method,
      contrast = contrast
    ), class = "causal_effect")
    
  } else {
     # Unadjusted
    Y <- data[[spec$outcome]]
    A <- data[[spec$treatment]]
    mu1 <- mean(Y[A == treated_val], na.rm=TRUE)
    mu0 <- mean(Y[A == control_val], na.rm=TRUE)
    ate <- mu1 - mu0
    
    structure(list(
      estimate = ate,
      mu1 = mu1,
      mu0 = mu0,
      type = "ATE (Unadjusted)",
      method = method,
      contrast = contrast
    ), class = "causal_effect")
  }
}

.estimate_effect_survival <- function(spec, kernel, method, contrast, ...) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    .msg_error("Package 'survival' needed for survival effect estimation.")
  }
  
  data <- spec$data
  treated_val <- contrast[1]
  control_val <- contrast[2]
  
  # Setup Surv object
  # If competing risks, we currently focus on the primary event of interest
  # censoring others.
  
  if (!is.null(spec$competing_event)) {
    # Competing Risks: Usually requires specialized analysis (CumInc)
    # For now, we use Cause-Specific Hazard approach (censoring competing events)
    # which is valid for estimating the specific hazard, but KM represents hypothetical worlds
    # where competing events are eliminated, which requires independence assumptions.
    # To keep it simple and robust, we treat competing events as censored (CSH).
    status <- ifelse(data[[spec$event]] == 1, 1, 0)
  } else {
    status <- data[[spec$event]]
  }
  
  surv_obj <- survival::Surv(data[[spec$time]], status)
  
  # Weights
  if (method %in% c("iptw", "aipw") && !is.null(kernel$weights)) {
     weights_vec <- kernel$weights
  } else {
     weights_vec <- rep(1, nrow(data))
  }
  
  # Weighted KM
  # Creates a temporary data frame to avoid scoping issues with survfit
  df_calc <- data
  df_calc$.surv <- surv_obj
  df_calc$.weights <- weights_vec
  
  formula <- as.formula(paste(".surv ~", spec$treatment))
  
  # call survfit. weights must be in the data frame or globally accessible
  km_fit <- survival::survfit(formula, data = df_calc, weights = .weights)
  
  # Extract estimates at horizon if RMST
  if (spec$estimand == "RMST") {
    if (is.null(spec$horizon)) .msg_error("Horizon required for RMST")
    
    rmst_res <- summary(km_fit, rmean = spec$horizon)
    
    # Validated RMST Extraction
    if (!inherits(rmst_res, "summary.survfit")) {
      .msg_error("Failed to generate valid summary from survival fit")
    }
    
    # The summary(fit, rmean=...) object contains a matrix named 'table'
    tbl <- rmst_res$table
    
    if (is.null(tbl)) {
      .msg_error("RMST table not found in survival summary")
    }
    
    # Identify relevant column for RMST (Restricted Mean Time)
    # Common names: "*rmean", "rmean"
    rmean_col_idx <- grep("rmean", colnames(tbl))
    if (length(rmean_col_idx) == 0) {
       # Fallback: if not named, it is usually the last or specific numeric column.
       # Standard columns: records, n.max, n.start, events, *rmean, *se(rmean)
       # We look for a column that could be it.
       # However, we should be strict.
       .msg_error("Could not identify 'rmean' column in survival output. Check survival package version.")
    } else {
       rmean_col <- rmean_col_idx[1]
    }

    # Helper to clean row names for matching
    clean_names <- function(x) gsub(".*=", "", x) 
    
    get_rmst <- function(val) {
      row_names <- rownames(tbl)
      val_str <- as.character(val)
      
      # Try exact match (e.g. "treatment=1")
      idx <- grep(paste0(spec$treatment, "=", val_str), row_names, fixed = TRUE)
      
      # Try loose match (e.g. "1")
      if (length(idx) == 0) {
         idx <- which(clean_names(row_names) == val_str)
      }
      
      if (length(idx) != 1) {
         return(NA) # Ambiguous or not found
      }
      
      return(tbl[idx, rmean_col])
    }
    
    rmst1 <- get_rmst(treated_val)
    rmst0 <- get_rmst(control_val)
    
    est <- rmst1 - rmst0
    
    structure(list(
      estimate = est,
      rmst1 = rmst1,
      rmst0 = rmst0,
      fit = km_fit,
      horizon = spec$horizon,
      type = "RMST Difference",
      method = method,
      contrast = contrast
    ), class = "causal_effect")
    
  } else {
    # Default to Survival Probability at specific time?
    # Or just return the fit
    structure(list(
      fit = km_fit,
      type = "Survival Curve",
      method = method,
      contrast = contrast
    ), class = "causal_effect")
  }
}

#' @export
print.causal_effect <- function(x, ...) {
  cat("\n-- Causal Effect Estimate ----------------------\n")
  cat("Method:   ", x$method, "\n")
  cat("Type:     ", x$type, "\n")
  cat("Contrast: ", paste(x$contrast, collapse = " vs "), "\n")
  
  if (!is.null(x$estimate)) {
    cat("Estimate: ", round(x$estimate, 4), "\n")
  }
  
  if (!is.null(x$horizon)) {
    cat("Horizon:  ", x$horizon, "\n")
  }
  
  invisible(x)
}
