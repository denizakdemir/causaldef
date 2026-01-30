# =============================================================================
# causal_spec() - Main Constructor
# =============================================================================

#' Create a Causal Problem Specification
#'
#' Defines the causal inference problem including treatment, outcome, covariates,
#' and optional diagnostic variables.
#'
#' @param data A data.frame or data.table containing the analysis data
#' @param treatment Character: name of the treatment variable
#' @param outcome Character: name of the outcome variable
#' @param covariates Character vector: names of adjustment covariates (NULL for none)
#' @param negative_control Character: name of negative control outcome (optional)
#' @param instrument Character: name of instrumental variable (optional)
#' @param estimand Character: target estimand, one of "ATE", "ATT", or "ATC"
#' @param outcome_type Character: type of outcome, one of "continuous", "binary", "count"
#' @param na.action Function: how to handle missing values (default: na.omit)
#'
#' @return Object of class "causal_spec".
#'   \strong{Note}: While `causal_spec` can describe categorical or continuous treatments, 
#'   current downstream estimation functions (`estimate_deficiency`, `estimate_effect`) 
#'   require **binary** treatments.
#'
#' @examples
#' # Create sample data
#' n <- 200
#' W <- rnorm(n)
#' A <- rbinom(n, 1, plogis(0.5 * W))
#' Y <- 1 + 2 * A + W + rnorm(n)
#' df <- data.frame(W = W, A = A, Y = Y)
#'
#' # Create causal specification
#' spec <- causal_spec(
#'   data = df,
#'   treatment = "A",
#'   outcome = "Y",
#'   covariates = "W"
#' )
#'
#' print(spec)
#'
#' @seealso [estimate_deficiency()], [nc_diagnostic()], [policy_regret_bound()]
#' @export
causal_spec <- function(data, treatment, outcome, covariates = NULL,
                        negative_control = NULL, instrument = NULL,
                        estimand = c("ATE", "ATT", "ATC"),
                        outcome_type = c("continuous", "binary", "count"),
                        na.action = na.omit) {
  
  # Input validation
  checkmate::assert_data_frame(data, min.rows = 10)
  checkmate::assert_string(treatment)
  checkmate::assert_string(outcome)
  checkmate::assert_character(covariates, null.ok = TRUE)
  checkmate::assert_string(negative_control, null.ok = TRUE)
  checkmate::assert_string(instrument, null.ok = TRUE)
  estimand <- match.arg(estimand)
  outcome_type <- match.arg(outcome_type)
  
  # Check variables exist
  all_vars <- c(treatment, outcome, covariates, negative_control, instrument)
  all_vars <- all_vars[!is.null(all_vars)]
  missing_vars <- setdiff(all_vars, names(data))
  
  if (length(missing_vars) > 0) {
    .msg_error(paste("Variables not found in data:", paste(missing_vars, collapse = ", ")))
  }
  
  # Handle missing values
  complete_data <- na.action(data[, all_vars, drop = FALSE])
  n_dropped <- nrow(data) - nrow(complete_data)
  
  if (n_dropped > 0) {
    .msg_warning(paste(n_dropped, "observations dropped due to missing values"))
  }
  
  # Infer treatment type
  treatment_type <- .infer_treatment_type(complete_data[[treatment]])
  
  # Check positivity (for binary treatment)
  if (treatment_type == "binary" && !is.null(covariates)) {
    .check_positivity(complete_data, treatment, covariates)
  }
  
  # Construct object
  spec <- new_causal_spec(
    data = complete_data,
    treatment = treatment,
    outcome = outcome,
    covariates = covariates,
    negative_control = negative_control,
    instrument = instrument,
    estimand = estimand,
    outcome_type = outcome_type,
    treatment_type = treatment_type
  )
  
  .msg_success(paste0("Created causal specification: n=", spec$n, ", ", 
                      length(covariates), " covariate(s)"))
  
  validate_causal_spec(spec)
}

# =============================================================================
# Helper Functions
# =============================================================================

#' Infer Treatment Type
#' @keywords internal
.infer_treatment_type <- function(treatment_var) {
  unique_vals <- unique(treatment_var)
  n_unique <- length(unique_vals)
  
  if (n_unique == 2) {
    return("binary")
  } else if (n_unique <= 10 && all(unique_vals == floor(unique_vals))) {
    return("categorical")
  } else {
    return("continuous")
  }
}

#' Check Positivity Assumption
#' @keywords internal
.check_positivity <- function(data, treatment, covariates, threshold = 0.01) {
  # Simple check: ensure treatment probabilities are bounded away from 0 and 1
  if (is.null(covariates) || length(covariates) == 0) {
    return(invisible(NULL))
  }
  
  # Fit propensity model
  formula_str <- paste(treatment, "~", paste(covariates, collapse = " + "))
  ps_model <- tryCatch(
    glm(as.formula(formula_str), data = data, family = binomial()),
    error = function(e) NULL
  )
  
  if (is.null(ps_model)) {
    return(invisible(NULL))
  }
  
  ps <- predict(ps_model, type = "response")
  
  # Check for extreme propensity scores
  if (any(ps < threshold) || any(ps > (1 - threshold))) {
    n_extreme <- sum(ps < threshold | ps > (1 - threshold))
    .msg_warning(
      paste(n_extreme, "observations have extreme propensity scores")
    )
  }
  
  invisible(NULL)
}
