# =============================================================================
# S3 Class Definitions for causaldef
# =============================================================================

# =============================================================================
# Class: causal_spec - Base causal specification
# =============================================================================

#' @title Create a New Causal Specification Object
#' @description Internal constructor for causal_spec class
#' @keywords internal
#' @noRd
new_causal_spec <- function(data, treatment, outcome, covariates = NULL,
                            negative_control = NULL, instrument = NULL,
                            estimand = "ATE", outcome_type = "continuous",
                            treatment_type = "binary", ...) {
  structure(
    list(
      data = data,
      treatment = treatment,
      outcome = outcome,
      covariates = covariates,
      negative_control = negative_control,
      instrument = instrument,
      estimand = estimand,
      outcome_type = outcome_type,
      n = nrow(data),
      treatment_type = treatment_type,
      ...
    ),
    class = "causal_spec"
  )
}

#' @title Validate Causal Specification
#' @description Validates a causal_spec object
#' @param x A causal_spec object
#' @return The validated object (invisibly)
#' @export
validate_causal_spec <- function(x) {
  checkmate::assert_class(x, "causal_spec")
  checkmate::assert_data_frame(x$data, min.rows = 10)
  checkmate::assert_string(x$treatment)
  
  if (inherits(x, "causal_spec_survival")) {
    checkmate::assert_string(x$time)
    checkmate::assert_string(x$event)
    checkmate::assert_subset(c(x$treatment, x$time, x$event, x$covariates),
                           names(x$data))
  } else {
    checkmate::assert_string(x$outcome)
    checkmate::assert_subset(c(x$treatment, x$outcome, x$covariates),
                           names(x$data))
  }
  invisible(x)
}

# =============================================================================
# Class: causal_spec_survival - Survival causal specification
# =============================================================================

#' @title Create a New Survival Causal Specification
#' @description Internal constructor for causal_spec_survival class
#' @keywords internal
#' @noRd
new_causal_spec_survival <- function(data, treatment, time, event,
                                     covariates = NULL, competing_event = NULL,
                                     estimand = "ATE", horizon = NULL,
                                     treatment_type = "binary", ...) {
  structure(
    list(
      data = data,
      treatment = treatment,
      time = time,
      event = event,
      covariates = covariates,
      competing_event = competing_event,
      estimand = estimand,
      horizon = horizon,
      n = nrow(data),
      treatment_type = treatment_type,
      n_events = sum(data[[event]] == 1),
      max_time = max(data[[time]]),
      ...
    ),
    class = c("causal_spec_survival", "causal_spec")
  )
}

# =============================================================================
# Class: deficiency - Le Cam deficiency estimates
# =============================================================================

#' @title Create a New Deficiency Object
#' @description Internal constructor for deficiency class
#' @keywords internal
#' @noRd
new_deficiency <- function(estimates, se = NULL, ci = NULL, method,
                           kernel = NULL, boot_samples = NULL, spec = NULL, ...) {
  structure(
    list(
      estimates = estimates,
      se = se,
      ci = ci,
      method = method,
      kernel = kernel,
      boot_samples = boot_samples,
      spec = spec,
      ...
    ),
    class = "deficiency"
  )
}

# =============================================================================
# Class: nc_diagnostic - Negative control diagnostic results
# =============================================================================

#' @title Create a New Negative Control Diagnostic Object
#' @description Internal constructor for nc_diagnostic class
#' @keywords internal
#' @noRd
new_nc_diagnostic <- function(delta_nc, delta_bound, falsified, p_value,
                              kappa, kernel, ...) {
  structure(
    list(
      delta_nc = delta_nc,
      delta_bound = delta_bound,
      falsified = falsified,
      p_value = p_value,
      kappa = kappa,
      kernel = kernel,
      ...
    ),
    class = "nc_diagnostic"
  )
}

# =============================================================================
# Class: policy_bound - Policy regret bounds
# =============================================================================

#' @title Create a New Policy Bound Object
#' @description Internal constructor for policy_bound class
#' @keywords internal
#' @noRd
new_policy_bound <- function(regret_bound, safety_floor, delta,
                             utility_range, obs_regret = NULL, 
                             all_estimates = NULL, ...) {
  structure(
    list(
      regret_bound = regret_bound,
      safety_floor = safety_floor,
      delta = delta,
      utility_range = utility_range,
      M = diff(utility_range),
      obs_regret = obs_regret,
      all_estimates = all_estimates,
      ...
    ),
    class = "policy_bound"
  )
}

# =============================================================================
# Class: confounding_frontier - Sensitivity analysis grid
# =============================================================================

#' @title Create a New Confounding Frontier Object
#' @description Internal constructor for confounding_frontier class
#' @keywords internal
#' @noRd
new_confounding_frontier <- function(grid, frontier = NULL, model, params) {
  structure(
    list(
      grid = grid,
      frontier = frontier,
      model = model,
      params = params
    ),
    class = "confounding_frontier"
  )
}

# =============================================================================
# Class: data_audit_report - Data integrity audit results
# =============================================================================

#' @title Create a New Data Audit Report Object
#' @description Internal constructor for data_audit_report class
#' @keywords internal
#' @noRd
new_data_audit_report <- function(issues, recommendations, summary_stats,
                                   treatment, outcome, spec = NULL) {
  structure(
    list(
      issues = issues,
      recommendations = recommendations,
      summary_stats = summary_stats,
      treatment = treatment,
      outcome = outcome,
      spec = spec
    ),
    class = "data_audit_report"
  )
}
