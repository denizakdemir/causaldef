# =============================================================================
# Survival Causal Specification
# =============================================================================

#' Create a Survival Causal Specification
#'
#' Defines a causal inference problem with time-to-event outcomes.
#'
#' @param data A data.frame containing survival data
#' @param treatment Character: name of the treatment variable
#' @param time Character: name of the time-to-event variable
#' @param event Character: name of the event indicator (1 = event, 0 = censored).
#'   If the column is a factor/character (e.g., \code{"Death"}, \code{"Relapse"},
#'   \code{"Censored"}), it will be coerced to a 0/1 indicator by treating
#'   \code{"Censored"} (case-insensitive) and \code{"0"} as censored and all other
#'   non-missing values as events. For competing risks, prefer
#'   [causal_spec_competing()].
#' @param covariates Character vector: names of adjustment covariates
#' @param competing_event Character: name of competing event indicator (optional)
#' @param estimand Character: target estimand ("ATE", "RMST", "HR")
#' @param horizon Numeric: time horizon for RMST (required if estimand = "RMST")
#' @param na.action Function: how to handle missing values
#'
#' @return Object of class "causal_spec_survival"
#'
#' @details
#' For survival outcomes, the deficiency measures how well we can simulate
#' the interventional survival distribution from observational data.
#'
#' @section Estimands:
#' \describe{
#'   \item{ATE}{Average treatment effect on survival probability at horizon}
#'   \item{RMST}{Restricted mean survival time difference}
#'   \item{HR}{Hazard ratio (log scale)}
#' }
#'
#' @examples
#' # Simulate survival data
#' n <- 200
#' W <- rnorm(n)
#' A <- rbinom(n, 1, plogis(0.5 * W))
#' time <- rexp(n, rate = exp(-0.5 * A + 0.3 * W))
#' event <- rbinom(n, 1, 0.8)
#' df <- data.frame(W = W, A = A, time = time, event = event)
#'
#' spec <- causal_spec_survival(
#'   data = df,
#'   treatment = "A",
#'   time = "time",
#'   event = "event",
#'   covariates = "W",
#'   estimand = "RMST",
#'   horizon = 5
#' )
#'
#' @seealso [causal_spec()], [estimate_deficiency()]
#' @export
causal_spec_survival <- function(data, treatment, time, event,
                                 covariates = NULL, competing_event = NULL,
                                 estimand = c("ATE", "RMST", "HR"),
                                 horizon = NULL, na.action = na.omit) {
  
  # Input validation
  checkmate::assert_data_frame(data, min.rows = 10)
  checkmate::assert_string(treatment)
  checkmate::assert_string(time)
  checkmate::assert_string(event)
  checkmate::assert_character(covariates, null.ok = TRUE)
  checkmate::assert_string(competing_event, null.ok = TRUE)
  estimand <- match.arg(estimand)
  
  # Check variables exist
  all_vars <- c(treatment, time, event, covariates, competing_event)
  all_vars <- all_vars[!is.null(all_vars)]
  missing_vars <- setdiff(all_vars, names(data))
  
  if (length(missing_vars) > 0) {
    .msg_error(paste("Variables not found in data:", paste(missing_vars, collapse = ", ")))
  }
  
  # RMST requires horizon
  if (estimand == "RMST" && is.null(horizon)) {
    .msg_error("RMST estimand requires horizon to be specified")
  }
  
  # Handle missing values
  complete_data <- na.action(data[, all_vars, drop = FALSE])
  n_dropped <- nrow(data) - nrow(complete_data)
  
  if (n_dropped > 0) {
    .msg_warning(paste(n_dropped, "observations dropped due to missing values"))
  }
  
  # Infer treatment type
  treatment_type <- .infer_treatment_type(complete_data[[treatment]])
  
  # Coerce event indicator to a 0/1 numeric vector if needed.
  event_vec <- complete_data[[event]]
  if (is.factor(event_vec) || is.character(event_vec)) {
    event_chr <- tolower(as.character(event_vec))
    censored <- is.na(event_chr) | event_chr %in% c("censored", "censor", "0", "")
    non_censor_levels <- sort(unique(event_chr[!censored]))
    
    if (length(non_censor_levels) > 1) {
      .msg_warning(paste0(
        "Event column '", event, "' has multiple non-censor levels (",
        paste(shQuote(non_censor_levels), collapse = ", "), "). ",
        "Coercing to a binary indicator where any non-censored level is treated as an event. ",
        "For competing risks, use causal_spec_competing() or create an explicit 0/1 event indicator."
      ))
    }
    
    complete_data[[event]] <- as.integer(!censored)
  } else if (is.logical(event_vec)) {
    complete_data[[event]] <- as.integer(event_vec)
  } else if (is.numeric(event_vec) || is.integer(event_vec)) {
    # Accept 0/1, and coerce positive integers to 1 with a warning.
    non_missing <- event_vec[!is.na(event_vec)]
    if (length(non_missing) > 0 && !all(non_missing %in% c(0, 1))) {
      if (all(non_missing >= 0 & non_missing == floor(non_missing))) {
        .msg_warning(paste0(
          "Event column '", event, "' is not binary. ",
          "Coercing to 0/1 by treating values > 0 as events."
        ))
        complete_data[[event]] <- as.integer(event_vec > 0)
      } else {
        .msg_error(paste0(
          "Event column '", event, "' must be a 0/1 indicator (or a non-negative integer code)."
        ))
      }
    }
  } else {
    .msg_error(paste0(
      "Unsupported type for event column '", event, "'. ",
      "Provide a 0/1 indicator (numeric/logical) or a factor/character with a 'Censored' level."
    ))
  }
  
  # Competing event (if provided) should be 0/1. Coerce similarly.
  if (!is.null(competing_event)) {
    ce_vec <- complete_data[[competing_event]]
    if (is.factor(ce_vec) || is.character(ce_vec)) {
      ce_chr <- tolower(as.character(ce_vec))
      ce_censored <- is.na(ce_chr) | ce_chr %in% c("censored", "censor", "0", "")
      complete_data[[competing_event]] <- as.integer(!ce_censored)
    } else if (is.logical(ce_vec)) {
      complete_data[[competing_event]] <- as.integer(ce_vec)
    } else if (is.numeric(ce_vec) || is.integer(ce_vec)) {
      complete_data[[competing_event]] <- as.integer(ce_vec > 0)
    } else {
      .msg_error(paste0(
        "Unsupported type for competing_event column '", competing_event, "'. ",
        "Provide a 0/1 indicator."
      ))
    }
  }
  
  # Construct object
  spec <- new_causal_spec_survival(
    data = complete_data,
    treatment = treatment,
    time = time,
    event = event,
    covariates = covariates,
    competing_event = competing_event,
    estimand = estimand,
    horizon = horizon,
    treatment_type = treatment_type
  )
  
  .msg_success(paste0("Created survival causal specification: n=", spec$n, ", ", spec$n_events, " events"))
  
  spec
}

#' @export
print.causal_spec_survival <- function(x, ...) {
  cat("\n-- Survival Causal Specification ", paste(rep("-", 40), collapse = ""), "\n\n", sep = "")
  cat("* Treatment:", x$treatment, "(", x$treatment_type, ")\n")
  cat("* Time:", x$time, "(max =", round(x$max_time, 2), ")\n")
  cat("* Event:", x$event, "(", x$n_events, "events)\n")
  cat("* Covariates:", paste(x$covariates, collapse = ", "), "\n")
  cat("* Sample size:", x$n, "\n")
  cat("* Estimand:", x$estimand, "\n")
  
  if (!is.null(x$horizon)) {
    cat("* Horizon:", x$horizon, "\n")
  }
  if (!is.null(x$competing_event)) {
    cat("* Competing event:", x$competing_event, "\n")
  }
  cat("\n")
  invisible(x)
}
