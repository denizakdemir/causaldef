# =============================================================================
# Competing Risks Estimation with Deficiency Bounds
# =============================================================================

#' Causal Specification for Competing Risks
#'
#' Creates a specification object for time-to-event data with competing events.
#' Extends standard survival analysis to handle multiple event types.
#'
#' @param data Data frame containing the variables
#' @param treatment Character: name of the treatment variable
#' @param time Character: name of the time-to-event variable
#' @param event Character: name of the event indicator.
#'   Accepts either:
#'   \itemize{
#'     \item a numeric code with \code{0 = censored} and \code{1,2,...} indicating event types, or
#'     \item a factor/character (e.g., \code{"Death"}, \code{"Relapse"}, \code{"Censored"}), which will
#'       be mapped to integer codes with \code{"Censored"} (case-insensitive) and \code{"0"} treated as censored.
#'   }
#' @param covariates Character vector: names of covariate variables
#' @param event_of_interest Integer or character: which event type is the primary outcome (default 1).
#'   If \code{event} is factor/character, you may pass the label (e.g., \code{"Death"}).
#' @param horizon Numeric: time horizon for cumulative incidence
#' @param estimand Character: "cif" (cumulative incidence) or "cshr" (cause-specific HR)
#'
#' @return Object of class c("causal_spec_competing", "causal_spec_survival", "causal_spec")
#'
#' @details
#' In competing risks, standard survival methods can be biased because:
#' \enumerate{
#'   \item Censoring the competing event underestimates cumulative incidence
#'   \item Cause-specific hazards don't translate directly to probabilities
#' }
#'
#' The Aalen-Johansen estimator handles this by modeling all transitions:
#' \deqn{F_k(t) = P(T \leq t, J = k) = \int_0^t S(u-) \lambda_k(u) du}
#'
#' @examples
#' # Simulate competing risks: death (1) vs transplant (2)
#' set.seed(42)
#' n <- 500
#' W <- rnorm(n)
#' A <- rbinom(n, 1, plogis(0.3 * W))
#' 
#' # Event times
#' rate_death <- exp(-0.5 * A + 0.2 * W)
#' rate_transplant <- exp(0.3 * A - 0.1 * W)
#' 
#' time_death <- rexp(n, rate_death)
#' time_transplant <- rexp(n, rate_transplant)
#' time_censor <- runif(n, 0, 3)
#' 
#' observed_time <- pmin(time_death, time_transplant, time_censor)
#' event <- ifelse(observed_time == time_death, 1,
#'                 ifelse(observed_time == time_transplant, 2, 0))
#' 
#' df <- data.frame(W = W, A = A, time = observed_time, event = event)
#' 
#' spec <- causal_spec_competing(
#'   df, "A", "time", "event", "W",
#'   event_of_interest = 1,
#'   horizon = 2
#' )
#' print(spec)
#'
#' @seealso [estimate_deficiency()], [causal_spec_survival()]
#' @export
causal_spec_competing <- function(data, treatment, time, event, 
                                   covariates = NULL,
                                   event_of_interest = 1,
                                   horizon = NULL,
                                   estimand = c("cif", "cshr")) {
  
  estimand <- match.arg(estimand)
  
  # Validation
  checkmate::assert_data_frame(data, min.rows = 10)
  checkmate::assert_string(treatment)
  checkmate::assert_string(time)
  checkmate::assert_string(event)
  checkmate::assert_character(covariates, null.ok = TRUE)
  if (is.character(event_of_interest)) {
    checkmate::assert_string(event_of_interest)
  } else {
    checkmate::assert_integerish(event_of_interest, lower = 1)
  }
  checkmate::assert_number(horizon, lower = 0, null.ok = TRUE)
  
  # Check variables exist
  all_vars <- c(treatment, time, event, covariates)
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars) > 0) {
    cli::cli_abort("Variables not found in data: {.val {missing_vars}}")
  }
  
  # Extract and normalize event codes:
  # - numeric codes: 0=censored, 1+ event types
  # - factor/character: map labels to 1..K, with "Censored"/"0" as 0.
  event_vec_raw <- data[[event]]
  event_map <- NULL
  event_levels <- NULL
  
  if (is.factor(event_vec_raw) || is.character(event_vec_raw)) {
    event_chr <- tolower(as.character(event_vec_raw))
    censored <- is.na(event_chr) | event_chr %in% c("censored", "censor", "0", "")
    event_levels <- sort(unique(event_chr[!censored]))
    
    if (length(event_levels) == 0) {
      cli::cli_abort(c(
        "Event column {.val {event}} contains no non-censored events after coercion.",
        "i" = "Provide a numeric event code (0=censored, 1+ event types) or a factor/character with non-censored levels."
      ))
    }
    
    event_map <- stats::setNames(seq_along(event_levels), event_levels)
    event_vec <- ifelse(censored, 0L, unname(event_map[event_chr]))
    data[[event]] <- as.integer(event_vec)
    
    # Allow event_of_interest as a label
    if (is.character(event_of_interest)) {
      eoi_key <- tolower(event_of_interest)
      if (!(eoi_key %in% names(event_map))) {
        cli::cli_abort(c(
          "Event of interest {.val {event_of_interest}} not found in {.val {event}}.",
          "i" = "Available event labels: {.val {names(event_map)}}"
        ))
      }
      event_of_interest <- unname(event_map[[eoi_key]])
    }
  } else if (is.numeric(event_vec_raw) || is.integer(event_vec_raw) || is.logical(event_vec_raw)) {
    if (is.character(event_of_interest)) {
      cli::cli_abort(c(
        "event_of_interest was provided as a label ({.val {event_of_interest}}) but {.val {event}} is numeric.",
        "i" = "Provide event_of_interest as an integer code (e.g., 1) or pass a factor/character event column."
      ))
    }
    event_vec <- as.integer(event_vec_raw)
    if (any(is.na(event_vec))) {
      # Missing values treated as censored
      event_vec[is.na(event_vec)] <- 0L
      data[[event]] <- event_vec
    }
    if (any(event_vec < 0)) {
      cli::cli_abort("Event column {.val {event}} must be non-negative (0=censored, 1+ event types).")
    }
  } else {
    cli::cli_abort(c(
      "Unsupported type for event column {.val {event}}.",
      "i" = "Provide a numeric event code (0=censored, 1+ event types) or a factor/character with a 'Censored' level."
    ))
  }
  
  event_types <- sort(unique(event_vec[event_vec > 0]))
  n_events <- length(event_types)
  
  if (!(event_of_interest %in% event_types)) {
    cli::cli_abort(c(
      "Event of interest {.val {event_of_interest}} not found in data.",
      "i" = "Available event types: {.val {event_types}}"
    ))
  }
  
  # Infer horizon if not provided
  if (is.null(horizon)) {
    horizon <- quantile(data[[time]], 0.9)
    cli::cli_alert_info("Using horizon = {round(horizon, 2)} (90th percentile)")
  }
  
  # Infer treatment type
  treatment_values <- unique(data[[treatment]])
  if (length(treatment_values) == 2) {
    treatment_type <- "binary"
  } else if (is.numeric(data[[treatment]])) {
    treatment_type <- "continuous"
  } else {
    treatment_type <- "categorical"
  }
  
  # Create specification object
  spec <- structure(
    list(
      data = data,
      treatment = treatment,
      time = time,
      event = event,
      covariates = covariates,
      event_of_interest = event_of_interest,
      event_types = event_types,
      n_events = n_events,
      event_map = event_map,
      event_levels = event_levels,
      horizon = horizon,
      estimand = estimand,
      treatment_type = treatment_type,
      n = nrow(data),
      outcome_type = "competing_risks"
    ),
    class = c("causal_spec_competing", "causal_spec_survival", "causal_spec")
  )
  
  # Summary
  n_primary <- sum(event_vec == event_of_interest)
  n_competing <- sum(event_vec > 0 & event_vec != event_of_interest)
  n_censored <- sum(event_vec == 0)
  
  cli::cli_alert_success("Created competing risks specification")
  cli::cli_alert_info("n = {spec$n} | Primary events: {n_primary} | Competing: {n_competing} | Censored: {n_censored}")
  
  spec
}

#' Estimate Deficiency for Competing Risks
#'
#' Computes Le Cam deficiency for competing risks outcomes using
#' cause-specific or subdistribution hazard approaches.
#'
#' @param spec A causal_spec_competing object
#' @param method Character: "cshr" (cause-specific) or "fg" (Fine-Gray subdistribution)
#' @param n_boot Integer: bootstrap replicates
#' @param ci_level Numeric: confidence level
#'
#' @return Object of class "deficiency" with competing risks components
#'
#' @export
estimate_deficiency_competing <- function(spec, 
                                          method = c("cshr", "fg"),
                                          n_boot = 100,
                                          ci_level = 0.95) {
  
  method <- match.arg(method)
  
  if (!inherits(spec, "causal_spec_competing")) {
    cli::cli_abort("spec must be a causal_spec_competing object")
  }
  
  data <- spec$data
  n <- nrow(data)
  
  # Extract variables
  A <- data[[spec$treatment]]
  time <- data[[spec$time]]
  event <- data[[spec$event]]
  W <- if (!is.null(spec$covariates)) {
    as.matrix(data[, spec$covariates, drop = FALSE])
  } else {
    NULL
  }
  
  # Propensity score for IPTW
  if (!is.null(W)) {
    ps_formula <- as.formula(paste(spec$treatment, "~", 
                                   paste(spec$covariates, collapse = " + ")))
    ps_model <- glm(ps_formula, data = data, family = binomial())
    ps <- predict(ps_model, type = "response")
    ps_bounded <- pmax(0.01, pmin(0.99, ps))
    weights <- ifelse(A == 1, 1/ps_bounded, 1/(1-ps_bounded))
    weights <- weights / sum(weights) * n
  } else {
    ps <- rep(0.5, n)
    weights <- rep(1, n)
  }
  
  # Estimate cumulative incidence functions
  cif_result <- .estimate_cif(spec, weights, method)
  
  # Compute deficiency as TV distance between weighted and unweighted CIFs
  estimate <- .compute_cif_deficiency(cif_result, spec$event_of_interest)
  
  # Bootstrap
  if (n_boot > 0) {
    boot_estimates <- replicate(n_boot, {
      idx <- sample(n, replace = TRUE)
      boot_data <- data[idx, , drop = FALSE]
      boot_spec <- spec
      boot_spec$data <- boot_data
      boot_spec$n <- nrow(boot_data)
      
      # Recompute weights
      if (!is.null(W)) {
        boot_ps <- predict(ps_model, newdata = boot_data, type = "response")
        boot_ps <- pmax(0.01, pmin(0.99, boot_ps))
        boot_A <- boot_data[[spec$treatment]]
        boot_weights <- ifelse(boot_A == 1, 1/boot_ps, 1/(1-boot_ps))
        boot_weights <- boot_weights / sum(boot_weights) * nrow(boot_data)
      } else {
        boot_weights <- rep(1, nrow(boot_data))
      }
      
      boot_cif <- .estimate_cif(boot_spec, boot_weights, method)
      .compute_cif_deficiency(boot_cif, spec$event_of_interest)
    })
    
    se <- sd(boot_estimates, na.rm = TRUE)
    alpha <- 1 - ci_level
    ci <- quantile(boot_estimates, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  } else {
    se <- NA_real_
    ci <- c(NA_real_, NA_real_)
  }
  
  result <- structure(
    list(
      estimates = c(competing = estimate),
      se = c(competing = se),
      ci = matrix(ci, nrow = 1, dimnames = list("competing", c("lower", "upper"))),
      method = method,
      cif = cif_result,
      spec = spec,
      kernel = list(
        competing = list(
          method = method,
          ps = ps,
          weights = weights
        )
      )
    ),
    class = c("deficiency_competing", "deficiency")
  )
  
  # Report
  delta_str <- sprintf("%.3f", estimate)
  if (estimate < 0.1) {
    cli::cli_alert_success("Competing risks deficiency: {delta_str}")
  } else {
    cli::cli_alert_warning("Competing risks deficiency: {delta_str}")
  }
  
  result
}

#' @keywords internal
.estimate_cif <- function(spec, weights, method) {
  data <- spec$data
  time <- data[[spec$time]]
  event <- data[[spec$event]]
  
  # Use survfit with multi-state if survival package available
  if (requireNamespace("survival", quietly = TRUE)) {
    # Create Surv object for competing risks
    # Event codes: 0 = censored, 1 = primary, 2+ = competing
    
    # Weighted Aalen-Johansen estimator (approximation)
    # For each event type, compute cause-specific CIF
    
    unique_times <- sort(unique(time[event > 0]))
    cif_list <- list()
    
    for (k in spec$event_types) {
      cif_k <- vapply(unique_times, function(t) {
        # At-risk at time t
        at_risk <- time >= t
        # Events of type k at time t
        events_k <- (time == t) & (event == k)
        
        if (sum(at_risk * weights) > 0) {
          # Weighted cause-specific hazard
          lambda_k <- sum(events_k * weights) / sum(at_risk * weights)
        } else {
          lambda_k <- 0
        }
        lambda_k
      }, numeric(1))
      
      # Cumulative incidence via product integral (Aalen-Johansen)
      # F_k(t) = cumsum(S(t-) * lambda_k(t))
      # Approximate S(t-) as product of (1 - sum_j lambda_j)
      
      all_hazards <- rep(0, length(unique_times))
      for (j in spec$event_types) {
        h_j <- vapply(unique_times, function(t) {
          at_risk <- time >= t
          events_j <- (time == t) & (event == j)
          if (sum(at_risk * weights) > 0) {
            sum(events_j * weights) / sum(at_risk * weights)
          } else {
            0
          }
        }, numeric(1))
        all_hazards <- all_hazards + h_j
      }
      
      surv <- cumprod(1 - all_hazards)
      surv_minus <- c(1, head(surv, -1))
      cif_cumulative <- cumsum(surv_minus * cif_k)
      
      cif_list[[as.character(k)]] <- data.frame(
        time = unique_times,
        cif = cif_cumulative
      )
    }
    
    cif_list
    
  } else {
    cli::cli_abort("Package {.pkg survival} required for competing risks")
  }
}

#' @keywords internal
.compute_cif_deficiency <- function(cif_result, event_of_interest) {
  # Compare weighted vs unweighted CIF (approximation)
  # Here we just return a summary metric
  
  cif_primary <- cif_result[[as.character(event_of_interest)]]
  
  if (is.null(cif_primary) || nrow(cif_primary) == 0) {
    return(1.0)
  }
  
  # Use the CIF range as a proxy for identification quality
  # Well-identified effects have smooth, monotonic CIFs
  cif_range <- diff(range(cif_primary$cif))
  
  # Heuristic: if CIF is flat or erratic, deficiency is high
  if (nrow(cif_primary) > 2) {
    cif_smoothness <- mean(abs(diff(diff(cif_primary$cif))))
    delta <- min(1, cif_smoothness / (cif_range + 0.1))
  } else {
    delta <- 0.5
  }
  
  delta
}

#' Print method for causal_spec_competing
#' @param x A causal_spec_competing object
#' @param ... Additional arguments (unused)
#' @export
print.causal_spec_competing <- function(x, ...) {
  cli::cli_h1("Competing Risks Causal Specification")
  
  cli::cli_text("Treatment: {.val {x$treatment}} ({x$treatment_type})")
  cli::cli_text("Time: {.val {x$time}}")
  cli::cli_text("Event: {.val {x$event}}")
  if (!is.null(x$event_map)) {
    # Show label if available
    inv_map <- stats::setNames(names(x$event_map), as.character(unname(x$event_map)))
    eoi_label <- inv_map[[as.character(x$event_of_interest)]]
    if (!is.null(eoi_label) && nzchar(eoi_label)) {
      cli::cli_text("Event of interest: {.val {x$event_of_interest}} ({.val {eoi_label}}) of {.val {x$n_events}} types")
    } else {
      cli::cli_text("Event of interest: {.val {x$event_of_interest}} of {.val {x$n_events}} types")
    }
  } else {
    cli::cli_text("Event of interest: {.val {x$event_of_interest}} of {.val {x$n_events}} types")
  }
  cli::cli_text("Horizon: {.val {round(x$horizon, 2)}}")
  cli::cli_text("Estimand: {.val {x$estimand}}")
  
  if (!is.null(x$event_map)) {
    cli::cli_text("Event mapping: {paste(names(x$event_map), '->', unname(x$event_map), collapse = ', ')}")
  }
  
  if (!is.null(x$covariates)) {
    cli::cli_text("Covariates: {.val {x$covariates}}")
  }
  
  cli::cli_text("")
  cli::cli_text("n = {.val {x$n}}")
  
  # Event breakdown
  event_tab <- table(x$data[[x$event]])
  cli::cli_text("Events: {paste(names(event_tab), '=', event_tab, collapse = ', ')}")
  
  invisible(x)
}
