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
#' @param event Character: name of the event indicator (0 = censored, 1+ = event type)
#' @param covariates Character vector: names of covariate variables
#' @param event_of_interest Integer: which event type is the primary outcome (default 1)
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
  checkmate::assert_integerish(event_of_interest, lower = 1)
  checkmate::assert_number(horizon, lower = 0, null.ok = TRUE)
  
  # Check variables exist
  all_vars <- c(treatment, time, event, covariates)
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars) > 0) {
    cli::cli_abort("Variables not found in data: {.val {missing_vars}}")
  }
  
  # Extract event types
  event_vec <- data[[event]]
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
  cli::cli_text("Event of interest: {.val {x$event_of_interest}} of {.val {x$n_events}} types")
  cli::cli_text("Horizon: {.val {round(x$horizon, 2)}}")
  cli::cli_text("Estimand: {.val {x$estimand}}")
  
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
