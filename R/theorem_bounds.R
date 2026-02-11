# =============================================================================
# Theorem-Oriented Bounds & Diagnostics (Manuscript-aligned utilities)
# =============================================================================

#' RKHS Rate Bound (Result 1)
#'
#' Computes the explicit rate bound stated in Result 1 (finite-sample RKHS rate).
#'
#' @param n Sample size.
#' @param beta Smoothness exponent (\eqn{\beta} in a \eqn{\beta}-H\"older class).
#' @param d_w Covariate dimension (\eqn{d_W}).
#' @param eta Positivity constant (\eqn{\eta>0}) such that \eqn{P(A=a|W=w) \ge \eta}.
#' @param xi Failure probability (\eqn{\xi\in(0,1)}).
#'
#' @return Numeric rate bound (up to universal constants).
#' @export
rkhs_rate_bound <- function(n, beta, d_w, eta, xi = 0.05) {
  checkmate::assert_integerish(n, lower = 1)
  checkmate::assert_number(beta, lower = 0)
  checkmate::assert_integerish(d_w, lower = 1)
  checkmate::assert_number(eta, lower = 0, upper = 1)
  checkmate::assert_number(xi, lower = 0, upper = 1)
  if (eta == 0) {
    return(Inf)
  }
  if (xi == 0) {
    xi <- .Machine$double.eps
  }

  n <- as.numeric(n)
  d_w <- as.numeric(d_w)

  n^(-beta / (2 * beta + d_w)) + (1 / eta) * sqrt(log(1 / xi) / n)
}

#' Policy Regret Bound with VC Term (Result 3)
#'
#' Extends `policy_regret_bound()` with the explicit VC-complexity penalty from
#' Result 3: \eqn{C M \sqrt{(VC(\Pi)\log n + \log(1/\xi))/n}}.
#'
#' @param deficiency A `deficiency` object or numeric \eqn{\delta}.
#' @param vc_dim VC dimension of the policy class.
#' @param n Sample size used to learn the policy.
#' @param xi Failure probability (\eqn{\xi}).
#' @param utility_range Numeric vector c(min, max) of utility bounds.
#' @param obs_regret Optional observed regret under the observational regime.
#' @param delta_mode Passed through to `policy_regret_bound()` when `deficiency` is an object.
#' @param C Universal constant in the bound (default 2).
#'
#' @return An object of class `policy_bound` with additional field `complexity_penalty`.
#' @export
policy_regret_bound_vc <- function(deficiency,
                                   vc_dim,
                                   n,
                                   xi = 0.05,
                                   utility_range = c(0, 1),
                                   obs_regret = NULL,
                                   delta_mode = c("point", "upper"),
                                   C = 2) {
  checkmate::assert_number(vc_dim, lower = 0)
  checkmate::assert_integerish(n, lower = 1)
  checkmate::assert_number(xi, lower = 0, upper = 1)
  checkmate::assert_number(C, lower = 0)

  base <- policy_regret_bound(
    deficiency = deficiency,
    utility_range = utility_range,
    obs_regret = obs_regret,
    delta_mode = delta_mode
  )

  if (xi == 0) {
    xi <- .Machine$double.eps
  }

  n <- as.numeric(n)
  vc_dim <- as.numeric(vc_dim)

  complexity_penalty <- C * base$M * sqrt((vc_dim * log(n) + log(1 / xi)) / n)
  base$complexity_penalty <- complexity_penalty

  if (!is.null(base$regret_bound)) {
    base$regret_bound <- base$regret_bound + complexity_penalty
  }

  class(base) <- unique(c("policy_bound_vc", class(base)))
  base
}

#' Wasserstein Deficiency (Linear Gaussian) (Result 6)
#'
#' Closed-form Wasserstein-1 deficiency for the two-point construction in Result 6:
#' \eqn{\delta_{W_1}(a) = |a \alpha \gamma| / (\alpha^2 + \sigma_A^2)}.
#'
#' @param alpha Confounding strength U -> A.
#' @param gamma Confounding strength U -> Y.
#' @param sigma_A Treatment noise standard deviation.
#' @param a Intervention level \eqn{do(A=a)}.
#'
#' @return Numeric deficiency value.
#' @export
wasserstein_deficiency_gaussian <- function(alpha, gamma, sigma_A = 1, a = 1) {
  checkmate::assert_number(alpha)
  checkmate::assert_number(gamma)
  checkmate::assert_number(sigma_A, lower = 0)
  checkmate::assert_number(a)
  denom <- alpha^2 + sigma_A^2
  if (!is.finite(denom) || denom <= 0) {
    return(0)
  }
  abs(a * alpha * gamma) / denom
}

#' Sharp Two-Point Bounds (Result 4)
#'
#' Computes matching lower/upper bounds (sharpness) in the linear Gaussian model,
#' for either TV-based (two-point construction underlying `thm:confounding_lb`) or Wasserstein-1-based
#' (Result 6 two-point construction) deficiencies.
#'
#' @param alpha Confounding strength U -> A.
#' @param gamma Confounding strength U -> Y.
#' @param sigma_A Treatment noise standard deviation.
#' @param sigma_Y Outcome noise standard deviation (TV case).
#' @param a Intervention level.
#' @param metric One of `"tv"` or `"wasserstein"`.
#'
#' @return List with fields: `lower`, `upper`, `ratio`, and `metric`.
#' @export
sharp_lower_bound <- function(alpha,
                             gamma,
                             sigma_A = 1,
                             sigma_Y = 1,
                             a = 1,
                             metric = c("tv", "wasserstein")) {
  metric <- match.arg(metric)
  checkmate::assert_number(alpha)
  checkmate::assert_number(gamma)
  checkmate::assert_number(sigma_A, lower = 0)
  checkmate::assert_number(a)

  if (metric == "wasserstein") {
    delta <- wasserstein_deficiency_gaussian(alpha = alpha, gamma = gamma, sigma_A = sigma_A, a = a)
    return(list(
      metric = "wasserstein",
      lower = delta,
      upper = delta,
      ratio = if (delta == 0) 1 else 1
    ))
  }

  checkmate::assert_number(sigma_Y, lower = 0)

  # Reuse the package's exact two-point TV computation (confounding_frontier).
  # .deficiency_gaussian returns (1/2)*TV between two do-laws in the Le Cam pair,
  # which is both a lower bound and (by midpoint kernel) an upper bound for the
  # two-point subexperiment.
  delta <- .deficiency_gaussian(alpha = alpha, gamma = gamma, sigma_A = sigma_A, sigma_Y = sigma_Y)

  list(
    metric = "tv",
    lower = delta,
    upper = delta,
    ratio = if (delta == 0) 1 else 1
  )
}

#' Partial Identification Set from a Delta Radius (Result 8)
#'
#' Converts a deficiency radius \eqn{\delta} into a conservative identified interval
#' for bounded functionals under TV.
#'
#' For a mean functional with \eqn{Y \in [y_{min}, y_{max}]}, TV control yields
#' \eqn{|E_Q[Y] - E_{Q_0}[Y]| \le 2\delta (y_{max}-y_{min})}.
#' For an ATE (difference of two means), the half-width doubles again.
#'
#' @param estimate Point estimate of the functional (numeric).
#' @param delta Deficiency radius \eqn{\delta \in [0,1]}.
#' @param estimand One of `"mean"` or `"ate"`.
#' @param outcome_range Numeric vector c(min, max) bounding the outcome.
#'
#' @return List with `lower`, `upper`, `estimate`, and `half_width`.
#' @export
partial_id_set <- function(estimate,
                           delta,
                           estimand = c("ate", "mean"),
                           outcome_range = c(0, 1)) {
  estimand <- match.arg(estimand)
  checkmate::assert_number(estimate)
  checkmate::assert_number(delta, lower = 0, upper = 1)
  checkmate::assert_numeric(outcome_range, len = 2)

  y_min <- outcome_range[1]
  y_max <- outcome_range[2]
  if (!is.finite(y_min) || !is.finite(y_max) || y_max <= y_min) {
    cli::cli_abort("outcome_range must be c(min, max) with min < max.")
  }
  span <- y_max - y_min

  factor <- if (estimand == "mean") 2 else 4
  half_width <- factor * delta * span

  out <- list(
    estimand = estimand,
    estimate = estimate,
    delta = delta,
    outcome_range = outcome_range,
    half_width = half_width,
    lower = estimate - half_width,
    upper = estimate + half_width
  )

  class(out) <- "partial_id_set"
  out
}

#' Overlap Diagnostic
#'
#' Summarizes overlap/positivity issues via propensity score distribution and
#' effective sample size under IPTW weights.
#'
#' @param spec A `causal_spec`.
#' @param trim Trimming threshold in (0, 0.5). Values outside \eqn{[trim, 1 - trim]}
#'   are flagged as extreme.
#'
#' @return List with propensity summary and effective sample size.
#' @export
overlap_diagnostic <- function(spec, trim = 0.01) {
  validate_causal_spec(spec)
  checkmate::assert_number(trim, lower = 0, upper = 0.49)

  if (is.null(spec$covariates) || length(spec$covariates) == 0) {
    cli::cli_abort("overlap_diagnostic() requires covariates in the causal_spec.")
  }
  if (spec$treatment_type != "binary") {
    cli::cli_abort("overlap_diagnostic() currently supports only binary treatments.")
  }

  data <- spec$data
  A <- data[[spec$treatment]]
  if (is.factor(A) || is.character(A)) {
    A <- as.numeric(as.factor(A)) - 1
  }

  ps_formula <- stats::as.formula(paste(spec$treatment, "~", paste(spec$covariates, collapse = " + ")))
  ps_model <- stats::glm(ps_formula, data = data, family = stats::binomial())
  ps <- stats::predict(ps_model, type = "response")

  ps_bounded <- pmax(1e-6, pmin(1 - 1e-6, ps))
  w <- A / ps_bounded + (1 - A) / (1 - ps_bounded)

  ess <- (sum(w)^2) / sum(w^2)
  extreme <- sum(ps < trim | ps > (1 - trim))
  kept <- sum(ps >= trim & ps <= (1 - trim))

  out <- list(
    trim = trim,
    n = nrow(data),
    extreme_n = extreme,
    kept_n = kept,
    ps_summary = stats::quantile(ps, probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1)),
    ess_iptw = ess
  )
  class(out) <- "overlap_diagnostic"
  out
}
