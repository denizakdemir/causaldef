# =============================================================================
# confounding_frontier() - Sensitivity Analysis
# =============================================================================

#' Map the Confounding Frontier
#'
#' Computes deficiency as a function of confounding strengths (\eqn{\alpha}, \eqn{\gamma})
#' for linear Gaussian models.
#'
#' @param spec A causal_spec object (optional, for parameter estimation)
#' @param alpha_range Numeric vector c(min, max): range of U->A confounding strength
#' @param gamma_range Numeric vector c(min, max): range of U->Y confounding strength
#' @param grid_size Integer: resolution per dimension
#' @param model Character: model type ("gaussian" for closed-form)
#'
#' @return Object of class "confounding_frontier" containing:
#'   \itemize{
#'     \item grid: Data frame with columns alpha, gamma, delta
#'     \item frontier: Subset of grid where delta \eqn{\approx} 0
#'     \item model: Model type used
#'     \item params: Parameters used in computation
#'   }
#'
#' @details
#' For the linear Gaussian structural causal model:
#' \deqn{U \sim N(0,1), \quad A = \alpha U + \varepsilon_A, \quad Y = \beta A + \gamma U + \varepsilon_Y}
#'
#' The manuscript's Confounding Lower Bound theorem (`thm:confounding_lb`) states the scaling of the causal deficiency:
#' \deqn{\delta \geq C \frac{|\alpha\gamma|}{\sqrt{(1 + \alpha^2/\sigma_A^2)(1 + \gamma^2/\sigma_Y^2)}}}
#' for some universal constant \eqn{C>0} (typically after rescaling so \eqn{|a|\le 1}).
#' This function computes a concrete, rigorous lower bound using the underlying
#' Le Cam two-point construction (i.e., \eqn{\delta \ge (1/2)\,\mathrm{TV}} between two
#' observationally indistinguishable parameterizations).
#'
#' This creates a "confounding frontier" - the boundary where identification
#' transitions from possible (\eqn{\delta=0}) to impossible (\eqn{\delta>0}).
#'
#' @references
#' Akdemir, D. (2026). Constraints on Causal Inference as Experiment Comparison.
#' DOI: 10.5281/zenodo.18367347. See `thm:confounding_lb` (Confounding Lower Bound).
#'
#' @section Interpretation:
#' \describe{
#'   \item{\eqn{\alpha}}{Strength of confounding path U -> A}
#'   \item{\eqn{\gamma}}{Strength of confounding path U -> Y}
#'   \item{\eqn{|\alpha\gamma|}}{Product determines the confounding bias}
#' }
#'
#' When alpha=0 OR gamma=0, the confounder has no effect and \eqn{\delta=0}.
#'
#' @examples
#' # Basic frontier
#' frontier <- confounding_frontier(
#'   alpha_range = c(-2, 2),
#'   gamma_range = c(-2, 2),
#'   grid_size = 25
#' )
#'
#' # With data-based parameter estimation
#' df <- data.frame(A = rnorm(100), Y = rnorm(100), W = rnorm(100))
#' spec <- causal_spec(df, "A", "Y", "W")
#' frontier <- confounding_frontier(spec, grid_size = 50)
#'
#' @seealso [estimate_deficiency()], [policy_regret_bound()]
#' @export
confounding_frontier <- function(spec = NULL,
                                 alpha_range = c(-2, 2),
                                 gamma_range = c(-2, 2),
                                 grid_size = 25,
                                 model = "gaussian") {
  checkmate::assert_numeric(alpha_range, len = 2)
  checkmate::assert_numeric(gamma_range, len = 2)
  checkmate::assert_integerish(grid_size, lower = 5)

  # Estimate noise parameters from data if spec provided
  if (!is.null(spec) && inherits(spec, "causal_spec")) {
    # Handle factor treatment (convert to numeric 0/1)
    treatment_var <- spec$data[[spec$treatment]]
    if (is.factor(treatment_var) || is.character(treatment_var)) {
      treatment_var <- as.numeric(as.factor(treatment_var)) - 1
    }
    sigma_A <- sd(treatment_var, na.rm = TRUE)

    # Handle outcome (for survival, use time; otherwise use outcome)
    if (inherits(spec, "causal_spec_survival")) {
      outcome_var <- spec$data[[spec$time]]
    } else {
      outcome_var <- spec$data[[spec$outcome]]
    }
    if (is.factor(outcome_var) || is.character(outcome_var)) {
      outcome_var <- as.numeric(as.factor(outcome_var)) - 1
    }
    sigma_Y <- sd(outcome_var, na.rm = TRUE)
  } else {
    sigma_A <- 1
    sigma_Y <- 1
  }

  # Create grid
  alpha_seq <- seq(alpha_range[1], alpha_range[2], length.out = grid_size)
  gamma_seq <- seq(gamma_range[1], gamma_range[2], length.out = grid_size)
  grid <- expand.grid(alpha = alpha_seq, gamma = gamma_seq)

  # Compute deficiency for each (alpha, gamma) pair
  grid$delta <- mapply(function(a, g) {
    .deficiency_gaussian(a, g, sigma_A, sigma_Y)
  }, grid$alpha, grid$gamma)

  # Identify frontier (delta approx 0)
  frontier <- grid[grid$delta < 0.01, ]

  # Compute Benchmarks (observed covariates)
  benchmarks <- NULL
  if (!is.null(spec) && !is.null(spec$covariates) && length(spec$covariates) > 0) {
    .msg_info("Computing benchmarks for observed covariates...")

    data <- spec$data
    A_name <- spec$treatment
    Y_name <- spec$outcome # or spec$time for survival, but linear approx usually uses outcome

    # Handle survival outcome for benchmarking by using event/time naïve approximation
    # or just skipping if not linear-amenable.
    # For now, we only benchmark if outcome is continuous/binary numeric
    if (!inherits(spec, "causal_spec_survival")) {
      benchmarks <- do.call(rbind, lapply(spec$covariates, function(w) {
        # Standardize W to mimic U ~ N(0,1)
        W_vec <- data[[w]]
        if (!is.numeric(W_vec)) {
          return(NULL)
        } # Skip factors for simple benchmarking

        W_std <- scale(W_vec)

        # alpha: A ~ W
        # A = alpha*U + ...
        m_a <- lm(data[[A_name]] ~ W_std)
        alpha_est <- coef(m_a)[2]

        # gamma: Y ~ A + W
        # Y = beta*A + gamma*U + ...
        m_y <- lm(data[[Y_name]] ~ data[[A_name]] + W_std)
        gamma_est <- coef(m_y)[3]

        data.frame(
          covariate = w,
          alpha = alpha_est,
          gamma = gamma_est,
          delta = .deficiency_gaussian(alpha_est, gamma_est, sigma_A, sigma_Y)
        )
      }))
    }
  }

  result <- new_confounding_frontier(
    grid = grid,
    frontier = frontier,
    model = model,
    params = list(
      alpha_range = alpha_range,
      gamma_range = gamma_range,
      sigma_A = sigma_A,
      sigma_Y = sigma_Y,
      grid_size = grid_size
    )
  )

  # Attach benchmarks if available
  if (!is.null(benchmarks)) {
    result$benchmarks <- benchmarks
  }

  .msg_success(paste0("Computed confounding frontier: ", grid_size, "x", grid_size, " grid"))

  result
}

# =============================================================================
# Internal: Gaussian Deficiency Formula (Confounding Lower Bound)
# =============================================================================

.tv_distance_normal <- function(mu1, sd1, mu2, sd2) {
  if (!is.finite(sd1) || !is.finite(sd2) || sd1 <= 0 || sd2 <= 0) {
    return(NA_real_)
  }
  if (identical(mu1, mu2) && identical(sd1, sd2)) {
    return(0)
  }

  # Solve for cutpoints where the densities intersect.
  a_q <- 1 / sd2^2 - 1 / sd1^2
  b_q <- -2 * mu2 / sd2^2 + 2 * mu1 / sd1^2
  c_q <- mu2^2 / sd2^2 - mu1^2 / sd1^2 + 2 * log(sd2 / sd1)

  if (abs(a_q) < 1e-12) {
    x <- -c_q / b_q
    cutpoints <- x
  } else {
    disc <- b_q^2 - 4 * a_q * c_q
    disc <- max(0, disc)
    sqrt_disc <- sqrt(disc)
    x1 <- (-b_q - sqrt_disc) / (2 * a_q)
    x2 <- (-b_q + sqrt_disc) / (2 * a_q)
    cutpoints <- sort(c(x1, x2))
  }

  endpoints <- c(-Inf, cutpoints, Inf)
  tv <- 0
  scale <- 10 * max(sd1, sd2)

  for (i in seq_len(length(endpoints) - 1)) {
    lo <- endpoints[i]
    hi <- endpoints[i + 1]

    if (is.infinite(lo) && lo < 0) {
      mid <- hi - scale
    } else if (is.infinite(hi) && hi > 0) {
      mid <- lo + scale
    } else {
      mid <- 0.5 * (lo + hi)
    }

    f1 <- stats::dnorm(mid, mean = mu1, sd = sd1)
    f2 <- stats::dnorm(mid, mean = mu2, sd = sd2)

    if (f1 >= f2) {
      p1 <- stats::pnorm(hi, mean = mu1, sd = sd1) - stats::pnorm(lo, mean = mu1, sd = sd1)
      p2 <- stats::pnorm(hi, mean = mu2, sd = sd2) - stats::pnorm(lo, mean = mu2, sd = sd2)
      tv <- tv + (p1 - p2)
    }
  }

  min(1, max(0, tv))
}

#' Compute Deficiency for Linear Gaussian SCM
#'
#' Computes a rigorous Le Cam two-point lower bound for the confounding lower bound
#' theorem in the manuscript (`thm:confounding_lb`) in the linear Gaussian setting
#' (with an implicit normalization \eqn{|a|\le 1}).
#'
#' @param alpha Numeric: U -> A confounding strength
#' @param gamma Numeric: U -> Y confounding strength
#' @param sigma_A Numeric: noise standard deviation for A
#' @param sigma_Y Numeric: noise standard deviation for Y
#'
#' @return Numeric: deficiency lower bound
#' @keywords internal
.deficiency_gaussian <- function(alpha, gamma, sigma_A, sigma_Y) {
  denom_A <- alpha^2 + sigma_A^2
  if (!is.finite(denom_A) || denom_A <= 0) {
    return(0)
  }

  # Two-point construction yields mean shift |Δμ| = |a|·|αγ|/(α²+σ_A²).
  # We take |a|=1 by default (consistent with the manuscript's normalization).
  delta_mu <- abs(alpha * gamma) / denom_A

  # With Var(Y) matched observationally, γ'² = γ²·σ_A²/(α²+σ_A²).
  gamma_prime_sq <- (gamma^2) * (sigma_A^2) / denom_A

  sd1 <- sqrt(gamma^2 + sigma_Y^2)
  sd2 <- sqrt(gamma_prime_sq + sigma_Y^2)

  tv <- .tv_distance_normal(0, sd1, delta_mu, sd2)
  if (!is.finite(tv)) {
    return(NA_real_)
  }

  # Le Cam two-point method: δ ≥ (1/2)·TV.
  0.5 * tv
}
