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
#' The deficiency satisfies:
#' \deqn{\delta \geq C \frac{|\alpha\gamma|}{\sqrt{(1 + \alpha^2/\sigma_A^2)(1 + \gamma^2/\sigma_Y^2)}}}
#'
#' This creates a "confounding frontier" - the boundary where identification
#' transitions from possible (\eqn{\delta=0}) to impossible (\eqn{\delta>0}).
#'
#' @references
#' Akdemir, D. (2026). Constraints on Causal Inference as Experiment Comparison.
#' DOI: 10.5281/zenodo.18367347. See Theorem 4.1 for derivation.
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

  # Theorem 4.1 constant
  C <- 0.25

  # Compute deficiency for each (alpha, gamma) pair
  grid$delta <- mapply(function(a, g) {
    .deficiency_gaussian(a, g, sigma_A, sigma_Y, C)
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
          delta = .deficiency_gaussian(alpha_est, gamma_est, sigma_A, sigma_Y, C)
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
      C = C,
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
# Internal: Gaussian Deficiency Formula (Theorem 4.1)
# =============================================================================

#' Compute Deficiency for Linear Gaussian SCM
#'
#' Implements Theorem 4.1: \eqn{\delta \ge C|\alpha\gamma| / \sqrt{(1+\alpha^2/\sigma_A^2)(1+\gamma^2/\sigma_Y^2)}}
#'
#' @param alpha Numeric: U -> A confounding strength
#' @param gamma Numeric: U -> Y confounding strength
#' @param sigma_A Numeric: noise standard deviation for A
#' @param sigma_Y Numeric: noise standard deviation for Y
#' @param C Numeric: constant (default 0.25)
#'
#' @return Numeric: deficiency lower bound
#' @keywords internal
.deficiency_gaussian <- function(alpha, gamma, sigma_A, sigma_Y, C = 0.25) {
  num <- abs(alpha * gamma)
  denom <- sqrt(1 + alpha^2 / sigma_A^2) * sqrt(1 + gamma^2 / sigma_Y^2)
  C * num / denom
}
