# =============================================================================
# policy_regret_bound() - Decision Theory Bounds
# =============================================================================

#' Compute Policy Regret Bounds
#'
#' Given deficiency \eqn{\delta}, computes worst-case bounds on policy regret.
#' The bound states: \eqn{Regret_{do} \le Regret_{obs} + 2M \cdot \delta}
#'
#' @param deficiency A deficiency object or numeric \eqn{\delta} value (between 0 and 1)
#' @param utility_range Numeric vector c(min, max) of utility bounds
#' @param obs_regret Numeric: observed regret from policy on observational data (optional)
#' @param policy_class Character: description of policy class (for reporting)
#'
#' @return Object of class "policy_bound" containing:
#'   \itemize{
#'     \item regret_bound: Upper bound on interventional regret (if obs_regret provided)
#'     \item safety_floor: Minimum unavoidable regret = 2*M*\eqn{\delta}
#'     \item delta: Deficiency value used
#'     \item utility_range: Range of utility function
#'     \item M: Utility range (max - min)
#'   }
#'
#' @details
#' For any policy \eqn{\pi} learned from observational data, the regret bound is:
#' \deqn{Regret_{do}(\pi) \leq Regret_{obs}(\pi) + 2M \cdot \delta}
#'
#' The "safety floor" \eqn{2M\delta} represents the irreducible regret penalty from
#' the information gap between observational and interventional experiments.
#' 
#' @references
#' Akdemir, D. (2026). Constraints on Causal Inference as Experiment Comparison.
#' DOI: 10.5281/zenodo.18367347. See Theorem 3.2 for proof.
#'
#' @section Implications for Safe AI:
#' If \eqn{\delta > 0} (unobserved confounding exists), no algorithm can guarantee
#' zero regret. The deficiency \eqn{\delta} is the "price of safety" for deploying
#' policies without randomized experimentation.
#'
#' @examples
#' # From a deficiency estimate
#' # From a deficiency estimate
#' df <- data.frame(W=rnorm(100), A=rbinom(100,1,0.5), Y=rnorm(100))
#' spec <- causal_spec(df, "A", "Y", "W")
#' def <- estimate_deficiency(spec, methods = "iptw", n_boot = 0)
#' bound <- policy_regret_bound(def, utility_range = c(0, 1))
#'
#' # From a numeric value
#' bound <- policy_regret_bound(
#'   deficiency = 0.1,
#'   utility_range = c(0, 100),
#'   obs_regret = 5
#' )
#' print(bound)
#'
#' @seealso [estimate_deficiency()], [confounding_frontier()]
#' @export
policy_regret_bound <- function(deficiency, utility_range = c(0, 1),
                                obs_regret = NULL, policy_class = NULL) {
  
  # Extract delta from deficiency object or use numeric value
  all_estimates <- NULL
  if (inherits(deficiency, "deficiency")) {
    delta <- unname(min(deficiency$estimates))  # Best method
    all_estimates <- deficiency$estimates
  } else {
    checkmate::assert_number(deficiency, lower = 0, upper = 1)
    delta <- deficiency
  }
  
  # Validate utility range
  checkmate::assert_numeric(utility_range, len = 2)
  if (utility_range[2] <= utility_range[1]) {
    .msg_error("utility_range must be specified as c(min, max) with min < max")
  }
  
  M <- diff(utility_range)
  
  # Theorem 3.2: safety floor = 2M*delta
  safety_floor <- 2 * M * delta
  
  # Full bound if obs_regret provided
  if (!is.null(obs_regret)) {
    checkmate::assert_number(obs_regret, lower = 0)
    regret_bound <- obs_regret + safety_floor
  } else {
    regret_bound <- NULL
  }
  
  result <- new_policy_bound(
    regret_bound = regret_bound,
    safety_floor = safety_floor,
    delta = delta,
    utility_range = utility_range,
    obs_regret = obs_regret,
    policy_class = policy_class,
    all_estimates = all_estimates
  )
  
  .msg_info(paste0("Safety floor: ", round(safety_floor, 4), 
                   " (minimum regret given delta = ", round(delta, 4), ")"))
  
  result
}
