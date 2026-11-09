# =============================================================================
# policy_regret_bound() - Decision Theory Bounds
# =============================================================================

#' Compute Policy Regret Bounds
#'
#' Given a deficiency \eqn{\delta}, computes the regret transfer penalty and the
#' minimax regret floor of the manuscript's Policy Regret Transfer and Minimax
#' Safety Floor theorems.
#'
#' @param deficiency A deficiency object or numeric \eqn{\delta} value (between 0 and 1)
#' @param utility_range Numeric vector c(min, max) of utility bounds
#' @param obs_regret Numeric: regret of the policy computed under the simulated
#'   (adjusted) interventional laws, i.e. the regret the analyst can evaluate from
#'   observational data after adjustment (optional)
#' @param policy_class Character: description of policy class (for reporting)
#' @param delta_mode Character: how to extract \eqn{\delta} from a deficiency object.
#'   \itemize{
#'     \item "point": use the minimum point estimate (default)
#'     \item "upper": use the minimum upper CI bound (more conservative; requires bootstrap)
#'   }
#' @param method Character: when `deficiency` is a `deficiency` object, which fitted
#'   method to use when translating the estimate into a regret bound. Strongly
#'   recommended whenever multiple methods were compared.
#' @param n_actions Integer: number of actions \eqn{|\mathcal{A}|} (default 2, a
#'   binary treatment).
#' @param policy_type Character: `"covariate"` (default) for policies that choose
#'   the action as a function of covariates, or `"constant"` for comparing fixed
#'   actions (e.g. treat-all vs treat-none), where the factor \eqn{|\mathcal{A}|}
#'   is not needed.
#'
#' @return Object of class "policy_bound" containing:
#'   \itemize{
#'     \item regret_bound: Upper bound on interventional regret (if obs_regret provided)
#'     \item value_error: Bound on the error of any policy's simulated value,
#'       \eqn{|\mathcal{A}| M\delta} (covariate policies) or \eqn{M\delta} (constant policies)
#'     \item transfer_penalty: Additive regret penalty \eqn{= 2 \times} value_error
#'     \item safety_floor: Backward-compatible alias for transfer_penalty
#'     \item minimax_floor: Minimax lower bound \eqn{= (M/4)\delta}
#'     \item transfer_multiplier, floor_multiplier: the factors multiplying
#'       \eqn{M\delta} in transfer_penalty and minimax_floor
#'     \item n_actions, policy_type: as supplied
#'     \item delta: Deficiency value used
#'     \item delta_selection: How the \eqn{\delta} value was selected
#'     \item utility_range: Range of utility function
#'     \item M: Utility range (max - min)
#'   }
#'
#' @details
#' Let \eqn{\delta} be the worst-case total-variation error of the kernels used to
#' simulate each interventional law. For any policy \eqn{\pi}, the simulated value
#' differs from the interventional value by at most \eqn{|\mathcal{A}| M\delta}
#' (\eqn{M\delta} for constant policies), and therefore
#' \deqn{Regret_{do}(\pi) \leq Regret_{sim}(\pi) + 2|\mathcal{A}| M \delta.}
#' The factor 2 cannot be removed in general. Conversely, for utilities in a range
#' of width \eqn{M}, every observational learner has worst-case regret at least
#' \eqn{(M/4)\delta} for some bounded utility (the "minimax floor"; the constant
#' 1/4 is sharp).
#'
#' Versions up to 0.2.1 reported a transfer penalty of \eqn{M\delta} and a minimax
#' floor of \eqn{(M/2)\delta}. Both constants were wrong and were corrected in 0.3.0
#' (see NEWS).
#'
#' When `deficiency` comes from [estimate_deficiency()] it is a propensity-score
#' TV proxy, not the exact deficiency, so the result is a plug-in heuristic.
#'
#' When a multi-method `deficiency` object is supplied, the safest workflow is to
#' pre-specify `method` before calling `policy_regret_bound()`. If `method` is omitted,
#' the current implementation falls back to the smallest available estimate (or upper
#' confidence limit), which is convenient for exploration but optimistic after
#' post-selection.
#' 
#' @references
#' Akdemir, D. (2026). Constraints on Causal Inference as Experiment Comparison.
#' DOI: 10.5281/zenodo.18367346. See `thm:policy_regret` (Policy Regret Transfer)
#' and `thm:safety_floor` (Minimax Safety Floor).
#'
#' @section Implications for Safe AI:
#' If \eqn{\delta > 0} (unobserved confounding exists), no algorithm can guarantee
#' zero regret. The deficiency \eqn{\delta} is the "price of safety" for deploying
#' policies without randomized experimentation.
#'
#' @examples
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
                                obs_regret = NULL, policy_class = NULL,
                                delta_mode = c("point", "upper"),
                                method = NULL, n_actions = 2,
                                policy_type = c("covariate", "constant")) {
  
  delta_mode <- match.arg(delta_mode)
  policy_type <- match.arg(policy_type)
  checkmate::assert_integerish(n_actions, lower = 2, len = 1)
  
  # Extract delta from deficiency object or use numeric value
  all_estimates <- NULL
  delta_method <- NULL
  delta_selection <- "numeric_input"
  delta_is_proxy <- FALSE
  if (inherits(deficiency, "deficiency")) {
    all_estimates <- deficiency$estimates
    delta_is_proxy <- identical(deficiency$metric, "ps_tv")
    available_methods <- names(deficiency$estimates)

    if (!is.null(method)) {
      checkmate::assert_choice(method, available_methods)
      candidate_methods <- method
      delta_selection <- "explicit_method"
    } else if (length(available_methods) == 1) {
      candidate_methods <- available_methods
      delta_selection <- "single_method"
    } else {
      candidate_methods <- available_methods
      delta_selection <- "optimistic_minimum"
      cli::cli_warn(c(
        "Multiple fitted methods are available but {.arg method} was not specified.",
        "i" = "Using the smallest available delta across methods is optimistic after model selection.",
        "i" = "For a pre-specified decision bound, call {.fn policy_regret_bound} with {.code method = '<chosen method>'}."
      ))
    }

    if (delta_mode == "upper") {
      if (!is.null(deficiency$ci) && !all(is.na(deficiency$ci)) && ncol(deficiency$ci) == 2) {
        upper <- deficiency$ci[candidate_methods, 2]
        if (is.null(names(upper))) {
          names(upper) <- candidate_methods
        }
        delta_method <- names(which.min(upper))
        delta <- unname(min(upper))
      } else {
        .msg_warning("No confidence intervals available; falling back to point estimates for delta.")
        subset_estimates <- deficiency$estimates[candidate_methods]
        delta_method <- names(which.min(subset_estimates))
        delta <- unname(min(subset_estimates))
      }
    } else {
      subset_estimates <- deficiency$estimates[candidate_methods]
      delta_method <- names(which.min(subset_estimates))
      delta <- unname(min(subset_estimates))
    }
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
  
  # Value error and regret transfer penalty (Policy Regret Transfer theorem):
  # |V_sim - V_do| <= k * M * delta with k = n_actions (covariate policies) or 1
  # (constant policies); regret penalty = 2 * value error.
  value_multiplier <- if (policy_type == "covariate") as.numeric(n_actions) else 1
  value_error <- value_multiplier * M * delta
  transfer_multiplier <- 2 * value_multiplier
  transfer_penalty <- transfer_multiplier * M * delta

  # Minimax floor (Minimax Safety Floor theorem; constant 1/4 is sharp)
  floor_multiplier <- 0.25
  minimax_floor <- floor_multiplier * M * delta

  # Backward-compatible alias
  safety_floor <- transfer_penalty
  
  # Full bound if obs_regret provided
  if (!is.null(obs_regret)) {
    checkmate::assert_number(obs_regret, lower = 0)
    regret_bound <- obs_regret + transfer_penalty
  } else {
    regret_bound <- NULL
  }
  
  result <- new_policy_bound(
    regret_bound = regret_bound,
    safety_floor = safety_floor,
    transfer_penalty = transfer_penalty,
    minimax_floor = minimax_floor,
    value_error = value_error,
    transfer_multiplier = transfer_multiplier,
    floor_multiplier = floor_multiplier,
    n_actions = as.integer(n_actions),
    policy_type = policy_type,
    delta = delta,
    utility_range = utility_range,
    obs_regret = obs_regret,
    policy_class = policy_class,
    all_estimates = all_estimates,
    delta_mode = delta_mode,
    delta_method = delta_method,
    delta_selection = delta_selection,
    delta_is_proxy = delta_is_proxy
  )
  
  .msg_info(paste0("Transfer penalty: ", round(transfer_penalty, 4), 
                   " (delta = ", round(delta, 4), ")"))
  
  result
}
