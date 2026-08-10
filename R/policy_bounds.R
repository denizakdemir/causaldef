# =============================================================================
# policy_regret_bound() - Decision Theory Bounds
# =============================================================================

#' Compute Policy Regret Bounds
#'
#' Given deficiency \eqn{\delta}, computes worst-case bounds on policy regret.
#' The bound states: \eqn{Regret_{do} \le Regret_{obs} + M \cdot \delta}
#'
#' @param deficiency A deficiency object or numeric \eqn{\delta} value (between 0 and 1)
#' @param utility_range Numeric vector c(min, max) of utility bounds
#' @param obs_regret Numeric: observed regret from policy on observational data (optional)
#' @param policy_class Character: description of policy class (for reporting)
#' @param delta_mode Character: how to extract \eqn{\delta} from a deficiency object.
#'   \itemize{
#'     \item "point": use the minimum point estimate (default)
#'     \item "upper": use the minimum upper CI bound (more conservative; requires bootstrap)
#'   }
#' @param method Character: when `deficiency` is a `deficiency` object, which fitted
#'   method to use when translating the estimate into a regret bound. Strongly
#'   recommended whenever multiple methods were compared.
#'
#' @return Object of class "policy_bound" containing:
#'   \itemize{
#'     \item regret_bound: Upper bound on interventional regret (if obs_regret provided)
#'     \item safety_floor: Backward-compatible alias for transfer_penalty
#'     \item transfer_penalty: Additive regret penalty = M*\eqn{\delta}
#'     \item minimax_floor: Minimax lower bound = (M/2)*\eqn{\delta}
#'     \item delta: Deficiency value used
#'     \item delta_selection: How the \eqn{\delta} value was selected
#'     \item utility_range: Range of utility function
#'     \item M: Utility range (max - min)
#'   }
#'
#' @details
#' For any policy \eqn{\pi} learned from observational data, the regret bound is:
#' \deqn{Regret_{do}(\pi) \leq Regret_{obs}(\pi) + M \cdot \delta}
#'
#' The "transfer penalty" \eqn{M\delta} quantifies worst-case regret inflation from
#' the information gap between observational and interventional experiments.
#'
#' When a multi-method `deficiency` object is supplied, the safest workflow is to
#' pre-specify `method` before calling `policy_regret_bound()`. If `method` is omitted,
#' the current implementation falls back to the smallest available estimate (or upper
#' confidence limit), which is convenient for exploration but optimistic after
#' post-selection.
#'
#' Separately, a minimax lower bound ("minimax floor") of \eqn{(M/2)\delta} holds for
#' bounded utilities (see manuscript theorem `thm:safety_floor`).
#' 
#' @references
#' Akdemir, D. (2026). Constraints on Causal Inference as Experiment Comparison.
#' DOI: 10.5281/zenodo.21877511. See `thm:policy_regret` (Policy Regret Transfer)
#' and `thm:safety_floor` (Minimax Safety Floor).
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
                                obs_regret = NULL, policy_class = NULL,
                                delta_mode = c("point", "upper"),
                                method = NULL) {
  
  delta_mode <- match.arg(delta_mode)
  
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
  
  # Regret transfer penalty (upper bound additive term)
  transfer_penalty <- M * delta

  # Minimax floor (lower bound; see manuscript)
  minimax_floor <- 0.5 * M * delta

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
