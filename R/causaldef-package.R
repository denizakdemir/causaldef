#' causaldef: Decision-Theoretic Causal Diagnostics via Le Cam Deficiency
#'
#' Theory-forward causal diagnostics organized around Le Cam deficiency.
#'
#' The package centers causal inference around the question: how much information
#' separates the observational study at hand from the interventional experiment we
#' would ideally like to run?
#'
#' \strong{Scientific contract}
#'
#' The package exposes four kinds of quantities:
#' \itemize{
#'   \item \strong{Theorem-backed utilities}: closed-form or theorem-aligned bounds such
#'   as [policy_regret_bound()], [policy_regret_bound_vc()], [confounding_frontier()],
#'   [sharp_lower_bound()], and [wasserstein_deficiency_gaussian()].
#'   \item \strong{Computable deficiency proxies}: [estimate_deficiency()] currently returns
#'   a propensity-score total-variation proxy (`metric = "ps_tv"`), not a generic
#'   nonparametric estimator of the exact Le Cam deficiency.
#'   \item \strong{Sensitivity diagnostics}: functions such as [nc_diagnostic()] combine
#'   observable diagnostics with user-supplied sensitivity parameters.
#'   \item \strong{Experimental heuristics}: some modules currently expose effect estimates
#'   together with heuristic risk scores or proxy summaries. These should be read as
#'   diagnostics unless the individual help page states a theorem-backed identification
#'   result.
#' }
#'
#' \strong{Recommended workflow}
#'
#' \enumerate{
#'   \item Specify the causal problem with [causal_spec()] or a survival variant.
#'   \item Compute a pre-specified diagnostic or proxy with [estimate_deficiency()] or
#'   a specialized module.
#'   \item Stress-test assumptions with sensitivity tools such as [nc_diagnostic()] or
#'   [confounding_frontier()].
#'   \item Translate the resulting information gap into decision consequences with
#'   [policy_regret_bound()].
#' }
#'
#' @references
#' Akdemir, D. (2026). Constraints on Causal Inference as Experiment Comparison: A Framework for Identification, Transportability, and Policy Learning. DOI: 10.5281/zenodo.18367347
#'
#' Le Cam, L., & Yang, G. L. (2000). Asymptotics in Statistics: Some Basic Concepts. Springer.
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats binomial coef confint complete.cases glm lm na.omit plogis
#'   pnorm predict qnorm quantile rbinom rnorm sd var weighted.mean as.formula
#'   chisq.test cov fitted model.matrix pf reorder setNames t.test
#' @importFrom graphics par plot abline legend text
#' @importFrom utils capture.output head
#' @importFrom checkmate assert_class assert_character assert_choice
#'   assert_data_frame assert_flag assert_integerish assert_number
#'   assert_numeric assert_string assert_subset
#' @importFrom cli cli_abort cli_alert_danger cli_alert_info cli_alert_success
#'   cli_alert_warning cli_progress_bar cli_progress_done cli_progress_update
#'   cli_warn
## usethis namespace: end
NULL

# Global variables to suppress R CMD check NOTEs for ggplot2 NSE
utils::globalVariables(c(
    "covariate", "alpha", "gamma", "delta",
    "variable", "shift_metric", "severity", "weight"
))
