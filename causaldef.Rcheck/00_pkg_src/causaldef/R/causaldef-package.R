#' @keywords internal
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
#' @importFrom data.table as.data.table data.table
## usethis namespace: end
NULL

# Global variables to suppress R CMD check NOTEs for data.table/ggplot2 NSE
utils::globalVariables(c(
    "covariate", "alpha", "gamma", "delta",
    "variable", "shift_metric", "severity", "weight"
))
