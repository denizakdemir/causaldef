# =============================================================================
# audit_data() - Automated Data Auditor
# =============================================================================

#' Audit Data for Causal Validity
#'
#' Automatically scans a dataset for causal validity issues using
#' negative control diagnostics (manuscript `thm:nc_bound`). Identifies potential
#' confounders, instruments, and negative control variables.
#'
#' @param data A data.frame containing the analysis data
#' @param treatment Character: name of the treatment variable
#' @param outcome Character: name of the outcome variable
#' @param covariates Character vector: specific covariates to audit (NULL = all non-treatment/outcome columns)
#' @param negative_controls Character vector: known negative control outcomes to validate
#' @param alpha Numeric: significance level for tests (default 0.05)
#' @param verbose Logical: print progress messages
#'
#' @return Object of class "data_audit_report" containing:
#'   \itemize{
#'     \item issues: data.frame with columns (variable, issue_type, p_value, recommendation)
#'     \item recommendations: character vector of action items
#'     \item summary_stats: list with n_vars_audited, n_issues, etc.
#'   }
#'
#' @details
#' The auditor tests each variable for:
#' \enumerate{
#'   \item Independence from treatment (correlation test)
#'   \item Correlation with outcome
#' }
#'
#' Based on these tests, variables are classified as:
#' \describe{
#'   \item{Potential Instrument}{Correlates with treatment but NOT outcome}
#'   \item{Confounder}{Correlates with BOTH treatment and outcome}
#'   \item{Safe Covariate}{No significant correlations}
#' }
#'
#' If negative_controls are specified, they are validated to ensure they
#' do not show spurious treatment effects (negative control logic).
#'
#' @references
#' Akdemir, D. (2026). Constraints on Causal Inference as Experiment Comparison.
#' DOI: 10.5281/zenodo.18367347. See `thm:nc_bound` (Negative Control Sensitivity Bound).
#'
#' @examples
#' # Create sample data with known structure
#' n <- 300
#' U <- rnorm(n)
#' W <- U + rnorm(n, sd = 0.5)  # Confounder
#' Z <- rnorm(n)                 # Potential instrument
#' A <- rbinom(n, 1, plogis(0.5 * Z + 0.3 * U))
#' Y <- 2 * A + 1.5 * U + rnorm(n)
#' df <- data.frame(W = W, Z = Z, A = A, Y = Y)
#'
#' # Run audit
#' report <- audit_data(df, treatment = "A", outcome = "Y")
#' print(report)
#'
#' @seealso [nc_diagnostic()], [causal_spec()]
#' @export
audit_data <- function(data, treatment, outcome, covariates = NULL,
                        negative_controls = NULL, alpha = 0.05,
                        verbose = TRUE) {

  # Input validation
  checkmate::assert_data_frame(data, min.rows = 30,
                                .var.name = "data (need >= 30 observations)")
  checkmate::assert_string(treatment)
  checkmate::assert_string(outcome)
  checkmate::assert_character(covariates, null.ok = TRUE)
  checkmate::assert_character(negative_controls, null.ok = TRUE)
  checkmate::assert_number(alpha, lower = 0.001, upper = 0.5)

  # Check that treatment and outcome exist
  if (!treatment %in% names(data)) {
    .msg_error(paste("Treatment variable", treatment, "not found in data"))
  }
  if (!outcome %in% names(data)) {
    .msg_error(paste("Outcome variable", outcome, "not found in data"))
  }

  # Get variables to audit
  all_vars <- names(data)
  exclude_vars <- c(treatment, outcome)
  
  if (!is.null(covariates)) {
    vars_to_audit <- intersect(covariates, all_vars)
    vars_to_audit <- setdiff(vars_to_audit, exclude_vars)
  } else {
    vars_to_audit <- setdiff(all_vars, exclude_vars)
  }

  if (verbose) {
    .msg_info(paste("Auditing", length(vars_to_audit), "variables..."))
  }

  # Get treatment and outcome vectors
  A <- data[[treatment]]
  Y <- data[[outcome]]

  # Convert factor/character treatment to numeric for correlation
  if (is.factor(A) || is.character(A)) {
    A_numeric <- as.numeric(as.factor(A))
  } else {
    A_numeric <- as.numeric(A)
  }
  Y_numeric <- as.numeric(Y)

  # Audit each variable
  issues_list <- list()

  for (var in vars_to_audit) {
    X <- data[[var]]

    # Skip non-numeric variables for now
    if (!is.numeric(X) && !is.factor(X) && !is.logical(X)) {
      next
    }

    # Convert to numeric for correlation tests
    if (is.factor(X) || is.logical(X)) {
      X_numeric <- as.numeric(X)
    } else {
      X_numeric <- X
    }

    # Test 1: Correlation with treatment
    tryCatch({
      cor_A <- stats::cor.test(X_numeric, A_numeric, method = "pearson")
      p_A <- cor_A$p.value
      r_A <- cor_A$estimate
    }, error = function(e) {
      p_A <- 1
      r_A <- 0
    })

    # Test 2: Correlation with outcome
    tryCatch({
      cor_Y <- stats::cor.test(X_numeric, Y_numeric, method = "pearson")
      p_Y <- cor_Y$p.value
      r_Y <- cor_Y$estimate
    }, error = function(e) {
      p_Y <- 1
      r_Y <- 0
    })

    # Classify based on correlations
    sig_A <- p_A < alpha
    sig_Y <- p_Y < alpha

    issue_type <- NULL
    recommendation <- NULL
    p_value <- NA

    if (sig_A && sig_Y) {
      # Correlated with both: Confounder
      issue_type <- "Confounder"
      recommendation <- "Include in adjustment set"
      p_value <- max(p_A, p_Y)
    } else if (sig_A && !sig_Y) {
      # Correlated with A but not Y: Potential Instrument
      issue_type <- "Potential Instrument"
      recommendation <- "Consider using as instrumental variable"
      p_value <- p_A
    } else if (!sig_A && sig_Y) {
      # Correlated with Y but not A: Predictor only
      issue_type <- "Outcome Predictor"
      recommendation <- "May improve precision if included"
      p_value <- p_Y
    } else {
      # No significant correlations: Safe
      issue_type <- "Safe"
      recommendation <- "No causal concerns detected"
      p_value <- max(p_A, p_Y)
    }

    # Check if this is a specified negative control
    is_nc <- var %in% negative_controls
    if (is_nc) {
      # For NC, we expect NO treatment correlation after adjustment
      # Here we just test raw correlation as a baseline
      if (sig_A) {
        issue_type <- "NC Violation"
        recommendation <- "Negative control shows treatment correlation - do not use as control"
        p_value <- p_A
      } else {
        issue_type <- "Valid Negative Control"
        recommendation <- "Can be used for negative control diagnostic"
        p_value <- p_A
      }
    }

    issues_list[[var]] <- data.frame(
      variable = var,
      issue_type = issue_type,
      p_value = p_value,
      r_treatment = r_A,
      r_outcome = r_Y,
      recommendation = recommendation,
      stringsAsFactors = FALSE
    )
  }

  # Combine issues
  if (length(issues_list) > 0) {
    issues <- do.call(rbind, issues_list)
    rownames(issues) <- NULL
  } else {
    issues <- data.frame(
      variable = character(),
      issue_type = character(),
      p_value = numeric(),
      r_treatment = numeric(),
      r_outcome = numeric(),
      recommendation = character(),
      stringsAsFactors = FALSE
    )
  }

  # Generate recommendations summary
  recommendations <- .generate_audit_recommendations(issues)

  # Summary statistics
  summary_stats <- list(
    n_vars_audited = length(vars_to_audit),
    n_issues = sum(!issues$issue_type %in% c("Safe", "Valid Negative Control")),
    n_confounders = sum(issues$issue_type == "Confounder"),
    n_instruments = sum(issues$issue_type == "Potential Instrument"),
    n_nc_violations = sum(issues$issue_type == "NC Violation"),
    treatment = treatment,
    outcome = outcome,
    alpha = alpha
  )

  if (verbose) {
    .msg_success(paste("Audit complete:", summary_stats$n_issues, "potential issues found"))
  }

  new_data_audit_report(
    issues = issues,
    recommendations = recommendations,
    summary_stats = summary_stats,
    treatment = treatment,
    outcome = outcome
  )
}

# =============================================================================
# Internal: Generate Recommendations
# =============================================================================

#' @keywords internal
.generate_audit_recommendations <- function(issues) {

  recommendations <- character()

  # Confounders
  confounders <- issues$variable[issues$issue_type == "Confounder"]
  if (length(confounders) > 0) {
    recommendations <- c(
      recommendations,
      paste0("CONFOUNDERS: Variables [", paste(confounders, collapse = ", "),
             "] correlate with both treatment and outcome - must adjust for these")
    )
  }

  # Instruments
  instruments <- issues$variable[issues$issue_type == "Potential Instrument"]
  if (length(instruments) > 0) {
    recommendations <- c(
      recommendations,
      paste0("INSTRUMENTS: Variables [", paste(instruments, collapse = ", "),
             "] correlate with treatment but not outcome - consider as IVs")
    )
  }

  # NC violations
  nc_violations <- issues$variable[issues$issue_type == "NC Violation"]
  if (length(nc_violations) > 0) {
    recommendations <- c(
      recommendations,
      paste0("NC VIOLATIONS: Variables [", paste(nc_violations, collapse = ", "),
             "] fail negative control test - do not use as controls")
    )
  }

  if (length(recommendations) == 0) {
    recommendations <- "No significant causal validity issues detected."
  }

  recommendations
}
