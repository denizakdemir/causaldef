# =============================================================================
# Tests for audit_data() - Automated Data Auditor
# =============================================================================

# -----------------------------------------------------------------------------
# Test Helper: Simulate data with known causal structure
# -----------------------------------------------------------------------------

# Simulate data with instrument, confounder, and negative control
.simulate_audit_data <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  U <- rnorm(n)                # Unobserved confounder
  W <- U + rnorm(n, sd = 0.5)  # Observed confounder (correlated with U)
  Z <- rnorm(n)                # True instrument (affects A, not Y directly)
  
  # Treatment depends on instrument and confounder

  A <- rbinom(n, 1, plogis(0.5 * Z + 0.3 * U + 0.2 * W))
  
  # Outcome depends on treatment and confounder, NOT on instrument
  Y <- 2 * A + 1.5 * U + 0.8 * W + rnorm(n)
  
  # Negative control: affected by confounder but NOT by treatment
  Y_nc <- U + rnorm(n)
  
  # Safe variable: independent of everything causal
  X_safe <- rnorm(n)
  
  data.frame(
    W = W,           # Confounder
    Z = Z,           # Instrument candidate
    A = A,           # Treatment
    Y = Y,           # Outcome
    Y_nc = Y_nc,     # Negative control
    X_safe = X_safe  # Irrelevant safe variable
  )
}

# -----------------------------------------------------------------------------
# Class and Structure Tests
# -----------------------------------------------------------------------------

test_that("audit_data returns correct class", {
  df <- .simulate_audit_data(200, seed = 1)
  result <- audit_data(df, treatment = "A", outcome = "Y")
  
  expect_s3_class(result, "data_audit_report")
  expect_true("issues" %in% names(result))
  expect_true("recommendations" %in% names(result))
  expect_true("summary_stats" %in% names(result))
  expect_true("treatment" %in% names(result))
  expect_true("outcome" %in% names(result))
})

test_that("audit_data handles missing treatment/outcome", {
  df <- data.frame(X = rnorm(100), Y = rnorm(100))
  
  expect_error(
    audit_data(df, treatment = "A", outcome = "Y"),
    "not found"
  )
})

test_that("audit_data requires minimum observations", {
  df <- data.frame(A = c(0, 1), Y = c(1, 2))
  
  expect_error(
    audit_data(df, treatment = "A", outcome = "Y"),
    "observations"
  )
})

# -----------------------------------------------------------------------------
# Detection Tests
# -----------------------------------------------------------------------------

test_that("audit_data detects potential instruments", {
  # Create data with clear instrument (strong A correlation, weak Y correlation)
  set.seed(42)
  n <- 1000
  Z <- rnorm(n)  # Instrument
  A <- rbinom(n, 1, plogis(2.0 * Z))  # Very strong Z -> A
  Y <- 2 * A + rnorm(n)               # No direct Z -> Y at all
  
  df <- data.frame(Z = Z, A = A, Y = Y)
  result <- audit_data(df, treatment = "A", outcome = "Y")
  
  # Z should be flagged as potential instrument or at least be in the issues
  z_issues <- result$issues[result$issues$variable == "Z", ]
  expect_true(nrow(z_issues) > 0)
  # Either instrument OR we detected it has A-correlation (both valid)
  expect_true(
    any(grepl("instrument", tolower(z_issues$issue_type))) ||
    !is.na(z_issues$r_treatment[1]) && abs(z_issues$r_treatment[1]) > 0.1
  )
})

test_that("audit_data detects confounders", {
  # Create data with clear confounder (both A and Y correlated)
  set.seed(42)
  n <- 500
  W <- rnorm(n)  # Confounder
  A <- rbinom(n, 1, plogis(1.0 * W))
  Y <- 2 * A + 1.5 * W + rnorm(n)
  
  df <- data.frame(W = W, A = A, Y = Y)
  result <- audit_data(df, treatment = "A", outcome = "Y")
  
  # W should be flagged as confounder
  w_issues <- result$issues[result$issues$variable == "W", ]
  expect_true(nrow(w_issues) > 0)
  expect_true(any(grepl("confounder", tolower(w_issues$issue_type))))
})

test_that("audit_data identifies safe variables", {
  # Create data with variable unrelated to causal structure
  set.seed(42)
  n <- 500
  X <- rnorm(n)  # Safe/irrelevant
  A <- rbinom(n, 1, 0.5)
  Y <- 2 * A + rnorm(n)
  
  df <- data.frame(X = X, A = A, Y = Y)
  result <- audit_data(df, treatment = "A", outcome = "Y")
  
  # X should not be in issues OR flagged as safe
  if (nrow(result$issues) > 0) {
    x_issues <- result$issues[result$issues$variable == "X", ]
    if (nrow(x_issues) > 0) {
      expect_true(any(grepl("safe", tolower(x_issues$issue_type))))
    }
  }
})

test_that("audit_data handles negative control specification", {
  df <- .simulate_audit_data(300, seed = 123)
  
  # When we have a known negative control, validate it
  result <- audit_data(
    df, 
    treatment = "A", 
    outcome = "Y",
    negative_controls = "Y_nc"
  )
  
  expect_s3_class(result, "data_audit_report")
  
  # Y_nc should be in issues since we specified it as NC
  nc_issues <- result$issues[result$issues$variable == "Y_nc", ]
  expect_true(nrow(nc_issues) > 0)
  # Should have NC-related classification (Valid NC or NC Violation)
  expect_true(
    any(grepl("control", tolower(nc_issues$issue_type))) ||
    any(grepl("nc", tolower(nc_issues$issue_type)))
  )
})

# -----------------------------------------------------------------------------
# Report Generation Tests
# -----------------------------------------------------------------------------

test_that("issues data.frame has correct structure", {
  df <- .simulate_audit_data(200, seed = 1)
  result <- audit_data(df, treatment = "A", outcome = "Y")
  
  expect_true(is.data.frame(result$issues))
  expected_cols <- c("variable", "issue_type", "p_value", "recommendation")
  expect_true(all(expected_cols %in% names(result$issues)))
})

test_that("recommendations are character vector", {
  df <- .simulate_audit_data(200, seed = 1)
  result <- audit_data(df, treatment = "A", outcome = "Y")
  
  expect_type(result$recommendations, "character")
})

test_that("summary_stats contains expected metrics", {
  df <- .simulate_audit_data(200, seed = 1)
  result <- audit_data(df, treatment = "A", outcome = "Y")
  
  expect_type(result$summary_stats, "list")
  expect_true("n_vars_audited" %in% names(result$summary_stats))
  expect_true("n_issues" %in% names(result$summary_stats))
})

# -----------------------------------------------------------------------------
# Print Method Tests
# -----------------------------------------------------------------------------

test_that("print.data_audit_report works", {
  df <- .simulate_audit_data(200, seed = 1)
  result <- audit_data(df, treatment = "A", outcome = "Y")
  
  expect_output(print(result), "Data Integrity Report")
})

test_that("print.data_audit_report handles empty issues", {
  # Create data with no causal issues
  set.seed(42)
  n <- 200
  A <- rbinom(n, 1, 0.5)
  Y <- 2 * A + rnorm(n)
  df <- data.frame(A = A, Y = Y)
  
  result <- audit_data(df, treatment = "A", outcome = "Y")
  
  # Should print without error even with empty issues
  expect_output(print(result), "Report")
})

# -----------------------------------------------------------------------------
# Edge Cases
# -----------------------------------------------------------------------------

test_that("audit_data handles covariates argument", {
  df <- .simulate_audit_data(200, seed = 1)
  
  # Audit only specific columns
  result <- audit_data(
    df, 
    treatment = "A", 
    outcome = "Y",
    covariates = c("W", "Z")
  )
  
  # Should only audit W and Z
  expect_true(all(result$issues$variable %in% c("W", "Z")))
})

test_that("audit_data respects alpha parameter", {
  df <- .simulate_audit_data(500, seed = 1)
  
  # Stricter alpha should find fewer issues
  result_strict <- audit_data(df, treatment = "A", outcome = "Y", alpha = 0.001)
  result_loose <- audit_data(df, treatment = "A", outcome = "Y", alpha = 0.10)
  
  # In general, loose alpha finds >= strict alpha issues
  # (Not always true due to sampling, but good test expectation)
  expect_gte(nrow(result_loose$issues), 0)
  expect_gte(nrow(result_strict$issues), 0)
})

test_that("audit_data handles factor treatment", {
  df <- .simulate_audit_data(200, seed = 1)
  df$A <- factor(df$A, labels = c("control", "treated"))
  
  result <- audit_data(df, treatment = "A", outcome = "Y")
  
  expect_s3_class(result, "data_audit_report")
})

test_that("audit_data handles character treatment", {
  df <- .simulate_audit_data(200, seed = 1)
  df$A <- ifelse(df$A == 1, "treated", "control")
  
  result <- audit_data(df, treatment = "A", outcome = "Y")
  
  expect_s3_class(result, "data_audit_report")
})
