# =============================================================================
# Utility Functions for causaldef
# =============================================================================

# Declare global variables used in ggplot2 aesthetics (avoids R CMD check NOTEs)
utils::globalVariables(c("alpha", "delta", "estimate", "lower", "upper", 
                         "method", "outcome", "Value", "Component", "safety",
                         "transfer", "minimax",
                         "delta_bound",
                         ".weights", "gamma"))

# Cross-compatible message functions
# Support both modern cli and older base R

.msg_info <- function(msg) {
  if (requireNamespace("cli", quietly = TRUE) && 
      exists("cli_alert_info", envir = asNamespace("cli"))) {
    cli::cli_alert_info(msg)
  } else {
    message("i ", msg)
  }
  invisible(NULL)
}

.msg_success <- function(msg) {
  if (requireNamespace("cli", quietly = TRUE) && 
      exists("cli_alert_success", envir = asNamespace("cli"))) {
    cli::cli_alert_success(msg)
  } else {
    message("v ", msg)
  }
  invisible(NULL)
}

.msg_warning <- function(msg) {
  if (requireNamespace("cli", quietly = TRUE) && 
      exists("cli_warn", envir = asNamespace("cli"))) {
    cli::cli_warn(msg)
  } else {
    warning(msg, call. = FALSE)
  }
  invisible(NULL)
}

.msg_error <- function(msg) {
  if (requireNamespace("cli", quietly = TRUE) && 
      exists("cli_abort", envir = asNamespace("cli"))) {
    cli::cli_abort(msg)
  } else {
    stop(msg, call. = FALSE)
  }
}

.msg_danger <- function(msg) {
  if (requireNamespace("cli", quietly = TRUE) && 
      exists("cli_alert_danger", envir = asNamespace("cli"))) {
    cli::cli_alert_danger(msg)
  } else {
    message("x ", msg)
  }
  invisible(NULL)
}

.msg_header <- function(msg) {
  if (requireNamespace("cli", quietly = TRUE) && 
      exists("cli_h1", envir = asNamespace("cli"))) {
    cli::cli_h1(msg)
  } else {
    cat("\n--", msg, paste(rep("-", max(0, 60 - nchar(msg))), collapse = ""), "\n\n")
  }
  invisible(NULL)
}

.msg_list <- function(items) {
  if (requireNamespace("cli", quietly = TRUE) && 
      exists("cli_ul", envir = asNamespace("cli"))) {
    cli::cli_ul(items)
  } else {
    cat(paste0("* ", items, "\n"), sep = "")
  }
  invisible(NULL)
}

.survival_runtime_supported <- function() {
  getRversion() >= "4.0.0" &&
    exists("deparse1", envir = baseenv(), inherits = FALSE)
}

.require_survival_runtime <- function(feature = "Survival estimation") {
  if (!requireNamespace("survival", quietly = TRUE)) {
    .msg_error(paste0("Package 'survival' needed for ", feature, "."))
  }

  if (!.survival_runtime_supported()) {
    .msg_error(paste0(
      feature,
      " requires R >= 4.0 or a survival installation compatible with base::deparse1. ",
      "Current runtime: R ", getRversion(), "."
    ))
  }

  invisible(TRUE)
}
