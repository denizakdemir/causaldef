# =============================================================================
# Utility Functions for causaldef
# =============================================================================

# Declare global variables used in ggplot2 aesthetics (avoids R CMD check NOTEs)
utils::globalVariables(c("alpha", "delta", "estimate", "lower", "upper", 
                         "method", "outcome", "Value", "Component", "safety", 
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
