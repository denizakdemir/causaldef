# =============================================================================
# Plot Methods for causaldef Classes
# =============================================================================

#' Plot Deficiency Estimates
#'
#' @param x A deficiency object
#' @param type Character: plot type ("bar", "forest")
#' @param ... Additional arguments passed to ggplot2 functions
#'
#' @return A ggplot2 object
#' @export
plot.deficiency <- function(x, type = c("bar", "forest"), ...) {
  type <- match.arg(type)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggplot2} is required for plotting")
  }
  
  df <- data.frame(
    method = factor(names(x$estimates), levels = names(x$estimates)),
    delta = x$estimates,
    se = x$se,
    stringsAsFactors = FALSE
  )
  
  if (!all(is.na(x$ci))) {
    df$lower <- x$ci[, 1]
    df$upper <- x$ci[, 2]
  }
  
  if (type == "bar") {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = method, y = delta, fill = method)) +
      ggplot2::geom_col(width = 0.6, show.legend = FALSE) +
      ggplot2::labs(
        title = "Le Cam Deficiency by Adjustment Method",
        x = "Method",
        y = expression(delta)
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      ) +
      # Traffic Light Thresholds
      ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", color = "darkgreen", alpha = 0.5) +
      ggplot2::geom_hline(yintercept = 0.10, linetype = "dashed", color = "red", alpha = 0.5) +
      ggplot2::annotate("text", x = 0.6, y = 0.05, label = "Excellent (<0.05)", color = "darkgreen", vjust = 1.2, size = 3) +
      ggplot2::annotate("text", x = 0.6, y = 0.10, label = "Caution (<0.10)", color = "red", vjust = 1.2, size = 3)
    
    if (!all(is.na(x$ci))) {
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = lower, ymax = upper),
        width = 0.2
      )
    }
  } else if (type == "forest") {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = delta, y = method)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
      ggplot2::labs(
        title = "Le Cam Deficiency: Forest Plot",
        x = expression(delta),
        y = "Method"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
      # Traffic Light Thresholds
      ggplot2::geom_vline(xintercept = 0.05, linetype = "dashed", color = "darkgreen", alpha = 0.5) +
      ggplot2::geom_vline(xintercept = 0.10, linetype = "dashed", color = "red", alpha = 0.5)
    
    if (!all(is.na(x$ci))) {
      p <- p + ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = lower, xmax = upper),
        height = 0.2
      )
    }
  }
  
  p
}

#' Plot Confounding Frontier
#'
#' @param x A confounding_frontier object
#' @param type Character: plot type ("heatmap", "contour")
#' @param threshold Numeric vector: threshold values for contour lines
#' @param ... Additional arguments
#'
#' @return A ggplot2 object
#' @export
plot.confounding_frontier <- function(x, type = c("heatmap", "contour"),
                                      threshold = c(0.05, 0.1, 0.2), ...) {
  type <- match.arg(type)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggplot2} is required for plotting")
  }
  
  grid <- x$grid
  
  if (type == "heatmap") {
    p <- ggplot2::ggplot(grid, ggplot2::aes(x = alpha, y = gamma, fill = delta)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(
        low = "white", mid = "orange", high = "darkred",
        midpoint = 0.2,
        name = expression(delta)
      ) +
      ggplot2::labs(
        title = "Confounding Frontier (Theorem 4.1)",
        x = expression(alpha ~ "(U" %->% "A)"),
        y = expression(gamma ~ "(U" %->% "Y)")
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.position = "right"
      ) +
      ggplot2::coord_fixed()
    
    # Add contour lines at thresholds
    if (!is.null(threshold)) {
      p <- p + ggplot2::geom_contour(
        ggplot2::aes(z = delta),
        breaks = threshold,
        color = "black",
        linetype = "dashed"
      )
    }
  } else if (type == "contour") {
    p <- ggplot2::ggplot(grid, ggplot2::aes(x = alpha, y = gamma, z = delta)) +
      ggplot2::geom_contour_filled() +
      ggplot2::labs(
        title = "Confounding Frontier Contours",
        x = expression(alpha ~ "(U" %->% "A)"),
        y = expression(gamma ~ "(U" %->% "Y)"),
        fill = expression(delta)
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
      ggplot2::coord_fixed()
  }
  
  # Overlay Benchmarks if available
  if (!is.null(x$benchmarks)) {
    # Ensure they are within plot range to avoid warnings
    bm <- x$benchmarks
    
    p <- p + 
      # Benchmark Points
      ggplot2::geom_point(
        data = bm,
        ggplot2::aes(x = alpha, y = gamma),
        inherit.aes = FALSE,
        color = "black", 
        fill = "yellow",
        shape = 21,
        size = 3
      ) +
      # Benchmark Labels
      ggplot2::geom_text(
        data = bm,
        ggplot2::aes(x = alpha, y = gamma, label = covariate),
        inherit.aes = FALSE,
        vjust = -1,
        size = 3,
        color = "black",
        fontface = "bold",
        check_overlap = TRUE
      ) +
      ggplot2::labs(
        caption = "Yellow points indicate observed covariate strength (if they were unobserved confounders)"
      )
  }
  
  p
}

#' Plot Policy Regret Bound
#'
#' @param x A policy_bound object
#' @param type Character: plot type
#' @param ... Additional arguments
#'
#' @return A ggplot2 object
#' @export
plot.policy_bound <- function(x, type = c("decomposition", "safety_curve"), ...) {
  type <- match.arg(type)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggplot2} is required for plotting")
  }
  
  if (type == "decomposition") {
    # Bar chart decomposing the regret bound
    components <- data.frame(
      Component = factor(c("Observed Regret", "Safety Floor"),
                        levels = c("Observed Regret", "Safety Floor")),
      Value = c(
        if (!is.null(x$obs_regret)) x$obs_regret else 0,
        x$safety_floor
      )
    )
    
    p <- ggplot2::ggplot(components, ggplot2::aes(x = 1, y = Value, fill = Component)) +
      ggplot2::geom_col(width = 0.5) +
      ggplot2::labs(
        title = "Policy Regret Decomposition",
        y = "Regret",
        x = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      ) +
      ggplot2::scale_fill_manual(values = c("steelblue", "coral"))
    
  } else if (type == "safety_curve") {
    # Safety floor as function of delta
    delta_seq <- seq(0, 1, length.out = 100)
    safety_seq <- 2 * x$M * delta_seq
    
    df <- data.frame(delta = delta_seq, safety = safety_seq)
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = delta, y = safety)) +
      ggplot2::geom_line(linewidth = 1.2, color = "darkred") +
      ggplot2::labs(
        title = "Safety Floor vs Deficiency",
        x = expression(delta),
        y = paste0("Safety Floor (M = ", x$M, ")")
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

    # Add points
    if (!is.null(x$all_estimates)) {
      points_df <- data.frame(
        method = names(x$all_estimates),
        delta = as.numeric(x$all_estimates),
        safety = 2 * x$M * as.numeric(x$all_estimates)
      )
      
      p <- p + ggplot2::geom_point(
        data = points_df,
        ggplot2::aes(x = delta, y = safety, color = method),
        size = 4
      ) + 
      # Add text labels
      ggplot2::geom_text(
         data = points_df,
         ggplot2::aes(x = delta, y = safety, label = method),
         vjust = -1, size = 3, show.legend = FALSE
      )
      
    } else {
      # Single point fallback
      p <- p + ggplot2::geom_point(
        data = data.frame(delta = x$delta, safety = x$safety_floor),
        size = 4, color = "blue"
      )
    }
  }
  
  p
}

#' Plot Causal Effect
#'
#' @param x A causal_effect object
#' @param ... Additional arguments passed to plot functions
#'
#' @return A plot object (base or ggplot)
#' @export
plot.causal_effect <- function(x, ...) {
  
  if (!is.null(x$fit) && inherits(x$fit, "survfit")) {
    # Survival Curve Plot
    # We use base R plot for survfit as it is robust
    
    # Extract method name for title
    method_title <- toupper(x$method)
    
    # Plot
    plot(x$fit, col = c("blue", "red"), lwd = 2,
         xlab = "Time", ylab = "Survival Probability",
         main = paste("Survival Curves (", method_title, ")", sep=""), ...)
    
    legend("topright", legend = c("Treated", "Control"), 
           col = c("blue", "red"), lwd = 2)
    
    # If RMST is involved, maybe add a vertical line?
    if (!is.null(x$horizon)) {
      abline(v = x$horizon, lty = 2, col = "gray")
      text(x$horizon, 0.1, paste("Horizon =", x$horizon), pos = 4)
    }
    
  } else {
    # Standard effect plot (Point estimate + CI if available)
    # Since we don't have CI always, just plot the point estimate on a forest-like plot
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      cli::cli_abort("Package {.pkg ggplot2} is required for plotting")
    }
    
    df <- data.frame(
      outcome = "Outcome",
      estimate = x$estimate,
      lower = if (!is.null(x$ci)) x$ci[1] else x$estimate,
      upper = if (!is.null(x$ci)) x$ci[2] else x$estimate
    )
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = outcome, y = estimate)) +
      ggplot2::geom_point(size = 4, color = "darkblue") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      ggplot2::labs(
        title = paste("Causal Effect:", x$type),
        subtitle = paste("Method:", x$method),
        y = "Estimate",
        x = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::coord_flip()
    
    if (!is.null(x$ci)) {
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = lower, ymax = upper),
        width = 0.1
      )
    }
    
    return(p)
  }
}
