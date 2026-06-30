# R/plots/plot_simulation.R
# Simulation performance visualizations (paper Figures 12-13).

#' Plot a performance metric versus noise level, by method and cutoff
#'
#' @param sim_df Data frame of simulation results (from run_full_simulation),
#'               with columns: method, sigma, scenario, and the metric column.
#' @param metric Character. Metric column to plot (e.g., "jaccard", "F1", "pearson").
#' @param facet_by Character. Column to facet by (default "scenario").
#' @return A ggplot object.
plot_metric_vs_noise <- function(sim_df, metric = "jaccard", facet_by = "scenario") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  if (!requireNamespace("dplyr", quietly = TRUE))   stop("dplyr is required.")

  summ <- dplyr::summarise(
    dplyr::group_by(sim_df, method, sigma, .data[[facet_by]]),
    mean_val = mean(.data[[metric]], na.rm = TRUE),
    se_val   = stats::sd(.data[[metric]], na.rm = TRUE) / sqrt(dplyr::n()),
    .groups  = "drop"
  )

  ggplot2::ggplot(summ, ggplot2::aes(x = sigma, y = mean_val, colour = method)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = mean_val - se_val,
                                      ymax = mean_val + se_val,
                                      fill = method), alpha = 0.15, colour = NA) +
    ggplot2::facet_wrap(facet_by) +
    ggplot2::labs(x = expression(sigma~"(noise level)"),
                  y = paste("Mean", metric),
                  colour = "Method", fill = "Method") +
    ggplot2::theme_bw()
}

#' Boxplot of a metric across methods at a fixed noise level
#'
#' @param sim_df Simulation results data frame.
#' @param metric Metric column to plot.
#' @param sigma_value Noise level to subset to.
#' @return A ggplot object.
plot_metric_box <- function(sim_df, metric = "jaccard", sigma_value = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  df <- if (is.null(sigma_value)) sim_df else sim_df[sim_df$sigma == sigma_value, ]

  ggplot2::ggplot(df, ggplot2::aes(x = method, y = .data[[metric]], fill = method)) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    ggplot2::labs(x = NULL, y = metric,
                  title = if (!is.null(sigma_value)) paste("sigma =", sigma_value) else NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
}
