# R/plots/plot_fabia_sparsity.R
# FABIA sparsity-parameter grid visualization (paper Figure 8).

#' Heatmap/tile plot of FABIA sparsity over the (spz, spl) grid at fixed alpha
#'
#' @param grid_df Data frame from \code{run_fabia_sparsity_grid()} with columns
#'                alpha, spz, spl, prop_nonzero_scores, prop_nonzero_<omic>...
#' @param value   Column to display (e.g., "prop_nonzero_scores").
#' @param alpha_value Optional alpha to subset to (defaults to all, faceted).
#' @return A ggplot object.
plot_fabia_sparsity <- function(grid_df, value = "prop_nonzero_scores", alpha_value = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  df <- if (is.null(alpha_value)) grid_df else grid_df[grid_df$alpha == alpha_value, ]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = factor(spz), y = factor(spl),
                                        fill = .data[[value]])) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::scale_fill_viridis_c(option = "magma", limits = c(0, 1)) +
    ggplot2::labs(x = "spz (score sparsity)", y = "spl (loading sparsity)",
                  fill = "Prop.\nnon-zero",
                  title = "FABIA sparsity grid") +
    ggplot2::theme_minimal()

  if (is.null(alpha_value)) p <- p + ggplot2::facet_wrap(~ alpha)
  p
}
