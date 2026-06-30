# R/plots/plot_scores.R
# Factor score visualizations (paper Figures 3, 10).

#' Scatter plot of two factor scores, coloured by a sample annotation
#'
#' @param scores  Matrix (k x n) or (n x k) of factor scores.
#' @param factors Integer vector of length 2: which factors to plot.
#' @param group   Optional factor/character vector (length n) for point colour.
#' @param method  Method name used in the plot title.
#' @return A ggplot object.
plot_factor_scatter <- function(scores, factors = c(1, 2), group = NULL, method = "") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  S <- as.matrix(scores)
  if (nrow(S) < ncol(S)) S <- t(S)   # coerce to (n x k)

  df <- data.frame(
    Factor_x = S[, factors[1]],
    Factor_y = S[, factors[2]],
    group    = if (is.null(group)) "sample" else as.factor(group)
  )

  ggplot2::ggplot(df, ggplot2::aes(Factor_x, Factor_y, colour = group)) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::labs(
      title = sprintf("%s factor scores", method),
      x = paste0("Factor ", factors[1]),
      y = paste0("Factor ", factors[2]),
      colour = NULL
    ) +
    ggplot2::theme_bw()
}

#' Save a factor-scatter plot to disk
#'
#' @param plot   ggplot object.
#' @param path   Output file path (.png/.pdf).
#' @param width,height Dimensions in inches.
save_plot <- function(plot, path, width = 7, height = 5) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = 300)
  invisible(path)
}
