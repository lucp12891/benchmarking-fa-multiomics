# R/plots/plot_loadings.R
# Factor loading visualizations: heatmaps and lollipop plots (paper Figures 4, 11).

#' Heatmap of factor loadings for one omics layer
#'
#' @param loadings Matrix (features x k) of loadings.
#' @param top_n    Show only the top_n features by maximum absolute loading.
#' @param title    Plot title.
#' @return Invisibly, the pheatmap object.
plot_loading_heatmap <- function(loadings, top_n = 50, title = "Loadings") {
  if (!requireNamespace("pheatmap", quietly = TRUE)) stop("pheatmap is required.")
  L <- as.matrix(loadings)
  rank_by <- order(apply(abs(L), 1, max), decreasing = TRUE)
  L <- L[head(rank_by, min(top_n, nrow(L))), , drop = FALSE]
  colnames(L) <- paste0("F", seq_len(ncol(L)))
  pheatmap::pheatmap(L, main = title, cluster_cols = FALSE,
                     color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                     silent = TRUE)
}

#' Lollipop plot of the top loadings for a single factor
#'
#' @param loadings Numeric vector of loadings (named with feature ids).
#' @param top_n    Number of largest-magnitude features to show.
#' @param title    Plot title.
#' @return A ggplot object.
plot_loading_lollipop <- function(loadings, top_n = 20, title = "Top loadings") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  v   <- loadings[order(abs(loadings), decreasing = TRUE)][seq_len(min(top_n, length(loadings)))]
  df  <- data.frame(feature = factor(names(v), levels = names(v)[order(v)]),
                    loading = as.numeric(v))

  ggplot2::ggplot(df, ggplot2::aes(x = loading, y = feature)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = loading,
                                       y = feature, yend = feature),
                          colour = "grey60") +
    ggplot2::geom_point(ggplot2::aes(colour = loading > 0), size = 3) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
                                 guide = "none") +
    ggplot2::labs(title = title, x = "Loading", y = NULL) +
    ggplot2::theme_bw()
}
