# R/plots/plot_jaccard.R
# Jaccard-index agreement curves and intersection plots (paper Figures 5-7).

#' Plot Jaccard index between method pairs across a sequence of cutoffs
#'
#' @param jaccard_df Data frame from \code{jaccard_across_cutoffs()} with columns
#'                   method1, method2, jaccard, cutoff.
#' @param title      Plot title.
#' @return A ggplot object.
plot_jaccard_curves <- function(jaccard_df, title = "Pairwise Jaccard agreement") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  jaccard_df$pair <- paste(jaccard_df$method1, jaccard_df$method2, sep = " vs ")

  ggplot2::ggplot(jaccard_df,
                  ggplot2::aes(x = cutoff, y = jaccard, colour = pair)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(title = title, x = "Quantile cutoff", y = "Jaccard index",
                  colour = "Method pair") +
    ggplot2::theme_bw()
}

#' Upset-style intersection plot of selected features across methods
#'
#' @param selected_list Named list of character vectors (selected feature ids per method).
#' @param title         Plot title.
#' @return A ggplot object (requires ggupset).
plot_feature_intersections <- function(selected_list, title = "Selected-feature intersections") {
  if (!requireNamespace("ggupset", quietly = TRUE)) stop("ggupset is required.")
  all_feats <- unique(unlist(selected_list))
  df <- data.frame(feature = all_feats)
  df$methods <- lapply(all_feats, function(f) {
    names(selected_list)[vapply(selected_list, function(s) f %in% s, logical(1))]
  })

  ggplot2::ggplot(df, ggplot2::aes(x = methods)) +
    ggplot2::geom_bar() +
    ggupset::scale_x_upset() +
    ggplot2::labs(title = title, x = NULL, y = "Shared features") +
    ggplot2::theme_bw()
}
