# scripts/04_generate_plots.R
# Reproduce all main-text figures from Osang'ir et al. (2026)

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

source("R/utils/metrics.R")
source("R/utils/thresholding.R")

FIGURE_DIR <- "figures"
dir.create(FIGURE_DIR, showWarnings = FALSE)

theme_bench <- function() {
  theme_classic(base_size = 11) +
    theme(
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold"),
      legend.position  = "right",
      plot.title       = element_text(face = "bold", size = 12)
    )
}

METHOD_COLORS <- c(
  FABIA = "#E63946",
  MOFA  = "#2A9D8F",
  GFA   = "#E9C46A",
  MFA   = "#264653"
)

# ─────────────────────────────────────────────
# Figure 3B — Factor 1 scores annotated by IGHV
# ─────────────────────────────────────────────
plot_factor_scores_cll <- function(aligned_scores, metadata) {
  df <- lapply(names(aligned_scores), function(m) {
    scores <- aligned_scores[[m]]
    data.frame(
      method = m,
      sample = seq_len(nrow(scores)),
      score  = scores[, 1],
      IGHV   = metadata$IGHV
    )
  })
  df <- do.call(rbind, df)
  df$method <- factor(df$method, levels = c("MOFA", "FABIA", "MFA", "GFA"))

  ggplot(df, aes(x = sample, y = score, colour = factor(IGHV))) +
    geom_point(size = 1.2, alpha = 0.8) +
    facet_wrap(~method, nrow = 2) +
    scale_colour_manual(
      values = c("0" = "#2A9D8F", "1" = "#E63946", "NA" = "grey70"),
      name   = "IGHV"
    ) +
    labs(title = "Factor 1 scores — CLL study",
         x = "Sample", y = "Factor 1 score") +
    theme_bench()
}

# ─────────────────────────────────────────────
# Figure 3C/D — Pairwise Pearson / cosine heatmap
# ─────────────────────────────────────────────
plot_pairwise_heatmap <- function(metric_df, metric = "pearson", title = "") {
  methods <- unique(c(metric_df$method1, metric_df$method2))
  mat     <- matrix(1, nrow = length(methods), ncol = length(methods),
                    dimnames = list(methods, methods))

  for (i in seq_len(nrow(metric_df))) {
    m1 <- metric_df$method1[i]; m2 <- metric_df$method2[i]
    v  <- metric_df[[metric]][i]
    mat[m1, m2] <- v; mat[m2, m1] <- v
  }

  df_long <- as.data.frame(as.table(mat))
  colnames(df_long) <- c("Method1", "Method2", "value")

  ggplot(df_long, aes(x = Method1, y = Method2, fill = value, label = round(value, 2))) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(size = 3.5, fontface = "bold") +
    scale_fill_gradient2(low = "white", high = "#2A9D8F", mid = "#E9C46A",
                         midpoint = 0.9, limits = c(0.7, 1),
                         name = metric) +
    labs(title = title) +
    theme_bench() +
    theme(axis.title = element_blank(),
          legend.position = "right")
}

# ─────────────────────────────────────────────
# Figures 5–7 — Jaccard curves across quantile cutoffs
# ─────────────────────────────────────────────
plot_jaccard_curves <- function(jaccard_df, title = "") {
  jaccard_df$pair <- paste0(jaccard_df$method1, "–", jaccard_df$method2)

  ggplot(jaccard_df, aes(x = cutoff, y = jaccard, colour = pair, group = pair)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.5) +
    scale_x_continuous(labels = scales::percent, limits = c(0.50, 0.95)) +
    scale_y_continuous(limits = c(0, 1), labels = scales::number_format(accuracy = 0.01)) +
    scale_colour_brewer(palette = "Dark2", name = "Method pair") +
    labs(title = title, x = "Quantile cutoff", y = "Jaccard index") +
    theme_bench()
}

# ─────────────────────────────────────────────
# Figure 13 — Simulation: Jaccard vs noise level
# ─────────────────────────────────────────────
plot_simulation_jaccard <- function(sim_results) {
  df <- sim_results %>%
    group_by(method, sigma, scenario) %>%
    summarise(mean_jaccard = mean(jaccard, na.rm = TRUE),
              se_jaccard   = sd(jaccard, na.rm = TRUE) / sqrt(n()),
              .groups = "drop")

  ggplot(df, aes(x = sigma, y = mean_jaccard,
                 colour = method, fill = method, group = method)) +
    geom_ribbon(aes(ymin = mean_jaccard - se_jaccard,
                    ymax = mean_jaccard + se_jaccard), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.8) +
    facet_wrap(~scenario, labeller = label_both) +
    scale_colour_manual(values = METHOD_COLORS) +
    scale_fill_manual(values = METHOD_COLORS) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "Factor loading recovery across noise levels",
         x = expression(sigma ~ "(noise)"), y = "Jaccard index (mean ± SE)",
         colour = "Method", fill = "Method") +
    theme_bench()
}

# ─────────────────────────────────────────────
# Run and save all figures
# ─────────────────────────────────────────────
generate_all_figures <- function() {
  cll_results  <- readRDS("results/cll_method_results.rds")
  rad_results  <- readRDS("results/radiation_method_results.rds")
  sim_results  <- readRDS("results/simulation/simulation_all_results.rds")
  cll_metadata <- readRDS("results/cll_metadata.rds")

  # Pairwise metrics
  scores_metrics <- compute_pairwise_metrics(
    cll_results$aligned_scores,
    lapply(cll_results$aligned_scores, threshold_tau3)
  )

  # Figure 3C
  p3c <- plot_pairwise_heatmap(scores_metrics, metric = "pearson",
                                title = "Factor 1 scores — Pearson correlation")
  ggsave(file.path(FIGURE_DIR, "fig3C_scores_pearson.pdf"), p3c, width = 5, height = 4)

  # Figure 3D
  p3d <- plot_pairwise_heatmap(scores_metrics, metric = "cosine",
                                title = "Factor 1 scores — Cosine similarity")
  ggsave(file.path(FIGURE_DIR, "fig3D_scores_cosine.pdf"), p3d, width = 5, height = 4)

  # Figure 13 — simulation
  p13 <- plot_simulation_jaccard(sim_results)
  ggsave(file.path(FIGURE_DIR, "fig13_simulation_jaccard.pdf"), p13, width = 10, height = 4)

  message("All figures saved to: ", FIGURE_DIR)
}

if (!interactive()) generate_all_figures()
