## --------- 1) Sign alignment to a reference column ---------
align_to_reference <- function(df, ref_col, target_cols = NULL,
                               cor_method = "pearson") {
  stopifnot(ref_col %in% names(df))
  num_ok <- vapply(df, is.numeric, logical(1))
  if (is.null(target_cols)) {
    target_cols <- setdiff(names(df)[num_ok], ref_col)
  } else {
    target_cols <- intersect(target_cols, names(df)[num_ok])
  }
  r_ref <- df[[ref_col]]
  flips <- setNames(logical(length(target_cols)), target_cols)
  cors  <- setNames(rep(NA_real_, length(target_cols)), target_cols)
  
  for (nm in target_cols) {
    ok   <- is.finite(r_ref) & is.finite(df[[nm]])
    if (sum(ok) >= 3) {
      r <- suppressWarnings(stats::cor(r_ref[ok], df[[nm]][ok], method = cor_method))
      cors[nm] <- r
      if (is.finite(r) && r < 0) {
        df[[nm]] <- -df[[nm]]
        flips[nm] <- TRUE
      }
    }
  }
  list(data = df, flipped = flips, corr = cors)
}

## --------- 2) Your existing cutoff + index + Jaccard (unchanged) ----------
varphi.data <- function(x, probs = 0.60, use_abs = FALSE) {
  x <- as.numeric(x); if (use_abs) x <- abs(x)
  stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
}

calculate_threshold_and_index <- function(data, probs = 0.60, use_abs = FALSE) {
  out <- data
  num_cols <- names(out)[vapply(out, is.numeric, logical(1))]
  for (nm in num_cols) {
    thr <- varphi.data(out[[nm]], probs = probs, use_abs = use_abs)
    out[[paste0(nm, "_index")]] <- as.integer(out[[nm]] >= thr)
  }
  out
}

jaccard_active <- function(x, y) {
  x <- as.integer(x != 0); y <- as.integer(y != 0)
  ok <- !is.na(x) & !is.na(y)
  active <- ok & ((x + y) >= 1L)
  if (!any(active)) return(NA_real_)
  inter <- sum(x[active] & y[active])
  denom <- sum(active)
  inter / denom
}

jaccard_similarity_table <- function(data) {
  idx_cols <- grep("_index$", names(data), value = TRUE)
  out <- list()
  k <- 0
  for (i in 1:(length(idx_cols) - 1)) for (j in (i + 1):length(idx_cols)) {
    c1 <- idx_cols[i]; c2 <- idx_cols[j]
    x <- data[[c1]]; y <- data[[c2]]
    sc <- sum(x == 1 & y == 1, na.rm = TRUE)
    dc <- sum(x != y & ((x + y) >= 1L), na.rm = TRUE)
    s1 <- sum(x == 1, na.rm = TRUE); s2 <- sum(y == 1, na.rm = TRUE)
    ji <- round(jaccard_active(x, y), 2)
    k <- k + 1
    out[[k]] <- data.frame(
      Column1 = c1, Column2 = c2, ji = ji,
      similar_count = sc, dissimilar_count = dc,
      sum_column1 = s1, sum_column2 = s2, stringsAsFactors = FALSE
    )
  }
  if (length(out)) do.call(rbind, out) else
    data.frame(Column1=character(), Column2=character(), ji=numeric(),
               similar_count=integer(), dissimilar_count=integer(),
               sum_column1=integer(), sum_column2=integer())
}


#====================================================
# Jaccard Index: Scores
#====================================================
scale_score_rad_tbl = data.frame(
  sample = scale_data_scores_cor$sample,
  F1_FABIA = (-1)*scale_data_scores_cor$`Score (FABIA)`,
  F1_MOFA = scale_data_scores_cor$`Score (MOFA)`,
  F1_MFA = (-1)*scale_data_scores_cor$`Score (MFA)`,
  F1_GFA = (-1)*scale_data_scores_cor$`Score (GFA)`
)

scale_score_rad_tbl = data.frame(
  sample = data_scores_cor$sample,
  F1_FABIA = (-1)*data_scores_cor$score_FABIA1,
  F1_MOFA = data_scores_cor$score_MOFA1,
  F1_MFA = -1*data_scores_cor$score_MFA1,
  F1_GFA = -1*data_scores_cor$score_GFA1
)

# flip all to match MOFA direction first:
score_aligned_rad <- align_to_reference(
  scale_score_rad_tbl, ref_col = "F1_FABIA",
  target_cols = c("F1_MOFA","F1_MFA","F1_GFA")
)$data

score_rad_tbl_index <- calculate_threshold_and_index(scale_score_rad_tbl, probs = 0.60, use_abs = FALSE)
jac_score_rad <- jaccard_similarity_table(score_rad_tbl_index)
jac_score_rad

# Upset plot - DNA methylation
xscore_rad <- list(
  A = which(score_rad_tbl_index$F1_FABIA_index == 1), 
  B = which(score_rad_tbl_index$F1_GFA_index == 1),
  C = which(score_rad_tbl_index$F1_MOFA_index == 1),
  D = which(score_rad_tbl_index$F1_MFA_index == 1)
)

m <- build_comb_mat(xscore_rad)                    # default mode = "distinct"
m_sz <- filter_by_size(m, min_size = 2)    # optional: keep intersections size >= 2
m_d2 <- filter_by_degree(m, degree = 2)    # optional: only 2-way overlaps

# 1) Publication style (counts + numbers)
plot_upset_pub(m, title = "UpSet — distinct mode")

# 2) Color by degree
plot_upset_by_degree(m, title = "")

#====================================================
# Jaccard Index: DNA methylation
#====================================================
mrna_tbl = data.frame(
  feature =df_mRNA_weights_tbl$feature,
  F1_FABIA = (-1)*df_mRNA_weights_tbl$loading_FABIA1,
  F1_MOFA = df_mRNA_weights_tbl$loading_MOFA1,
  F1_MFA = (-1)*df_mRNA_weights_tbl$loading_MFA1,
  F1_GFA = (-1)*df_mRNA_weights_tbl$loading_GFA1
)

# flip all to match MOFA direction first:
mrna_aligned <- align_to_reference(
  mrna_tbl, ref_col = "F1_FABIA",
  target_cols = c("F1_MOFA","F1_MFA","F1_GFA")
)$data

# now threshold on the (signed) aligned scores (no abs):
mrna_tbl_index <- calculate_threshold_and_index(mrna_aligned, probs = 0.60, use_abs = FALSE)
jac_mrna <- jaccard_similarity_table(mrna_tbl_index)
jac_mrna

# Upset plot - DNA methylation
xmrna <- list(
  A = which(mrna_tbl_index$F1_FABIA_index == 1), 
  B = which(mrna_tbl_index$F1_GFA_index == 1),
  C = which(mrna_tbl_index$F1_MOFA_index == 1),
  D = which(mrna_tbl_index$F1_MFA_index == 1)
)

m <- build_comb_mat(xmrna)                    # default mode = "distinct"
m_sz <- filter_by_size(m, min_size = 2)    # optional: keep intersections size >= 2
m_d2 <- filter_by_degree(m, degree = 2)    # optional: only 2-way overlaps

# 1) Publication style (counts + numbers)
plot_upset_pub(m, title = "UpSet — distinct mode")

# 2) Color by degree
plot_upset_by_degree(m, title = "")

#====================================================
# Jaccard Index: Drugs profiles
#====================================================
proteins_tbl = data.frame(
  feature = df_Proteomics_tbl$feature,
  F1_FABIA = (-1)*df_Proteomics_tbl$loading_FABIA1,
  F1_MOFA = df_Proteomics_tbl$loading_MOFA1,
  F1_MFA = (-1)*df_Proteomics_tbl$loading_MFA1,
  F1_GFA = (-1)*df_Proteomics_tbl$loading_GFA1
)

# flip all to match MOFA direction first:
proteins_aligned <- align_to_reference(
  proteins_tbl, ref_col = "F1_FABIA",
  target_cols = c("F1_MOFA","F1_MFA","F1_GFA")
)$data

# now threshold on the (signed) aligned scores (no abs):
proteins_tbl_index <- calculate_threshold_and_index(proteins_aligned, probs = 0.60, use_abs = FALSE)
jac_proteins <- jaccard_similarity_table(proteins_tbl_index)
jac_proteins

# Upset plot - Drugs Profiles
xproteins <- list(
  A = which(proteins_tbl_index$F1_FABIA_index == 1), 
  B = which(proteins_tbl_index$F1_GFA_index == 1),
  C = which(proteins_tbl_index$F1_MOFA_index == 1),
  D = which(proteins_tbl_index$F1_MFA_index == 1)
)

m <- build_comb_mat(xproteins)                    # default mode = "distinct"
m_sz <- filter_by_size(m, min_size = 2)    # optional: keep intersections size >= 2
m_d2 <- filter_by_degree(m, degree = 2)    # optional: only 2-way overlaps

# 1) Publication style (counts + numbers)
plot_upset_pub(m, title = "UpSet — distinct mode")

# 2) Color by degree
plot_upset_by_degree(m, title = "")

#===============================================================================
# Overal plots jaccard index
#===============================================================================

# ----- sweep helpers -----

# cutoffs (one row per numeric column per p)
cutoffs_over_probs <- function(df, probs_grid, use_abs = FALSE) {
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  rows <- list(); k <- 0
  for (p in probs_grid) {
    thr <- sapply(num_cols, function(nm) varphi.data(df[[nm]], probs = p, use_abs = use_abs))
    rows[[k <- k + 1]] <- data.frame(
      probs   = rep(p, length(num_cols)),
      column  = num_cols,
      cutoff  = as.numeric(thr),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# Jaccard tables (pairwise over *_index columns) for each p, stacked
jaccard_over_probs <- function(df, probs_grid, use_abs = FALSE) {
  out <- lapply(probs_grid, function(p) {
    idx_df <- calculate_threshold_and_index(df, probs = p, use_abs = use_abs)
    ji_tbl <- jaccard_similarity_table(idx_df)
    if (nrow(ji_tbl)) {
      ji_tbl$probs <- p
      ji_tbl
    } else NULL
  })
  do.call(rbind, Filter(Negate(is.null), out))
}

library(dplyr)
library(ggplot2)

# pull method name from something like "F1_FABIA_index" -> "FABIA"
.extract_method <- function(x) {
  x <- sub("_index$", "", x)
  x <- sub("^F[0-9]+_", "", x)
  x
}

plot_jaccard_by_prob <- function(jaccard_df, facet = FALSE, title = "Jaccard vs quantile cutoff") {
  df <- jaccard_df %>%
    mutate(
      m1 = .extract_method(Column1),
      m2 = .extract_method(Column2)
    ) %>%
    rowwise() %>%
    mutate(combo = paste(sort(c(m1, m2)), collapse = "–")) %>%  # order-insensitive label
    ungroup() %>%
    mutate(
      probs = as.numeric(probs),
      ji    = as.numeric(ji)
    ) %>%
    group_by(combo, probs) %>%                      # just in case there are duplicates per (combo, probs)
    summarise(ji = mean(ji, na.rm = TRUE), .groups = "drop")
  
  base <- ggplot(df, aes(x = probs, y = ji, color = combo, group = combo)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.05))) +
    scale_x_continuous(breaks = sort(unique(df$probs))) +
    labs(x = "Quantile cutoff (probs.)", y = "Jaccard index", color = "Combination", title = title) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          legend.position = if (facet) "none" else "right",
          plot.title = element_text(face = "bold"))
  
  if (facet) {
    base + facet_wrap(~ combo, ncol = 3)
  } else {
    base
  }
}

#===============================================================================
# Overal plots jaccard index: Scores
#===============================================================================
# grid: 0.50, 0.55, ..., 0.95, 0.975
probs_grid <- c(seq(0.50, 0.95, by = 0.05), 0.95)

# 1) thresholds used at each p (long table)
score_cutoffs_df_rad <- cutoffs_over_probs(score_aligned_rad, probs_grid, use_abs = FALSE)

# 2) Jaccard at each p (long table; one row per column pair × p)
score_jaccard_df_rad <- jaccard_over_probs(score_aligned_rad, probs_grid, use_abs = FALSE)

# 3) Single panel with six colored lines
p1_rad <- plot_jaccard_by_prob(score_jaccard_df_rad, facet = FALSE,
                           title = " ")
print(p1_rad)

# Or faceted (one small plot per combination)
p2_rad <- plot_jaccard_by_prob(score_jaccard_df_rad, facet = TRUE,
                           title = "Jaccard across cutoffs (faceted)")
print(p2_rad)

#===============================================================================
# Overal plots jaccard index: DNA Methylation
#===============================================================================
# grid: 0.50, 0.55, ..., 0.95, 0.975
probs_grid <- c(seq(0.50, 0.95, by = 0.05), 0.95)

# 1) thresholds used at each p (long table)
mrna_cutoffs_df <- cutoffs_over_probs(mrna_aligned, probs_grid, use_abs = FALSE)

# 2) Jaccard at each p (long table; one row per column pair × p)
mrna_jaccard_df <- jaccard_over_probs(mrna_aligned, probs_grid, use_abs = FALSE)

# 3) Single panel with six colored lines
p_mrna <- plot_jaccard_by_prob(mrna_jaccard_df, facet = FALSE,
                           title = " ")
print(p_mrna)

#===============================================================================
# Overal plots jaccard index: Drugs Profiles
#===============================================================================
# grid: 0.50, 0.55, ..., 0.95, 0.975
probs_grid <- c(seq(0.50, 0.95, by = 0.05), 0.95)

# 1) thresholds used at each p (long table)
proteins_cutoffs_df <- cutoffs_over_probs(proteins_aligned, probs_grid, use_abs = FALSE)

# 2) Jaccard at each p (long table; one row per column pair × p)
proteins_jaccard_df <- jaccard_over_probs(proteins_aligned, probs_grid, use_abs = FALSE)

# 3) Single panel with six colored lines
p_proteins <- plot_jaccard_by_prob(proteins_jaccard_df, facet = FALSE,
                                 title = " ")
print(p_proteins)
