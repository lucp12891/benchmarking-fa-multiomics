#===========================================
# Libraries
#===========================================
library(fabia); library(MOFA2); library(MOFAdata);library(MOFAdata);library(data.table); 
library(gridExtra); library(ggplot2); library(tidyverse); library(dplyr); library(GFA);
library(dplyr); library(reshape2); library(patchwork); library(grid)
#===========================================
# Create Benchmark Datasets
#===========================================
set.seed(21)
benchmark_cll_test <- run_factor_pipeline(
  list(drugs = as.matrix(drugs), methylation = as.matrix(methylation)),
  methods = c("fabia","mofa","gfa","mfa"),
  k = 2,
  scale_features = TRUE)

# ============================================================
# merged <- merge_scores_with_metadata(benchmark_cll, CLL_metadata, sample_col = "sample")
# ============================================================
merged <- merge_scores_with_metadata(benchmark_cll_test, CLL_metadata, sample_col = "sample")
# Keep original matrices; add scores_with_meta to each method
benchmark_cll_test <- attach_scores_metadata(benchmark_cll_test, CLL_metadata, sample_col = "sample")

#============================================
# Create score scatter plot by IGHV
#============================================
p_none <- plot_scores_scatter_grid(
  benchmark_cll_test, factor_idx = 1, group_col = "IGHV",
  methods = c("MOFA","FABIA","MFA","GFA"),
  panel_tags = NULL,
  tag_position = "none"
)
print(p_none)

#============================================
# Create scores dataframes
#============================================

# Scores
fabia_score_df <- data.frame(benchmark_cll_test$FABIA$scores); fabia_score_df$sample <- rownames(fabia_score_df)
# rename columns that start with "F"
colnames(fabia_score_df) <- sub("^F(\\d+)$", "F\\1 (FABIA)", colnames(fabia_score_df))

mofa_score_df <- data.frame(benchmark_cll_test$MOFA$scores); mofa_score_df$sample <- rownames(mofa_score_df)
# rename columns that start with "F"
colnames(mofa_score_df) <- sub("^F(\\d+)$", "F\\1 (MOFA)", colnames(mofa_score_df))

mfa_score_df <- data.frame(benchmark_cll_test$MFA$scores); mfa_score_df$sample <- rownames(mfa_score_df)
# rename columns that start with "F"
colnames(mfa_score_df) <- sub("^F(\\d+)$", "F\\1 (MFA)", colnames(mfa_score_df))

gfa_score_df <- data.frame(benchmark_cll_test$GFA$scores); gfa_score_df$sample <- rownames(gfa_score_df)
# rename columns that start with "F"
colnames(gfa_score_df) <- sub("^F(\\d+)$", "F\\1 (GFA)", colnames(gfa_score_df))

names(fabia_score_df)
names(mofa_score_df)
names(mfa_score_df)
names(gfa_score_df)

#============================================
# Create loading dataframes
#============================================
# LOADINGS
fabia_loading_df <- data.frame(benchmark_cll_test$FABIA$loading); fabia_loading_df$feature <- rownames(fabia_loading_df)
# rename columns that start with "F"
colnames(fabia_loading_df) <- sub("^F(\\d+)$", "F\\1 (FABIA)", colnames(fabia_loading_df))

mofa_loading_df <- data.frame(benchmark_cll_test$MOFA$loading); mofa_loading_df$feature <- rownames(mofa_loading_df)
# rename columns that start with "F"
colnames(mofa_loading_df) <- sub("^F(\\d+)$", "F\\1 (MOFA)", colnames(mofa_loading_df))

mfa_loading_df <- data.frame(benchmark_cll_test$MFA$loading); mfa_loading_df$feature <- rownames(mfa_loading_df)
# rename columns that start with "F"
colnames(mfa_loading_df) <- sub("^F(\\d+)$", "F\\1 (MFA)", colnames(mfa_loading_df))

gfa_loading_df <- data.frame(benchmark_cll_test$GFA$loading); gfa_loading_df$feature <- rownames(gfa_loading_df)
# rename columns that start with "F"
colnames(gfa_loading_df) <- sub("^F(\\d+)$", "F\\1 (GFA)", colnames(gfa_loading_df))

names(fabia_loading_df)
names(mofa_loading_df)
names(mfa_loading_df)
names(gfa_loading_df)


# Helper function: split loading df by prefix
split_by_prefix <- function(df, prefix_dr = "^drugs", prefix_met = "^methylation") {
  df_dr <- df[grepl(prefix_dr, df$feature), , drop = FALSE]
  df_met <- df[grepl(prefix_met, df$feature), , drop = FALSE]
  list(drugs = df_dr, methyl = df_met)
}

# Apply to each method
fabia_split <- split_by_prefix(fabia_loading_df)
mofa_split  <- split_by_prefix(mofa_loading_df)
mfa_split   <- split_by_prefix(mfa_loading_df)
gfa_split   <- split_by_prefix(gfa_loading_df)

#### ============================================================
#### 1) Ground truth extraction (scores + per-omic loadings)
####    - Scores rows:  sample_1, sample_2, ...
####    - Loadings rows: "rna::omic1_feature_i", "prot::omic2_feature_i", ...
####    - Per-omic data.frames also carry a `feature` column for merges.
#### ============================================================

# ============================================================
# Robust benchmark for factor methods with differing #columns
# ============================================================
bench_cll_test <- benchmark_factor_methods(benchmark_cll_test)
