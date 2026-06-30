features_omic1 <- paste0("omic1_feature", seq_len(length(true_load_omic1_FACTOR1)))
features_omic2 <- paste0("omic2_feature", seq_len(length(true_load_omic2_FACTOR1)))

features_omic1
features_omic2

# Construct factor 1 loading data frame
df_factor1_loading <- data.frame(
  feature = c(features_omic1, features_omic2),
  beta1   = c(true_load_omic1_FACTOR1, rep(NA, length(true_load_omic2_FACTOR1))),
  delta1  = c(rep(NA, length(true_load_omic1_FACTOR1)), true_load_omic2_FACTOR1),
  true_loading_F1 = c(true_load_omic1_FACTOR1, true_load_omic2_FACTOR1)
)

plot(df_factor1_loading$true_loading_F1)
plot(df_factor1_loading$beta1)
plot(df_factor1_loading$delta1)

# Construct factor 2 loading data frame
df_factor2_loading <- data.frame(
  feature = c(features_omic1, features_omic2),
  beta2   = c(true_load_omic1_FACTOR2, rep(NA, length(true_load_omic2_FACTOR2))),
  delta2  = c(rep(NA, length(true_load_omic1_FACTOR2)), true_load_omic2_FACTOR2),
  true_loading_F2 = c(true_load_omic1_FACTOR2, true_load_omic2_FACTOR2)
)


# TRUE SCORES
samples <- paste0("sample_", seq_len(length(true_score_FACTOR1)))

# Construct factor 1 loading data frame
df_factor_score <- data.frame(
  sample = samples,
  alpha1 = c(true_score_FACTOR1, rep(NA, length(true_score_FACTOR1))),
  alpha2  = c(true_score_FACTOR2, rep(NA, length(true_score_FACTOR2)))
)
df_factor_score <- na.omit(df_factor_score)

plot(df_factor_score$alpha1)
plot(df_factor_score$alpha2)

data.x <- as.numeric(simulation_result[["dataset_output"]][["variance_10_iteration_1"]][["gfa_result"]][["scores"]][["original_scores"]][["score_GFA2"]])
sample <- paste0("sample_", seq_len(length(true_score_FACTOR1)))

full_data <- data.frame(sample, data.x)

full_data <- merge(full_data, df_factor_score, by = "sample")

plot(full_data$data.x, full_data$alpha2

plot(full_data$data.x)
plot(result_score$score_MFA2)
plot(result_load_all$loading_GFA1)
plot(result_load_all$loading_GFA2)
plot(result_score$score_MFA1)
plot(result_load_all$loading_GFA1)
data.x <- as.numeric(simulation_result[["dataset_output"]][["variance_10_iteration_1"]][["gfa_result"]][["scores"]][["original_scores"]][["score_GFA2"]])
sample <- paste0("sample_", seq_len(length(true_score_FACTOR1)))
