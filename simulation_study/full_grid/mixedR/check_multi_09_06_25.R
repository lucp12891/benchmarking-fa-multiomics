setwd("/user/leuven/364/vsc36498")
set.seed(60667) 

# Create a specific environment to hold the variables
global_env <- new.env() # shared_env

# Define the main simulation function
simulate_MOmicsR <- function(n_features_one, 
                             n_features_two, 
                             n_samples, 
                             var_sigma, 
                             num_biclusters = 2, 
                             num_iterations) {
  # Assign parameters to global environment
  global_env$n_features_one <- n_features_one
  global_env$n_features_two <- n_features_two
  
  # Initialize result lists
  jaccard_results <- list()
  jaccard_comparison_results <- list()
  per_measures_results <- list()
  dataset_output <- list()
  dataset_output_corrected <- list()
  datasets <- list()
  
  iter <- num_iterations
  # Loop through each value of sigma
  for (sigma in var_sigma) {
    message(paste("Processing sigma =", sigma))
    
# Load packages and suppress message
packages <- c("ggplot2", "dplyr", "tidyr", "viridis", "fabia", "FactoMineR", "readr", "readxl", "tidyverse",
              "stringr", "basilisk", "MOFA2", "data.table", "GFA", "pROC", "mclust", "zoo")

for (pkg in packages) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# set seed

#----------------------------- simulate data ----------------------------
for (i in 1:iter) {
  current_seed <- 695#sample(1e6, 1)
  set.seed(current_seed)
  
simulated_data <- multiple_factor(
  n_features_one = 4000,#n_features_one,
  n_features_two = 3000,#n_features_two,
  n_samples = n_samples,
  sigmas = sigma,
  iterations = i,
  n_factors = num_biclusters
)


  
# extract data
current_data <- simulated_data[[paste0("iteration_", i)]]$concatenated_datasets[[1]]
true_score_FACTOR1 <- simulated_data[[paste0("iteration_", i)]][["list_alphas"]][["alpha1"]]
true_score_FACTOR2 <- simulated_data[[paste0("iteration_", i)]][["list_alphas"]][["alpha2"]]
true_load_omic1_FACTOR1 <- simulated_data[[paste0("iteration_", i)]][["list_betas"]][["beta1"]]
true_load_omic1_FACTOR2 <- simulated_data[[paste0("iteration_", i)]][["list_betas"]][["beta2"]]
true_load_omic2_FACTOR1 <- simulated_data[[paste0("iteration_", i)]][["list_deltas"]][["delta1"]]
true_load_omic2_FACTOR2 <- rnorm(n_features_two, 0, 0.05)

indices_features.OMIC1.A <- simulated_data[[paste0("iteration_", i)]]$indices_features.1[[1]]
indices_features.OMIC1.B <- simulated_data[[paste0("iteration_", i)]]$indices_features.1[[2]]
indices_features.OMIC2.A <- simulated_data[[paste0("iteration_", i)]]$indices_features.2[[1]]
indices_samples.1A <- simulated_data[[paste0("iteration_", i)]]$indices_samples[[1]]
indices_samples.2B <- simulated_data[[paste0("iteration_", i)]]$indices_samples[[2]]

print(paste("Processing iteration:", i))

# Run factorization methods
print(paste("FABIA MODELLING IN PROGRESS..."))
#run_with_seed({
fabia_result <- func_fabia(current_data, num_biclusters)
#plot(fabia_result[["scores"]][["scores"]][["score_FABIA1"]])
#plot(fabia_result[["scores"]][["scores"]][["score_FABIA2"]])
print(paste("MOFA MODELLING IN PROGRESS..."))
mofa_result <- func_mofa(current_data, num_biclusters)
#plot(mofa_result[["scores"]][["scores"]][["score_MOFA1"]])
#plot(mofa_result[["scores"]][["scores"]][["score_MOFA2"]])
print(paste("MFA MODELLING IN PROGRESS...")) 
mfa_result <- func_mfa(current_data, num_biclusters)
#plot(mfa_result[["scores"]][["factor_scores"]][["score_MFA1"]])
#plot(mfa_result[["scores"]][["factor_scores"]][["score_MFA2"]])
print(paste("GFA MODELLING IN PROGRESS..."))
gfa_result <- func_gfa(current_data, num_biclusters)
#plot(gfa_result[["scores"]][["factor_scores"]][["score_GFA1"]])
#plot(mfa_result[["scores"]][["factor_scores"]][["score_MFA2"]])
#}, seed = 60667)

random.data <- data.frame(t(current_data)); random.data$feature = rownames(random.data) 
random.data$ID <- ifelse(grepl("omic1_", random.data$feature),
                         sub(".*omic1_feature_(\\d+)", "\\1", random.data$feature),
                         ifelse(grepl("omic2_", random.data$feature),
                                sub(".*omic2_feature_(\\d+)", "\\1", random.data$feature),NA))

random.data$ID <- as.numeric(random.data$ID)
# Create dataview i.e that will help separate the two merged dataset
random.data <- random.data %>%
  mutate(dataview = ifelse(grepl("omic1_", feature), "omic.one", ifelse(grepl("omic2_", feature), "omic.two", NA)))

# Filtering the dataset
random.data.a <- random.data %>% filter(dataview == "omic.one")
random.data.b <- random.data %>% filter(dataview == "omic.two")

# Assigning 'indices_features.1' to 'random.data.a'
signal_features.omic1_F1 = simulated_data[[paste0("iteration_", i)]]$indices_features.1[1]
signal_features.omic1_F2 = simulated_data[[paste0("iteration_", i)]]$indices_features.1[2]
signal_features.omic2_F1 = simulated_data[[paste0("iteration_", i)]]$indices_features.2

# FACTOR 1 LOADING
# Omic One
signal_features.omic1_F1_vector <- unlist(signal_features.omic1_F1)  # Convert to atomic vector
signal_features.omic1_F1_sorted <- sort(signal_features.omic1_F1_vector)

# Omic Two
signal_features.omic2_F1_vector <- unlist(signal_features.omic2_F1)  # Convert to atomic vector
signal_features.omic2_F1_sorted <- sort(signal_features.omic2_F1_vector)

in_range.omic1_F1 <- random.data.a$ID %in% signal_features.omic1_F1_sorted

simulated_features_omic1_F1 <- data.frame(feature = random.data.a$feature, signal_a = in_range.omic1_F1)

in_range.omic2_F1 <- random.data.b$ID %in% signal_features.omic2_F1_sorted

simulated_features_omic2_F1 <- data.frame(feature = random.data.b$feature, signal_a = in_range.omic2_F1)

simulated_features_sigma_F1 = rbind(simulated_features_omic1_F1, simulated_features_omic2_F1)

# important signals
#simulated_features_sigma_F1
#simulated_features_omic1_F1
#simulated_features_omic2_F1

# FACTOR 2 LOADING
# Omic One
signal_features.omic1_F2_vector <- unlist(signal_features.omic1_F2)  # Convert to atomic vector
signal_features.omic1_F2_sorted <- sort(signal_features.omic1_F2_vector)

#signal_features.omic2_F2_vector <- unlist(signal_features.omic2_F2)  # Convert to atomic vector
#signal_features.omic2_F2_sorted <- sort(signal_features.omic2_F2_vector)

in_range.omic1_F2 <- random.data.a$ID %in% signal_features.omic1_F2_sorted

simulated_features_omic1_F2 <- data.frame(feature = random.data.a$feature, signal_b = in_range.omic1_F2)
simulated_features_omic2_F2 <- data.frame(feature = random.data.b$feature, signal_b = FALSE) # There is no signal for F2 OMIC2

simulated_features_sigma_b = rbind(simulated_features_omic1_F2, simulated_features_omic2_F2)

# FACTOR 1 SCORES
indices_samples.F1 = simulated_data[[paste0("iteration_", i)]]$indices_samples[1]
indices_samples.F2 = simulated_data[[paste0("iteration_", i)]]$indices_samples[2]

indices_samples.F1_vector <- unlist(indices_samples.F1)  # Convert to atomic vector
indices_samples.F1_sorted <- sort(indices_samples.F1_vector)
indices_samples.F2_vector <- unlist(indices_samples.F2)  # Convert to atomic vector
indices_samples.F2_sorted <- sort(indices_samples.F2_vector)

# samples
# Assuming indices_samples is a vector
indices_samples1 <- paste("sample_", unlist(simulated_data[[paste0("iteration_", i)]]$indices_samples[1]), sep = "")
indices_samples2 <- paste("sample_", unlist(simulated_data[[paste0("iteration_", i)]]$indices_samples[2]), sep = "")


in_range_sample_a <- rownames(current_data) %in% indices_samples1
simulated_samples_sigma_a <- data.frame(sample = rownames(current_data),
                                        signal_a = in_range_sample_a)

in_range_sample_b <- rownames(current_data) %in% indices_samples2
simulated_samples_sigma_b <- data.frame(sample = rownames(current_data),
                                        signal_b = in_range_sample_b)

# Extract loadings
# OMIC ONE
fabia_df <- extract_id_from_feature(fabia_result[["weights"]][["omic.one_weights"]])
fabia_aligned <- align_loading_columns_by_truth(fabia_df, true_load_omic1_FACTOR1)
mofa_df  <- extract_id_from_feature(mofa_result[["weights"]][["omic.one_weights"]])
mofa_aligned <- align_loading_columns_by_truth(mofa_df, true_load_omic1_FACTOR1)
mfa_df   <- extract_id_from_feature(mfa_result[["weights"]][["omic.one_weights"]])
#mfa_aligned <- align_loading_columns_by_truth(mfa_df, true_load_omic1_FACTOR1)
gfa_df   <- extract_id_from_feature(gfa_result[["weights"]][["omic.one_weights"]])
#gfa_aligned <- align_loading_columns_by_truth(gfa_df, true_load_omic1_FACTOR1)

fabia_sub <- fabia_aligned[, c("ID", "loading_FABIA1")]
mofa_sub  <- mofa_aligned[, c("ID", "loading_MOFA1", "feature", "dataview")]
mfa_sub   <- mfa_df[, c("ID", "loading_MFA1")]
gfa_sub   <- gfa_df[, c("ID", "loading_GFA1", "signal_a", "signal_b")]

merged_loadings_omic1_F1 <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE),
                          list(fabia_sub, mofa_sub, mfa_sub, gfa_sub))

#Convert all numeric columns to absolute values
abs.merged_loadings_omic1_F1 <- merged_loadings_omic1_F1
abs.merged_loadings_omic1_F1$loading_FABIA1= abs(abs.merged_loadings_omic1_F1$loading_FABIA1)
abs.merged_loadings_omic1_F1$loading_MOFA1 = abs(abs.merged_loadings_omic1_F1$loading_MOFA1)
abs.merged_loadings_omic1_F1$loading_MFA1 = abs(abs.merged_loadings_omic1_F1$loading_MFA1)
abs.merged_loadings_omic1_F1$loading_GFA1 = abs(abs.merged_loadings_omic1_F1$loading_GFA1)

fabia_sub <- fabia_aligned[, c("ID", "loading_FABIA2")]
mofa_sub  <- mofa_aligned[, c("ID", "loading_MOFA2", "feature", "dataview")]
mfa_sub   <- mfa_df[, c("ID", "loading_MFA2")]
gfa_sub   <- gfa_df[, c("ID", "loading_GFA2", "signal_a", "signal_b")]

merged_loadings_omic1_F2 <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE),
                                list(fabia_sub, mofa_sub, mfa_sub, gfa_sub))

#Convert all numeric columns to absolute values
abs.merged_loadings_omic1_F2 <- merged_loadings_omic1_F2 
abs.merged_loadings_omic1_F2$loading_FABIA2= abs(abs.merged_loadings_omic1_F2$loading_FABIA2)
abs.merged_loadings_omic1_F2$loading_MOFA2 = abs(abs.merged_loadings_omic1_F2$loading_MOFA2)
abs.merged_loadings_omic1_F2$loading_MFA2 = abs(abs.merged_loadings_omic1_F2$loading_MFA2)
abs.merged_loadings_omic1_F2$loading_GFA2 = abs(abs.merged_loadings_omic1_F2$loading_GFA2)

# OMIC TWO
fabia_df_omic2 <- extract_id_from_feature(fabia_result[["weights"]][["omic.two_weights"]])
fabia_aligned <- align_loading_columns_by_truth(fabia_df_omic2, true_load_omic2_FACTOR1)
mofa_df_omic2  <- extract_id_from_feature(mofa_result[["weights"]][["omic.two_weights"]])
mofa_aligned <- align_loading_columns_by_truth(mofa_df_omic2, true_load_omic2_FACTOR1)
mfa_df_omic2   <- extract_id_from_feature(mfa_result[["weights"]][["omic.two_weights"]])
#mfa_aligned <- align_loading_columns_by_truth(mfa_df_omic2, true_load_omic2_FACTOR1)
gfa_df_omic2   <- extract_id_from_feature(gfa_result[["weights"]][["omic.two_weights"]])
#gfa_aligned <- align_loading_columns_by_truth(gfa_df_omic2, true_load_omic2_FACTOR1)

fabia_sub_omic2 <- fabia_aligned[, c("ID", "loading_FABIA1")]
mofa_sub_omic2  <- mofa_aligned[, c("ID", "loading_MOFA1", "feature", "dataview")]
mfa_sub_omic2   <- mfa_df_omic2[, c("ID", "loading_MFA1")]
gfa_sub_omic2   <- gfa_df_omic2[, c("ID", "loading_GFA1", "signal_a")]

merged_loadings_omic2_F1 <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE),
                                list(fabia_sub_omic2, mofa_sub_omic2, mfa_sub_omic2, gfa_sub_omic2))

#Convert all numeric columns to absolute values
abs.merged_loadings_omic2_F1 <- merged_loadings_omic2_F1 
abs.merged_loadings_omic2_F1$loading_FABIA1= abs(abs.merged_loadings_omic2_F1$loading_FABIA1)
abs.merged_loadings_omic2_F1$loading_MOFA1 = abs(abs.merged_loadings_omic2_F1$loading_MOFA1)
abs.merged_loadings_omic2_F1$loading_MFA1 = abs(abs.merged_loadings_omic2_F1$loading_MFA1)
abs.merged_loadings_omic2_F1$loading_GFA1 = abs(abs.merged_loadings_omic2_F1$loading_GFA1)

fabia_sub_omic2 <- fabia_aligned[, c("ID", "loading_FABIA2")]
mofa_sub_omic2  <- mofa_aligned[, c("ID", "loading_MOFA2", "feature", "dataview")]
mfa_sub_omic2   <- mfa_df_omic2[, c("ID", "loading_MFA2")]
gfa_sub_omic2   <- gfa_df_omic2[, c("ID", "loading_GFA2")]

merged_loadings_omic2_F2 <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE),
                                   list(fabia_sub_omic2, mofa_sub_omic2, mfa_sub_omic2, gfa_sub_omic2))

#Convert all numeric columns to absolute values
abs.merged_loadings_omic2_F2 <- merged_loadings_omic2_F2 
abs.merged_loadings_omic2_F2$loading_FABIA2= abs(abs.merged_loadings_omic2_F2$loading_FABIA2)
abs.merged_loadings_omic2_F2$loading_MOFA2 = abs(abs.merged_loadings_omic2_F2$loading_MOFA2)
abs.merged_loadings_omic2_F2$loading_MFA2 = abs(abs.merged_loadings_omic2_F2$loading_MFA2)
abs.merged_loadings_omic2_F2$loading_GFA2 = abs(abs.merged_loadings_omic2_F2$loading_GFA2)

# SCORES
# FACTOR 1
fabia_df <- extract_id_from_sample (fabia_result[["scores"]][["scores"]])
fabia_aligned <- align_scores_columns_by_truth(fabia_df, true_score_FACTOR1)
mofa_df  <- extract_id_from_sample (mofa_result[["scores"]][["scores"]])
mofa_aligned <- align_scores_columns_by_truth(mofa_df, true_score_FACTOR1)
mfa_df   <- extract_id_from_sample (mfa_result[["scores"]][["factor_scores"]])
#mfa_aligned <- align_scores_columns_by_truth(mfa_df, true_score_FACTOR1)
gfa_df   <- extract_id_from_sample (gfa_result[["scores"]][["factor_scores"]])
#gfa_aligned <- align_scores_columns_by_truth(gfa_df, true_score_FACTOR1)

fabia_sub <- fabia_aligned[, c("ID", "score_FABIA1")]
mofa_sub  <- mofa_aligned[, c("ID", "score_MOFA1", "sample")]
mfa_sub   <- mfa_df[, c("ID", "score_MFA1")]
gfa_sub   <- gfa_df[, c("ID", "score_GFA1", "signal_a", "signal_b")]

merged_scores_F1 <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE),
                                   list(fabia_sub, mofa_sub, mfa_sub, gfa_sub))

#Convert all numeric columns to absolute values
abs.merged_scores_F1 <- merged_scores_F1 
abs.merged_scores_F1$score_FABIA1= abs(abs.merged_scores_F1$score_FABIA1)
abs.merged_scores_F1$score_MOFA1 = abs(abs.merged_scores_F1$score_MOFA1)
abs.merged_scores_F1$score_MFA1 = abs(abs.merged_scores_F1$score_MFA1)
abs.merged_scores_F1$score_GFA1 = abs(abs.merged_scores_F1$score_GFA1)

# FACTOR 2
fabia_sub <- fabia_aligned[, c("ID", "score_FABIA2")]
mofa_sub  <- mofa_aligned[, c("ID", "score_MOFA2", "sample")]
mfa_sub   <- mfa_df[, c("ID", "score_MFA2")]
gfa_sub   <- gfa_df[, c("ID", "score_GFA2", "signal_a", "signal_b")]

merged_scores_F2 <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE),
                           list(fabia_sub, mofa_sub, mfa_sub, gfa_sub))

#Convert all numeric columns to absolute values
abs.merged_scores_F2 <- merged_scores_F2 
abs.merged_scores_F2$score_FABIA2= abs(abs.merged_scores_F2$score_FABIA2)
abs.merged_scores_F2$score_MOFA2 = abs(abs.merged_scores_F2$score_MOFA2)
abs.merged_scores_F2$score_MFA2 = abs(abs.merged_scores_F2$score_MFA2)
abs.merged_scores_F2$score_GFA2 = abs(abs.merged_scores_F2$score_GFA2)

# Merge loadings
# Manually remove overlapping columns *before* joining
abs.merged_loadings_omic1_F2_clean <- abs.merged_loadings_omic1_F2 %>%
  select(-any_of(names(abs.merged_loadings_omic1_F1)[names(abs.merged_loadings_omic1_F1) != "feature"]))

# Join with no conflicts
FeatureData.Omic1 <- left_join(abs.merged_loadings_omic1_F1, abs.merged_loadings_omic1_F2_clean, by = "feature")

abs.merged_loadings_omic2_F2_clean <- abs.merged_loadings_omic2_F2 %>%
  select(-any_of(names(abs.merged_loadings_omic2_F1)[names(abs.merged_loadings_omic2_F1) != "feature"]))

# Join with no conflicts
FeatureData.Omic2 <- left_join(abs.merged_loadings_omic2_F1, abs.merged_loadings_omic2_F2_clean, by = "feature")
FeatureData.Omic2$signal_b = FALSE
# Scores
abs.merged_scores_F2_clean <- abs.merged_scores_F2 %>%
  select(-any_of(names(abs.merged_scores_F1)[names(abs.merged_scores_F1) != "sample"]))

# Join with no conflicts
SampleData <- left_join(abs.merged_scores_F1, abs.merged_scores_F2_clean, by = "sample")
# --------------------- varphi_functions ---------------------------
for (varphi in seq_along(varphi_functions)) {
varphi_function <- names(varphi_functions)[[varphi]]

# Apply the function to FeatureData.a.one
## ----------------- FeatureData.a.one ----------------
if (varphi_function == "varphi.five"){
  # Apply the function to Data_Ft.a
  FeatureData.a.one <- FeatureData.Omic1 %>%
    mutate(
      fabia_index_a = ifelse(
        abs(loading_FABIA1) >= do.call(get(varphi_function), list(
          loading_or_score = loading_FABIA1,
          omic = 1, factor = 1, type = "loading"
        )), 1, 0),
      
      mofa_index_a = ifelse(
        abs(loading_MOFA1) >= do.call(get(varphi_function), list(
          loading_or_score = loading_MOFA1,
          omic = 1, factor = 1, type = "loading"
        )), 1, 0),
      
      mfa_index_a = ifelse(
        abs(loading_MFA1) >= do.call(get(varphi_function), list(
          loading_or_score = loading_MFA1,
          omic = 1, factor = 1, type = "loading"
        )), 1, 0),
      
      gfa_index_a = ifelse(
        abs(loading_GFA1) >= do.call(get(varphi_function), list(
          loading_or_score = loading_GFA1,
          omic = 1, factor = 1, type = "loading"
        )), 1, 0),
      
      signal_index_a = ifelse(signal_a == 'TRUE', 1, 0),
      
      fabia_index_b = ifelse(
        abs(loading_FABIA2) >= do.call(get(varphi_function), list(
          loading_or_score = loading_FABIA2,
          omic = 1, factor = 2, type = "loading"
        )), 1, 0),
      
      mofa_index_b = ifelse(
        abs(loading_MOFA2) >= do.call(get(varphi_function), list(
          loading_or_score = loading_MOFA2,
          omic = 1, factor = 2, type = "loading"
        )), 1, 0),
      
      mfa_index_b = ifelse(
        abs(loading_MFA2) >= do.call(get(varphi_function), list(
          loading_or_score = loading_MFA2,
          omic = 1, factor = 2, type = "loading"
        )), 1, 0),
      
      gfa_index_b = ifelse(
        abs(loading_GFA2) >= do.call(get(varphi_function), list(
          loading_or_score = loading_GFA2,
          omic = 1, factor = 2, type = "loading"
        )), 1, 0),
      
      signal_index_b = ifelse(signal_b == 'TRUE', 1, 0)
    )
  
}else{
  # Apply the function to Data_Ft.a
  FeatureData.a.one <- FeatureData.Omic1 %>%
    mutate(
      fabia_index_a = ifelse(abs(loading_FABIA1) >= get(varphi_function)(loading_FABIA1), 1, 0),
      mofa_index_a = ifelse(abs(loading_MOFA1) >= get(varphi_function)(loading_MOFA1), 1, 0),
      mfa_index_a = ifelse(abs(loading_MFA1) >= get(varphi_function)(loading_MFA1), 1, 0),
      gfa_index_a = ifelse(abs(loading_GFA1) >= get(varphi_function)(loading_GFA1), 1, 0),
      signal_index_a = ifelse(signal_a == 'TRUE', 1, 0),
      fabia_index_b = ifelse(abs(loading_FABIA2) >= get(varphi_function)(loading_FABIA2), 1, 0),
      mofa_index_b = ifelse(abs(loading_MOFA2) >= get(varphi_function)(loading_MOFA2), 1, 0),
      mfa_index_b = ifelse(abs(loading_MFA2) >= get(varphi_function)(loading_MFA2), 1, 0),
      gfa_index_b = ifelse(abs(loading_GFA2) >= get(varphi_function)(loading_GFA2), 1, 0),
      signal_index_b = ifelse(signal_b == 'TRUE', 1, 0)
    ) #%>% arrange(as.numeric(ID))
}

FeatureData.a.one <- FeatureData.a.one  %>% arrange(as.numeric(ID))

## ----------------- FeatureData.b.two ----------------
if (varphi_function == "varphi.five"){
  # Apply the function to Data_Ft.a
  FeatureData.b.two <- FeatureData.Omic2 %>%
    mutate(
      fabia_index_a = ifelse(
        abs(loading_FABIA1) >= do.call(get(varphi_function), list(
          loading_or_score = loading_FABIA1,
          omic = 2, factor = 2, type = "loading"
        )), 1, 0),
      
      mofa_index_a = ifelse(
        abs(loading_MOFA1) >= do.call(get(varphi_function), list(
          loading_or_score = loading_MOFA1,
          omic = 2, factor = 2, type = "loading"
        )), 1, 0),
      
      mfa_index_a = ifelse(
        abs(loading_MFA1) >= do.call(get(varphi_function), list(
          loading_or_score = loading_MFA1,
          omic = 2, factor = 2, type = "loading"
        )), 1, 0),
      
      gfa_index_a = ifelse(
        abs(loading_GFA1) >= do.call(get(varphi_function), list(
          loading_or_score = loading_GFA1,
          omic = 2, factor = 2, type = "loading"
        )), 1, 0),
      
      signal_index_a = ifelse(signal_a == 'TRUE', 1, 0),
      
      fabia_index_b = ifelse(
        abs(loading_FABIA2) >= do.call(get(varphi_function), list(
          loading_or_score = loading_FABIA2,
          omic = 1, factor = 2, type = "loading"
        )), 1, 0),
      
      mofa_index_b = ifelse(
        abs(loading_MOFA2) >= do.call(get(varphi_function), list(
          loading_or_score = loading_MOFA2,
          omic = 1, factor = 2, type = "loading"
        )), 1, 0),
      
      mfa_index_b = ifelse(
        abs(loading_MFA2) >= do.call(get(varphi_function), list(
          loading_or_score = loading_MFA2,
          omic = 1, factor = 2, type = "loading"
        )), 1, 0),
      
      gfa_index_b = ifelse(
        abs(loading_GFA2) >= do.call(get(varphi_function), list(
          loading_or_score = loading_GFA2,
          omic = 1, factor = 2, type = "loading"
        )), 1, 0),
      
      signal_index_b = ifelse(signal_b == 'TRUE', 1, 0)
    )
  
}else{
  # Apply the function to Data_Ft.a
  FeatureData.b.two <- FeatureData.Omic2 %>%
    mutate(
      fabia_index_a = ifelse(abs(loading_FABIA1) >= get(varphi_function)(loading_FABIA1), 1, 0),
      mofa_index_a = ifelse(abs(loading_MOFA1) >= get(varphi_function)(loading_MOFA1), 1, 0),
      mfa_index_a = ifelse(abs(loading_MFA1) >= get(varphi_function)(loading_MFA1), 1, 0),
      gfa_index_a = ifelse(abs(loading_GFA1) >= get(varphi_function)(loading_GFA1), 1, 0),
      signal_index_a = ifelse(signal_a == 'TRUE', 1, 0),
      fabia_index_b = ifelse(abs(loading_FABIA2) >= get(varphi_function)(loading_FABIA2), 1, 0),
      mofa_index_b = ifelse(abs(loading_MOFA2) >= get(varphi_function)(loading_MOFA2), 1, 0),
      mfa_index_b = ifelse(abs(loading_MFA2) >= get(varphi_function)(loading_MFA2), 1, 0),
      gfa_index_b = ifelse(abs(loading_GFA2) >= get(varphi_function)(loading_GFA2), 1, 0),
      signal_index_b = ifelse(signal_b == 'TRUE', 1, 0)
    ) #%>% arrange(as.numeric(ID))
}

FeatureData.b.two <- FeatureData.b.two %>% arrange(as.numeric(ID))

## ----------------- SampleData2 ----------------
if (varphi_function == "varphi.five"){
  SampleData2 <- SampleData %>%
    mutate(
      fabia_index_a = ifelse(
        abs(score_FABIA1) >= do.call(get(varphi_function), list(
          loading_or_score = score_FABIA1,
          factor = 1, type = "score"
        )), 1, 0),
      
      mofa_index_a = ifelse(
        abs(score_MOFA1) >= do.call(get(varphi_function), list(
          loading_or_score = score_MOFA1,
          factor = 1, type = "score"
        )), 1, 0),
      
      mfa_index_a = ifelse(
        abs(score_MFA1) >= do.call(get(varphi_function), list(
          loading_or_score = score_MFA1,
          factor = 1, type = "score"
        )), 1, 0),
      
      gfa_index_a = ifelse(
        abs(score_GFA1) >= do.call(get(varphi_function), list(
          loading_or_score = score_GFA1,
          factor = 1, type = "score"
        )), 1, 0),
      
      signal_index_a = ifelse(signal_a == 'TRUE', 1, 0),
      
      fabia_index_b = ifelse(
        abs(score_FABIA2) >= do.call(get(varphi_function), list(
          loading_or_score = score_FABIA2,
          factor = 2, type = "score"
        )), 1, 0),
      
      mofa_index_b = ifelse(
        abs(score_MOFA2) >= do.call(get(varphi_function), list(
          loading_or_score = score_MOFA2,
          factor = 2, type = "score"
        )), 1, 0),
      
      mfa_index_b = ifelse(
        abs(score_MFA2) >= do.call(get(varphi_function), list(
          loading_or_score = score_MFA2,
          factor = 2, type = "score"
        )), 1, 0),
      
      gfa_index_b = ifelse(
        abs(score_GFA2) >= do.call(get(varphi_function), list(
          loading_or_score = score_GFA2,
          factor = 2, type = "score"
        )), 1, 0),
      
      signal_index_b = ifelse(signal_b == 'TRUE', 1, 0)
    )
  
}else{
  # Apply the function to Data_Ft.a
  SampleData2 <- SampleData %>%
    mutate(
      fabia_index_a = ifelse(abs(score_FABIA1) >= get(varphi_function)(score_FABIA1), 1, 0),
      mofa_index_a = ifelse(abs(score_MOFA1) >= get(varphi_function)(score_MOFA1), 1, 0),
      mfa_index_a = ifelse(abs(score_MFA1) >= get(varphi_function)(score_MFA1), 1, 0),
      gfa_index_a = ifelse(abs(score_GFA1) >= get(varphi_function)(score_GFA1), 1, 0),
      signal_index_a = ifelse(signal_a == 'TRUE', 1, 0),
      fabia_index_b = ifelse(abs(score_FABIA2) >= get(varphi_function)(score_FABIA2), 1, 0),
      mofa_index_b = ifelse(abs(score_MOFA2) >= get(varphi_function)(score_MOFA2), 1, 0),
      mfa_index_b = ifelse(abs(score_MFA2) >= get(varphi_function)(score_MFA2), 1, 0),
      gfa_index_b = ifelse(abs(score_GFA2) >= get(varphi_function)(score_GFA2), 1, 0),
      signal_index_b = ifelse(signal_b == 'TRUE', 1, 0)
    ) #%>% arrange(ID.x)
}

SampleData3 <- SampleData2 %>% arrange(as.numeric(ID))

# List of methods and corresponding indices
methods <- c("fabia", "mofa", "mfa", "gfa")

ji_omic_one_a <- sapply(methods, function(method) {
  method_index <- FeatureData.a.one[[paste0(method, "_index_a")]]
  jaccard.index.sim(method_index, FeatureData.a.one$signal_index_a)
}, simplify = TRUE, USE.NAMES = TRUE)

ji_omic_one_b <- sapply(methods, function(method) {
  method_index <- FeatureData.a.one[[paste0(method, "_index_b")]]
  jaccard.index.sim(method_index, FeatureData.a.one$signal_index_b)
}, simplify = TRUE, USE.NAMES = TRUE)

ji_omic_two_a <- sapply(methods, function(method) {
  method_index <- FeatureData.b.two[[paste0(method, "_index_a")]]
  jaccard.index.sim(method_index, FeatureData.b.two$signal_index_a)
}, simplify = TRUE, USE.NAMES = TRUE)

ji_samples_a <- sapply(methods, function(method) {
  method_index <- SampleData3[[paste0(method, "_index_a")]]
  jaccard.index.sim(method_index, SampleData3$signal_index_a)
}, simplify = TRUE, USE.NAMES = TRUE)

ji_omic_one_b <- sapply(methods, function(method) {
  method_index <- FeatureData.a.one[[paste0(method, "_index_b")]]
  jaccard.index.sim(method_index, FeatureData.a.one$signal_index_b)
}, simplify = TRUE, USE.NAMES = TRUE)

ji_omic_two_b <- sapply(methods, function(method) {
  method_index <- FeatureData.b.two[[paste0(method, "_index_b")]]
  jaccard.index.sim(method_index, FeatureData.b.two$signal_index_b)
}, simplify = TRUE, USE.NAMES = TRUE)

ji_samples_b <- sapply(methods, function(method) {
  method_index <- SampleData3[[paste0(method, "_index_b")]]
  jaccard.index.sim(method_index, SampleData3$signal_index_b)
}, simplify = TRUE, USE.NAMES = TRUE)

# Organize the results for the current method
jaccard_loading <- list(ji_omic_one_a=ji_omic_one_a, ji_omic_two_a=ji_omic_two_a, ji_samples_a=ji_samples_a,
                        ji_omic_one_b=ji_omic_one_b, ji_omic_two_b=ji_omic_two_b, ji_samples_b=ji_samples_b)

# Create a unique name for each dataset
#jaccard_results <- list()
dataset_name <- sprintf("variance_%d_iteration_%d_varphi_%d", sigma, i, varphi) 
jaccard_results[[dataset_name]] <- jaccard_loading

# COMPARISONS BY METHODS
methods <- c("fabia", "mofa", "mfa", "gfa")  # Define the list of methods


# JACCARD INDEX FACTOR 1
# OMIC ONE
# Calculate pairwise Jaccard indices dynamically
ji_omic_one_pairs_a <- combn(methods, 2, function(method_pair) {
  method1_index <- FeatureData.a.one[[paste0(method_pair[1], "_index_a")]]
  method2_index <- FeatureData.a.one[[paste0(method_pair[2], "_index_a")]]
  jaccard.index.sim(method1_index, method2_index)
}, simplify = TRUE)

# Assign names to the results
names(ji_omic_one_pairs_a) <- combn(methods, 2, function(method_pair) {
  paste(method_pair[1], method_pair[2], "ji_phi", sep = "_")
})

# Convert to a named vector or data frame for further use
ji_omic_one_pairs_a <- as.data.frame(as.list(ji_omic_one_pairs_a))

# OMIC TWO
# Calculate pairwise Jaccard indices dynamically
ji_omic_two_pairs_a <- combn(methods, 2, function(method_pair) {
  method1_index <- FeatureData.b.two[[paste0(method_pair[1], "_index_a")]]
  method2_index <- FeatureData.b.two[[paste0(method_pair[2], "_index_a")]]
  jaccard.index.sim(method1_index, method2_index)
}, simplify = TRUE)

# Assign names to the results
names(ji_omic_two_pairs_a) <- combn(methods, 2, function(method_pair) {
  paste(method_pair[1], method_pair[2], "ji_phi", sep = "_")
})

# Convert to a named vector or data frame for further use
ji_omic_two_pairs_a <- as.data.frame(as.list(ji_omic_two_pairs_a))

# SAMPLES
# Calculate pairwise Jaccard indices dynamically
ji_smp_pairs_a <- combn(methods, 2, function(method_pair) {
  method1_index <- SampleData3[[paste0(method_pair[1], "_index_a")]]
  method2_index <- SampleData3[[paste0(method_pair[2], "_index_a")]]
  jaccard.index.sim(method1_index, method2_index)
}, simplify = TRUE)

# Assign names to the results
names(ji_smp_pairs_a) <- combn(methods, 2, function(method_pair) {
  paste(method_pair[1], method_pair[2], "ji_phi", sep = "_")
})

# Convert to a named vector or data frame for further use
ji_smp_pairs_a <- as.data.frame(as.list(ji_smp_pairs_a))

# JACCARD INDEX FACTOR 2
# OMIC ONE
# Calculate pairwise Jaccard indices dynamically
ji_omic_one_pairs_b <- combn(methods, 2, function(method_pair) {
  method1_index <- FeatureData.a.one[[paste0(method_pair[1], "_index_b")]]
  method2_index <- FeatureData.a.one[[paste0(method_pair[2], "_index_b")]]
  jaccard.index.sim(method1_index, method2_index)
}, simplify = TRUE)

# Assign names to the results
names(ji_omic_one_pairs_b) <- combn(methods, 2, function(method_pair) {
  paste(method_pair[1], method_pair[2], "ji_phi", sep = "_")
})

# Convert to a named vector or data frame for further use
ji_omic_one_pairs_b <- as.data.frame(as.list(ji_omic_one_pairs_b))

# OMIC TWO
# Calculate pairwise Jaccard indices dynamically
ji_omic_two_pairs_b <- combn(methods, 2, function(method_pair) {
  method1_index <- FeatureData.b.two[[paste0(method_pair[1], "_index_b")]]
  method2_index <- FeatureData.b.two[[paste0(method_pair[2], "_index_b")]]
  jaccard.index.sim(method1_index, method2_index)
}, simplify = TRUE)

# Assign names to the results
names(ji_omic_two_pairs_b) <- combn(methods, 2, function(method_pair) {
  paste(method_pair[1], method_pair[2], "ji_phi", sep = "_")
})

# Convert to a named vector or data frame for further use
ji_omic_two_pairs_b <- as.data.frame(as.list(ji_omic_two_pairs_b))

# SAMPLES
# Calculate pairwise Jaccard indices dynamically
ji_smp_pairs_b <- combn(methods, 2, function(method_pair) {
  method1_index <- SampleData3[[paste0(method_pair[1], "_index_b")]]
  method2_index <- SampleData3[[paste0(method_pair[2], "_index_b")]]
  jaccard.index.sim(method1_index, method2_index)
}, simplify = TRUE)

# Assign names to the results
names(ji_smp_pairs_b) <- combn(methods, 2, function(method_pair) {
  paste(method_pair[1], method_pair[2], "ji_phi", sep = "_")
})

# Convert to a named vector or data frame for further use
ji_smp_pairs_b <- as.data.frame(as.list(ji_smp_pairs_b))

# Organize the results for the current method
jaccard_comparison <- list(#ji_true_pairs_a=ji_true_pairs_a, 
  ji_omic_one_pairs_a=ji_omic_one_pairs_a, 
  ji_omic_two_pairs_a=ji_omic_two_pairs_a, ji_smp_pairs_a=ji_smp_pairs_a,
  #ji_true_pairs_b=ji_true_pairs_b, 
  ji_omic_one_pairs_b=ji_omic_one_pairs_b, 
  ji_omic_two_pairs_b=ji_omic_two_pairs_b, ji_smp_pairs_b=ji_smp_pairs_b)

# Create a unique name for each dataset
dataset_name <- sprintf("variance_%d_iteration_%d_varphi_%d", sigma, i, varphi) 
jaccard_comparison_results[[dataset_name]] <- jaccard_comparison


# COMPARISON METRICS: SENSITIVITY, SPECIFICITY; ACCURACY
methods_a <- c("fabia_index_a", "mofa_index_a", "mfa_index_a", "gfa_index_a") # factor 1
methods_b <- c("fabia_index_b", "mofa_index_b", "mfa_index_b", "gfa_index_b")

metric_samples_a <- CompareMetrics(data = SampleData3, signal_index = "signal_index_a", methods = methods_a)
metric_samples_b <- CompareMetrics(data = SampleData3, signal_index = "signal_index_b", methods = methods_b)
metric_omic.one_a <- CompareMetrics(data = FeatureData.a.one, signal_index = "signal_index_a", methods = methods_a)
metric_omic.one_b <- CompareMetrics(data = FeatureData.a.one, signal_index = "signal_index_b", methods = methods_b)
metric_omic.two_a <- CompareMetrics(data = FeatureData.b.two, signal_index = "signal_index_a", methods = methods_a)
metric_omic.two_b <- CompareMetrics(data = FeatureData.b.two, signal_index = "signal_index_b", methods = methods_b)

# Organize the results for the current method
per_measures <- list(metric_samples_a = metric_samples_a, 
                     metric_samples_b = metric_samples_b, 
                     metric_omic.one_a = metric_omic.one_a,
                     metric_omic.one_b = metric_omic.one_b, 
                     metric_omic.two_a = metric_omic.two_a,
                     metric_omic.two_b = metric_omic.two_b)

# Create a unique name for each dataset
dataset_name <- sprintf("variance_%d_iteration_%d_varphi_%d", sigma, i, varphi) 
per_measures_results[[dataset_name]] <- per_measures

method_output <- list(
  simulated_data = simulated_data,
  fabia_result = fabia_result,
  mofa_result = mofa_result,
  mfa_result = mfa_result,
  gfa_result = gfa_result
) 

dataset <- sprintf("variance_%d_iteration_%d", sigma, i)
dataset_output[[dataset]] <- method_output
#} iter

method_output_corrected <- list(
  merged_scores_F1 = merged_scores_F1,
  merged_scores_F2 = merged_scores_F2,
  merged_loadings_omic1_F1 = merged_loadings_omic1_F1,
  merged_loadings_omic1_F2 = merged_loadings_omic1_F2,
  merged_loadings_omic2_F1 = merged_loadings_omic2_F1,
  merged_loadings_omic2_F2 = merged_loadings_omic2_F2
) 

dataset <- sprintf("variance_%d_iteration_%d", sigma, i)
dataset_output_corrected[[dataset]] <- method_output_corrected

}

 }


}
# ---------------- simulation results ------------------------
  # Combine results into a list
  results <- list(
    jaccard_results = jaccard_results,
    jaccard_comparison_results = jaccard_comparison_results,
    per_measures_results = per_measures_results,#,
    dataset_output = dataset_output,
    dataset_output_corrected = dataset_output_corrected
  )
return(results)

}

#  -------------------------- Simulate multi-factor: Define the actual simulation parameters -------------------------- #

gfa_seed <- sample(1e6, 1)
set.seed(gfa_seed)

sigmas <- 20  # Define sigma or the range of sigmas

# Loop over sigmas
for (sigma in sigmas) {
  # Loop over methods for the current sigma
  #for (method in methods) {
    # Generate a unique result name
    result_name <- paste0("sim_MomicResults_multifactor_sigma_", sigma)
    
    # Run the simulation for the current method and sigma
    simulation_result <- simulate_MOmicsR(
      n_features_one = 4000, 
      n_features_two = 3000, 
      n_samples = 100, 
      var_sigma = sigma, 
      num_biclusters = 2,
      num_iteration = 1) 
    
    # Save the result to a file
    saveRDS(
      sim_result,
      #file = paste0("/user/leuven/364/vsc36498/", result_name, ".rds")
      #file = paste0("C:/Users/Lenovo/Downloads/", result_name, ".rds")
      file = paste0("C:/Users/bosangir/Downloads/", result_name, ".rds")
    )                                  
  #}
  
  # Print progress message
  cat("Completed simulations for sigma =", sigma, "\n")
} 
gfa_seed
# Print completion message
cat("All simulations completed and saved!\n")

