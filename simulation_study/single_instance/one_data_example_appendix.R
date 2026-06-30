set.seed(645) #644

# Load the required libraries
library(ggplot2)
library(tidyr)
library(viridis)
library(fabia)
library(FactoMineR)
library(readr)
library(readxl)
library(tidyverse)
library(dplyr)
library(stringr)
library(basilisk)
library(MOFA2)
library(data.table)
library(GFA)
library(pROC)
library(mclust)
library(zoo)

# -------------------------- SIMULATION: Define the functions to simulate data -----------------------------
# shared.factor function
shared.factor <- function(n_features_one, sig_omic_one = c(3650, n_features_one), n_features_two, sig_omic_two = c(n_features_two, 2), n_samples, sigmas, iterations) {
  # FABIA check function
  func_fabia <- function(data, num) {
    #set.seed(123)
    fabia_object <- fabia(as.matrix(data), p = 1, alpha = 0.01, cyc = 1000, spl = 0.5, spz = 0.5, 
                          random = 1.0, center = 2, norm = 2, lap = 1.0, nL = 1)
    fab_loading <- fab_loading(fabia_object, num)
    fab_score <- fab_score(fabia_object, num)
    X <- fabia_object@X
    fabia_result <- list(fab_loading = fab_loading, fab_score = fab_score, X = X)
    return(fabia_result)
  }
  
  # Initialize lists to store omic data for all iterations and sigmas
  all_omic_data <- list()
  
  for (iter in 1:iterations) {  # Iterate through the specified number of iterations
    omic.one <- list()
    omic.two <- list()
    
    for (k in 1:length(sigmas)) {  # Iterate through the list of sigmas
      n_s <- n_samples
      valid_data <- FALSE  # Flag for FABIA validation
      
      while (!valid_data) {
        # Generate first OMIC data
        s_sig_s <- ceiling(n_s / 5.3)
        s_sig_e <- ceiling(n_s / 2.8)
        alpha <- rnorm(n_s, 0, 1.0)
        assigned_indices_samples <- sample(s_sig_s:s_sig_e, length(s_sig_s:s_sig_e))
        alpha[assigned_indices_samples] <- rnorm(length(s_sig_s:s_sig_e), 4.5, 0.5)
        d_sig_s <- sig_omic_one[1]#ceiling(n_features_one / 1.11)
        d_sig_e <- ceiling(n_features_one / 1)
        beta <- rnorm(n_features_one, 0, 1.0)
        assigned_indices_features <- sample(d_sig_s:d_sig_e, length(d_sig_s:d_sig_e))
        beta[assigned_indices_features] <- rnorm(length(d_sig_s:d_sig_e), 5, 0.5)
        data.1 <- alpha %*% t(beta)
        eps <- rnorm(n_s * n_features_one, 0, sigmas[k])
        omic1_data <- matrix(data.1, n_s, n_features_one) + matrix(eps, n_s, n_features_one)
        colnames(omic1_data) <- paste0('omic1_feature_', seq_len(n_features_one))
        rownames(omic1_data) <- paste0('sample_', seq_len(n_s))
        
        # Generate second OMIC data
        #gamma <- rnorm(n_s, 0, 1.5)
        assigned_indices_samples2 <- sample(s_sig_s:s_sig_e, length(s_sig_s:s_sig_e))
        gamma = alpha#[assigned_indices_samples2] <- rnorm(length(s_sig_s:s_sig_e), 4.5, 0.25)
        d_sig2_s <- 1#ceiling(n_features_two / 1)
        d_sig2_e <- 335#ceiling(n_features_two + 334)# / sample(8:9, 1))
        delta <- rnorm(n_features_two, 0, 1.0)
        assigned_indices_features2 <- sample(d_sig2_s:d_sig2_e, length(d_sig2_s:d_sig2_e))
        delta[assigned_indices_features2] <- rnorm(length(d_sig2_s:d_sig2_e), 5.5, 0.5)
        data.2 <- gamma %*% t(delta)
        eps2 <- rnorm(n_s * n_features_two, 0, sigmas[k])
        omic2_data <- matrix(data.2, n_s, n_features_two) + matrix(eps2, n_s, n_features_two)
        colnames(omic2_data) <- paste0('omic2_feature_', seq_len(n_features_two))
        rownames(omic2_data) <- paste0('sample_', seq_len(n_s))
        
        # Concatenate datasets
        concatenated_data <- cbind(omic1_data, omic2_data)
        
        # Run FABIA and check conditions
        fabia_result <- func_fabia(concatenated_data, num = 1)
        fab_loading <- fabia_result$fab_loading$loading_FABIA
        fab_score <- fabia_result$fab_score$score_FABIA
        
        if (!any(is.na(fab_loading)) && mean(fab_loading, na.rm = TRUE) != 0 &&
            !any(is.na(fab_score)) && mean(fab_score, na.rm = TRUE) != 0) {
          valid_data <- TRUE
        } else {
          message("FABIA validation failed, regenerating concatenated_data...")
        }
      }
      
      omic.one[[k]] <- omic1_data
      omic.two[[k]] <- omic2_data
    }
    
    simulated_datasets <- list(object.one = omic.one, object.two = omic.two)
    concatenated_datasets <- list()
    for (i in seq_along(simulated_datasets$object.one)) {
      concatenated_datasets[[i]] <- cbind(simulated_datasets$object.one[[i]], simulated_datasets$object.two[[i]])
    }
    
    all_omic_data[[paste0("iteration_", iter)]] <- list(
      concatenated_datasets = concatenated_datasets,
      indices_features.1 = assigned_indices_features,
      indices_samples.1 = assigned_indices_samples,
      indices_features.2 = assigned_indices_features2,
      indices_samples.2 = assigned_indices_samples2,
      sample_sig_start = s_sig_s,
      sample_sig_end = s_sig_e,
      feature_sig_start = d_sig_s,
      feature_sig_end = d_sig_e,
      feature_sig2_start = d_sig2_s,
      feature_sig2_end = d_sig2_e,
      list_alphas = alpha, list_gammas = gamma, list_betas = beta, list_deltas = delta
    )
  }
  
  return(all_omic_data)
}

# -------------- Actual Simulation ------------------

sigma  = 10
iter = 1

set.seed(123)
simulated_data <- shared.factor(
  n_features_one = 4000,
  n_features_two = 3000,
  n_samples = 100,
  sigmas = sigma,
  iterations = iter 
)

# -------------- Raw Data Visualization -----------------
plot_simData(sim_object = simulated_data, data = "merged", type = "heatmap")

# scores
alpha_score = simulated_data$iteration_1$list_alphas

# data
sample_id <- 1:100

# Define signal range
sample_sig_start = simulated_data$iteration_1$sample_sig_start 
sample_sig_end = simulated_data$iteration_1$sample_sig_end 

# Assign group labels
group <- ifelse(sample_id >= sample_sig_start & sample_id <= sample_sig_end, "signal", "noise")

# Combine into a dataframe
df_alpha_factor <- data.frame(sample_id, alpha_score, group)

# Plot without horizontal line
ggplot(df_alpha_factor, aes(x = sample_id, y = alpha_score, color = group)) +
  geom_point(size = 2.1, alpha = 0.95) +
  scale_color_manual(values = c("signal" = "#E69F00", "noise" = "#56B4E9")) +
  labs(title = "", # Alpha Factor Scores by Sample ID (Scatter)
       x = "Sample ID",
       y = "Factor Score",
       color = "Group") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )

#  Weights

feature_id1 = 1:4000
omic1_weight = simulated_data[["iteration_1"]][["list_betas"]]
# Assign group labels
indices_features.1 = simulated_data[["iteration_1"]][["indices_features.1"]] 
group1 <- ifelse(feature_id1 %in% indices_features.1, "signal", "noise")
table(group1)
# Combine into a dataframe
df_weight_omic1 <- data.frame(feature_id1, omic1_weight, group1)

# Plot without horizontal line
ggplot(df_weight_omic1, aes(x = feature_id1, y = omic1_weight, color = group1)) +
  geom_point(size = 2.1, alpha = 0.95) +
  scale_color_manual(values = c("signal" = "#E69F00", "noise" = "#56B4E9")) +
  labs(title = "", # Alpha Factor Scores by Sample ID (Scatter)
       x = "Feature ID",
       y = "Loadings",
       color = "Group") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )

feature_id2 = 1:3000
omic2_weight = simulated_data[["iteration_1"]][["list_deltas"]]
# Assign group labels
indices_features.2 = simulated_data[["iteration_1"]][["indices_features.2"]] 
group2 <- ifelse(feature_id2 %in% indices_features.2, "signal", "noise")
table(group2)
# Combine into a dataframe
df_weight_omic2 <- data.frame(feature_id2, omic2_weight, group2)

# Plot without horizontal line
ggplot(df_weight_omic2, aes(x = feature_id2, y = omic2_weight, color = group2)) +
  geom_point(size = 2.1, alpha = 0.95) +
  scale_color_manual(values = c("signal" = "#E69F00", "noise" = "#56B4E9")) +
  labs(title = "", # Alpha Factor Scores by Sample ID (Scatter)
       x = "Feature ID",
       y = "Loadings",
       color = "Group") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )

# ---------- Multi-omics Analysis --------------------
# MOFA, FABIA, MFA, GFA
## -------------- MOFA ----------------------
library(MOFAdata) 
library(MOFA2); 
library(data.table); library(gridExtra); library(ggplot2); library(tidyverse); library(dplyr); 
library(GFA); library(fabia)

# Access variables from the shared environment
n_features_one <- 4000
n_features_two <- 3000

first_omic <- simulated_data$iteration_1$concatenated_datasets[[1]][, 1:n_features_one]
second_omic <- simulated_data$iteration_1$concatenated_datasets[[1]][, (n_features_one + 1):(n_features_one + n_features_two)]

sim_data <- list(`Omic 1` = as.matrix(t(first_omic)), `Omic 2` = as.matrix(t(second_omic)))

# Create the MOFA object and train the model
MOFAobject <- create_mofa(sim_data)
MOFAobject

# Plot data overview
plot_data_overview(MOFAobject)

# Define MOFA options
## Data options
data_opts <- get_default_data_options(MOFAobject)
data_opts

## Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 2
model_opts

## Training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
train_opts

## Train the MOFA model
MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts
)
outfile_object_sim = paste0(getwd(),"model_object_sim.hdf5")
MOFAobject.sim.trained <- run_mofa(MOFAobject, outfile_object_sim)

# Overview of the trained MOFA model
## Slots - The MOFA object consists of multiple slots where relevant data and information is stored. For descriptions, you can read the documentation using ?MOFA. The most important slots are:
# data: input data used to train the model (features are centered at zero mean)
# samples_metadata: sample metadata information
# expectations: expectations of the posterior distributions for the Weights and the Factors
slotNames(MOFAobject.sim.trained)

## Add sample metadata to the model
samples_metadata = samples_metadata(MOFAobject.sim.trained)

## Correlation between factors # A good sanity check is to verify that the Factors are largely uncorrelated. 
plot_factor_cor(MOFAobject.sim.trained)

# Explore Factors
(MOFAobject.sim.trained@expectations[["Z"]][["group1"]])

## Plot variance decomposition
### Variance decomposition by Factor
plot_variance_explained(MOFAobject.sim.trained, max_r2=15)

### Total variance explained per view
plot_variance_explained(MOFAobject.sim.trained, plot_total = T)[[2]]

# Factor scores
factor_mofa <- get_factors(MOFAobject.sim.trained, factors = 'all', as.data.frame = T)
factor_mofa <- factor_mofa %>%
  spread(factor, value)
factor_mofa2 = factor_mofa[, !(colnames(factor_mofa) %in% c('group'))]
colnames(factor_mofa2) =  c('sample', 'score_MOFA1', 'score_MOFA2')
factor_mofa2 = factor_mofa2[, -3]

# visualize
factor_mofa2$sample_id <- as.numeric(gsub("sample_", "", factor_mofa2$sample))
# Define signal range
sample_sig_start = simulated_data$iteration_1$sample_sig_start 
sample_sig_end = simulated_data$iteration_1$sample_sig_end 

# Assign group labels
group <- ifelse(sample_id >= sample_sig_start & sample_id <= sample_sig_end, "signal", "noise")
sample_id = 1:100
group_data <- data.frame(sample_id, group)

# Combine into a dataframe
df_mofa_factor <- merge(factor_mofa2, group_data, by = "sample_id")

# Plot without horizontal line
ggplot(df_mofa_factor, aes(x = sample_id, y = score_MOFA1, color = group)) +
  geom_point(size = 2.1, alpha = 0.95) +
  scale_color_manual(values = c("signal" = "#E69F00", "noise" = "#56B4E9")) +
  labs(title = "", # Alpha Factor Scores by Sample ID (Scatter)
       x = "Sample ID",
       y = "Factor Score",
       color = "Group") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )

#  Weights


# Weights/Loading 
loading_mofa<- get_weights(MOFAobject.sim.trained, factors = 'all', as.data.frame = T)
loading_mofa2 <- loading_mofa %>%
  spread(factor, value)
colnames(loading_mofa2) =  c('feature', 'view', 'loading_MOFA1', 'loading_MOFA2')

omic1 <- loading_mofa2 %>%
  filter(view == 'Omic 1') %>%
  select(feature, view, loading_MOFA1)
omic1$feature_id1 <- as.numeric(gsub("omic1_feature_", "", omic1$feature))

feature_id1 = 1:4000
# Assign group labels
indices_features.1 = simulated_data[["iteration_1"]][["indices_features.1"]] 
group1 <- ifelse(feature_id1 %in% indices_features.1, "signal", "noise")

group_omic1 = data.frame(feature_id1, group1)

# Combine into a dataframe
df_weight_omic1 <- merge(omic1, group_omic1, by = 'feature_id1')

# Plot without horizontal line
ggplot(df_weight_omic1, aes(x = feature_id1, y = loading_MOFA1, color = group1)) +
  geom_point(size = 2.1, alpha = 0.95) +
  scale_color_manual(values = c("signal" = "#E69F00", "noise" = "#56B4E9")) +
  labs(title = "", # Alpha Factor Scores by Sample ID (Scatter)
       x = "Feature ID",
       y = "Loadings",
       color = "Group") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )


omic2 <- loading_mofa2 %>%
  filter(view == 'Omic 2') %>%
  select(feature, view, loading_MOFA1)
omic2$feature_id2 <- as.numeric(gsub("omic2_feature_", "", omic2$feature))

feature_id2 = 1:3000
# Assign group labels
indices_features.2 = simulated_data[["iteration_1"]][["indices_features.2"]] 
group2 <- ifelse(feature_id2 %in% indices_features.2, "signal", "noise")

group_omic2 = data.frame(feature_id2, group2)

# Combine into a dataframe
df_weight_omic2 <- merge(omic2, group_omic2, by = 'feature_id2')

# Plot without horizontal line
ggplot(df_weight_omic2, aes(x = feature_id2, y = loading_MOFA1, color = group2)) +
  geom_point(size = 2.1, alpha = 0.95) +
  scale_color_manual(values = c("signal" = "#E69F00", "noise" = "#56B4E9")) +
  labs(title = "", # Alpha Factor Scores by Sample ID (Scatter)
       x = "Feature ID",
       y = "Loadings",
       color = "Group") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )


## -------------- FABIA ----------------------

# merge the two datasets
sim_data_fabia <- rbind(t(first_omic), t(second_omic))

# FABIA MODEL
set.seed(123)
FABIAobject.sim.trained <- fabia(sim_data_fabia,
                                 p = 2,           # number of hidden factors = number of biclusters, default = 5
                                 alpha = 0.01,     # sparseness loadings (0 - 1.0); default = 0.1
                                 cyc = 1000,      # number of iterations; default = 500
                                 spl = 0.5,       # sparseness prior loadings (0 - 2.0); default = 0 (Laplace)
                                 spz = 0.5,       # sparseness factors (0.5 - 2.0); default = 0.5 (Laplace)
                                 random = 1.0,    # random initialization of loadings in [-random,random]; default = 1.0.
                                 center = 2,      # data centering: 1 (mean), 2 (median), > 2 (mode), 0 (no); default = 2
                                 norm = 2,        # data normalization: 1 (0.75-0.25 quantile), >1 (var=1), 0 (no); default = 1
                                 lap = 1.0,       # minimal value of the variational parameter, default = 1
                                 nL = 1           # maximal number of biclusters at which a row element can participate; default = 0 (no limit)
)

fabia.sim.scaledData <- FABIAobject.sim.trained@X; 
fabia.sim.scaledData <- as.matrix(fabia.sim.scaledData)

# samples
samples_info.sim = data.frame(rownames(t(fabia.sim.scaledData))); colnames(samples_info.sim) <- c("samples")

# fabia scores
data_factors_fabia.sim = data.frame(t(FABIAobject.sim.trained@Z)) 
data_factors_fabia.sim$samples = rownames(data_factors_fabia.sim)
colnames(data_factors_fabia.sim) <- c("score_FABIA1", "score_FABIA2", "sample")
df_fabia_factor <- data_factors_fabia.sim
df_fabia_factor$sample_id <- as.numeric(gsub("sample_", "", df_fabia_factor$sample))
df_fabia_factor = merge(df_fabia_factor, group_data, by = "sample_id")

# Plot without horizontal line
ggplot(df_fabia_factor, aes(x = sample_id, y = (-1)*(score_FABIA1), color = group)) +
  geom_point(size = 2.1, alpha = 0.95) +
  scale_color_manual(values = c("signal" = "#E69F00", "noise" = "#56B4E9")) +
  labs(title = "", # Alpha Factor Scores by Sample ID (Scatter)
       x = "Sample ID",
       y = "Factor Score",
       color = "Group") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )

# Weights/Loading


# features
features_info.sim = data.frame(rownames(fabia.sim.scaledData)) 
colnames(features_info.sim) <- c("features") 
features_info.sim$ID <- as.integer(row.names(features_info.sim))
features_info.sim$view <- NA  # Create a new column with NAs
features_info.sim$view[1:4000] <- "Omic 1"  # Assign 'x' to rows 1 to 50
features_info.sim$view[4001:7000] <- 'Omic 2'  # Assign 'y' to rows 51 to 100


loading_fabia = data.frame(FABIAobject.sim.trained@L)
loading_fabia$feature = rownames(loading_fabia)
colnames(loading_fabia) <- c("loading_FABIA1", "loading_FABIA2","features")

loading_fabia = merge(loading_fabia, features_info.sim, by = "features")

omic1_fabia <- loading_fabia %>%
  filter(view == 'Omic 1') %>%
  select(features, view, loading_FABIA1, ID)
omic1_fabia$feature_id1 = omic1_fabia$ID 

# Assign group labels
indices_features.1 = simulated_data[["iteration_1"]][["indices_features.1"]] 
group1 <- ifelse(feature_id1 %in% indices_features.1, "signal", "noise")

group_omic1 = data.frame(feature_id1, group1)

# Combine into a dataframe
df_weight_omic1_fabia <- merge(omic1_fabia, group_omic1, by = 'feature_id1')

# Plot without horizontal line
ggplot(df_weight_omic1_fabia, aes(x = feature_id1, y = (-1)*(loading_FABIA1), color = group1)) +
  geom_point(size = 2.1, alpha = 0.95) +
  scale_color_manual(values = c("signal" = "#E69F00", "noise" = "#56B4E9")) +
  labs(title = "", # Alpha Factor Scores by Sample ID (Scatter)
       x = "Feature ID",
       y = "Loadings",
       color = "Group") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )


omic2_fabia <- loading_fabia %>%
  filter(view == 'Omic 2') %>%
  select(features, view, loading_FABIA1, ID)
omic2_fabia$feature_id2 <- as.numeric(gsub("omic2_feature_", "", omic2_fabia$features))

# Assign group labels
indices_features.2 = simulated_data[["iteration_1"]][["indices_features.2"]] 
group2 <- ifelse(feature_id2 %in% indices_features.2, "signal", "noise")

group_omic2 = data.frame(feature_id2, group2)

# Combine into a dataframe
df_weight_omic2_fabia <- merge(omic2_fabia, group_omic2, by = 'feature_id2')

# Plot without horizontal line
ggplot(df_weight_omic2_fabia, aes(x = feature_id2, y = (-1)*(loading_FABIA1), color = group2)) +
  geom_point(size = 2.1, alpha = 0.95) +
  scale_color_manual(values = c("signal" = "#E69F00", "noise" = "#56B4E9")) +
  labs(title = "", # Alpha Factor Scores by Sample ID (Scatter)
       x = "Feature ID",
       y = "Loadings",
       color = "Group") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )

## -------------- MFA ----------------------
set.seed(123)
simdata_mfa <- rbind(t(first_omic), t(second_omic)) 
simdata_mfa <- as.data.frame(t(simdata_mfa))

MFAobject.sim.trained <- FactoMineR::MFA(simdata_mfa, group = c(4000, 3000), type = c("s", "s"),
                                         name.group = c("Omic 1", "Omic 2"),
                                         graph = FALSE)

# Factor scores
factor_mfa = data.frame(MFAobject.sim.trained$ind$coord)
colnames(factor_mfa) <- c("score_MFA1", "score_MFA2", "score_MFA3", "score_MFA4", "score_MFA5")
factor_mfa$sample = rownames(factor_mfa)
factor_mfa2 = factor_mfa[,-c(2:5)]

# visualize
factor_mfa2$sample_id <- as.numeric(gsub("sample_", "", factor_mfa2$sample))
# Define signal range
sample_sig_start = simulated_data$iteration_1$sample_sig_start 
sample_sig_end = simulated_data$iteration_1$sample_sig_end 

# Assign group labels
group <- ifelse(sample_id >= sample_sig_start & sample_id <= sample_sig_end, "signal", "noise")
sample_id = 1:100
group_data <- data.frame(sample_id, group)

# Combine into a dataframe
df_mfa_factor <- merge(factor_mfa2, group_data, by = "sample_id")

# Plot without horizontal line
ggplot(df_mfa_factor, aes(x = sample_id, y = score_MFA1, color = group)) +
  geom_point(size = 2.1, alpha = 0.95) +
  scale_color_manual(values = c("signal" = "#E69F00", "noise" = "#56B4E9")) +
  labs(title = "", # Alpha Factor Scores by Sample ID (Scatter)
       x = "Sample ID",
       y = "Factor Score",
       color = "Group") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )

# Weights/Loading 
loading_mfa = data.frame(MFAobject.sim.trained$quanti.var$coord)
colnames(loading_mfa) <- c("loading_MFA1", "loading_MFA2", "loading_MFA3", "loading_MFA4", "loading_MFA5")
loading_mfa$features = rownames(loading_mfa)
loading_mfa2 = loading_mfa[, -c(2:5)]

# Weights/Loading 

loading_mfa2 <- merge(loading_mfa2, features_info.sim, by = "features")

omic1_mfa <- loading_mfa2 %>%
  filter(view == 'Omic 1') %>%
  select(feature, view, loading_MFA1)
omic1_mfa$feature_id1 <- as.numeric(gsub("omic1_feature_", "", omic1_mfa$feature))

feature_id1 = 1:4000
# Assign group labels
indices_features.1 = simulated_data[["iteration_1"]][["indices_features.1"]] 
group1 <- ifelse(feature_id1 %in% indices_features.1, "signal", "noise")

group_omic1 = data.frame(feature_id1, group1)

# Combine into a dataframe
df_weight_omic1_mfa <- merge(omic1_mfa, group_omic1, by = 'feature_id1')

# Plot without horizontal line
ggplot(df_weight_omic1_mfa, aes(x = feature_id1, y = loading_MFA1, color = group1)) +
  geom_point(size = 2.1, alpha = 0.95) +
  scale_color_manual(values = c("signal" = "#E69F00", "noise" = "#56B4E9")) +
  labs(title = "", # Alpha Factor Scores by Sample ID (Scatter)
       x = "Feature ID",
       y = "Loadings",
       color = "Group") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )


omic2_mfa <- loading_mfa2 %>%
  filter(view == 'Omic 2') %>%
  select(feature, view, loading_MFA1)
omic2_mfa$feature_id2 <- as.numeric(gsub("omic2_feature_", "", omic2_mfa$feature))

feature_id2 = 1:3000
# Assign group labels
indices_features.2 = simulated_data[["iteration_1"]][["indices_features.2"]] 
group2 <- ifelse(feature_id2 %in% indices_features.2, "signal", "noise")

group_omic2 = data.frame(feature_id2, group2)

# Combine into a dataframe
df_weight_omic2_mfa <- merge(omic2_mfa, group_omic2, by = 'feature_id2')

# Plot without horizontal line
ggplot(df_weight_omic2_mfa, aes(x = feature_id2, y = loading_MFA1, color = group2)) +
  geom_point(size = 2.1, alpha = 0.95) +
  scale_color_manual(values = c("signal" = "#E69F00", "noise" = "#56B4E9")) +
  labs(title = "", # Alpha Factor Scores by Sample ID (Scatter)
       x = "Feature ID",
       y = "Loadings",
       color = "Group") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )

## -------------------- GFA-----------------
fabia.scaledData.sim <- FABIAobject.sim.trained@X
gfa_sim.data <- as.data.frame(fabia.scaledData.sim)
merged_GFA_sim.data = list(t(gfa_sim.data))
model_option <- getDefaultOpts()
model_option$iter.max <- 1000
model_option$iter.burnin <- 10
GFAobject.sim.trained <- gfa(t(merged_GFA_sim.data), K= 1, opts=model_option)

# Factor scores
factor_gfa = as.data.frame(GFAobject.sim.trained$X)
colnames(factor_gfa) <- c("score_GFA1")
factor_gfa$sample = rownames(factor_gfa)

# visualize
factor_gfa$sample_id <- as.numeric(gsub("sample_", "", factor_gfa$sample))

# Define signal range
sample_sig_start = simulated_data$iteration_1$sample_sig_start 
sample_sig_end = simulated_data$iteration_1$sample_sig_end 

# Assign group labels
group <- ifelse(sample_id >= sample_sig_start & sample_id <= sample_sig_end, "signal", "noise")
sample_id = 1:100
group_data <- data.frame(sample_id, group)

# Combine into a dataframe
df_gfa_factor <- merge(factor_gfa, group_data, by = "sample_id")

# Plot without horizontal line
ggplot(df_gfa_factor, aes(x = sample_id, y = score_GFA1, color = group)) +
  geom_point(size = 2.1, alpha = 0.95) +
  scale_color_manual(values = c("signal" = "#E69F00", "noise" = "#56B4E9")) +
  labs(title = "", # Alpha Factor Scores by Sample ID (Scatter)
       x = "Sample ID",
       y = "Factor Score",
       color = "Group") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )

# Weights/Loading 
loading_gfa = as.data.frame(GFAobject.sim.trained$W)
colnames(loading_gfa) <- c("loading_GFA1")
loading_gfa$features = rownames(loading_gfa)

loading_gfa <- merge(loading_gfa, features_info.sim, by = "features")

omic1_gfa <- loading_gfa %>%
  filter(view == 'Omic 1') %>%
  select(feature, view, loading_GFA1)
omic1_gfa$feature_id1 <- as.numeric(gsub("omic1_feature_", "", omic1_gfa$feature))

feature_id1 = 1:4000
# Assign group labels
indices_features.1 = simulated_data[["iteration_1"]][["indices_features.1"]] 
group1 <- ifelse(feature_id1 %in% indices_features.1, "signal", "noise")

group_omic1 = data.frame(feature_id1, group1)

# Combine into a dataframe
df_weight_omic1_gfa <- merge(omic1_gfa, group_omic1, by = 'feature_id1')

# Plot without horizontal line
ggplot(df_weight_omic1_gfa, aes(x = feature_id1, y = loading_GFA1, color = group1)) +
  geom_point(size = 2.1, alpha = 0.95) +
  scale_color_manual(values = c("signal" = "#E69F00", "noise" = "#56B4E9")) +
  labs(title = "", # Alpha Factor Scores by Sample ID (Scatter)
       x = "Feature ID",
       y = "Loadings",
       color = "Group") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )


omic2_gfa <- loading_gfa %>%
  filter(view == 'Omic 2') %>%
  select(feature, view, loading_GFA1)
omic2_gfa$feature_id2 <- as.numeric(gsub("omic2_feature_", "", omic2_gfa$feature))

feature_id2 = 1:3000

# Assign group labels
indices_features.2 = simulated_data[["iteration_1"]][["indices_features.2"]] 
group2 <- ifelse(feature_id2 %in% indices_features.2, "signal", "noise")

group_omic2 = data.frame(feature_id2, group2)

# Combine into a dataframe
df_weight_omic2_gfa <- merge(omic2_gfa, group_omic2, by = 'feature_id2')

# Plot without horizontal line
ggplot(df_weight_omic2_gfa, aes(x = feature_id2, y = loading_GFA1, color = group2)) +
  geom_point(size = 2.1, alpha = 0.95) +
  scale_color_manual(values = c("signal" = "#E69F00", "noise" = "#56B4E9")) +
  labs(title = "", # Alpha Factor Scores by Sample ID (Scatter)
       x = "Feature ID",
       y = "Loadings",
       color = "Group") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )

# -------------- Factor & Loading
# Ensure all data frames have a common column for joining, e.g., "sample"
mofa_factor <- df_mofa_factor
fabia_factor <- df_fabia_factor[c("sample_id", "score_FABIA1")]
mfa_factor <- df_mfa_factor[c("sample_id", "score_MFA1")]
gfa_factor <- df_gfa_factor[c("sample_id", "score_GFA1")]

# ----------- Main Analysis ---------------
## ----------- Factor Analysis ---------------
# Inner join on "sample" column
factors_data <- Reduce(function(x, y) merge(x, y, by = "sample_id", all = FALSE), 
                         list(gfa_factor, mfa_factor, mofa_factor, fabia_factor))

# Visuaization of scores using correlation matrix plot
library(GGally); library(psych)

# Define new variable names
new_var_names <- c("Score (FABIA)", "Score (MOFA)", "Score (MFA)", "Score (GFA)")
colnames(data_scores_cor)[c(7,5,3,2)] <- new_var_names

# Plot with improved styling
pairs.panels(
  data_scores_cor[new_var_names],
  method = "pearson",          # Use Pearson correlation
  hist.col = "lightblue",      # Softer color for histograms
  breaks = 10,
  density = TRUE,              # Add density lines
  ellipses = TRUE,             # Add confidence ellipses
  lm = TRUE,                   # Add regression lines
  #smoother = TRUE,             # Add loess smoother
  stars = TRUE,                # Significance stars
  cex = 0.6,                   # Point size
  pch = 21,                    # Better point shape
  bg = "black",              # Background fill for points
  main = "",
  font.main = 2,               # Bold title
  col = "black"                # Point border color
)
# ------------------------------------- VISUALIZATION OF SCORES -------------------------------------------

# scaling
scale_matrix <- function(x) {
  min_x <- min(x, na.rm = TRUE)
  max_x <- max(x, na.rm = TRUE)
  return(2 * (x - min_x) / (max_x - min_x) - 1)
}

scale_data_scores_cor <- data_scores_cor[c(7,5,3,2)]
scale_data_scores_cor$`Score (FABIA)` = scale_matrix(scale_data_scores_cor$`Score (FABIA)`)
scale_data_scores_cor$`Score (MOFA)` = scale_matrix(scale_data_scores_cor$`Score (MOFA)`)
scale_data_scores_cor$`Score (MFA)` = scale_matrix(scale_data_scores_cor$`Score (MFA)`)
scale_data_scores_cor$`Score (GFA)` = scale_matrix(scale_data_scores_cor$`Score (GFA)`)

scale_data_scores_cor1 <- na.omit(scale_data_scores_cor)


ggplot(scale_data_scores_cor1, aes(x = sample_id, y = `Score (GFA)`, color = `Score (GFA)`)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # <- this adds the line
  scale_color_gradientn(colors = viridis::viridis(100)) +
  labs(x = "Sample ID", y = "Factor Score", color = NULL) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )

ggplot(scale_data_scores_cor1, aes(x = sample, y = `Score (MOFA)`, color = `Score (MOFA)`)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # <- this adds the line
  scale_color_gradientn(colors = viridis::viridis(100)) +
  labs(x = "Samples", y = "Factor1 score", color = NULL) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),       # remove x-axis text
    axis.ticks.x = element_blank(),      # remove x-axis ticks
    axis.text.y = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )
data_scores1 = na.omit(data_scores)
ggplot(data_scores1, aes(x = sample, y = score_MFA1, color = score_MFA1)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # <- this adds the line
  scale_color_gradientn(colors = viridis::viridis(100)) +
  labs(x = "Samples", y = "Factor1 score", color = NULL) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),       # remove x-axis text
    axis.ticks.x = element_blank(),      # remove x-axis ticks
    axis.text.y = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )
library(ggrepel)
library(ggplot2)

varphi.data <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    varphi <- quantile(numeric_values, probs = 0.60)[[1]]
    return(varphi)
  } else {
    return(NA)  # Return NA if not all values are numeric
  }
}
# First plot
plot1 = ggplot(data_scores_cor, aes(x = `Score (MOFA)`, y = `Score (FABIA)`)) + 
  geom_point(color = "dodgerblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -varphi.data(abs(data_scores_cor$`Score (FABIA)`)) , linetype = "dashed", color = "goldenrod4") +
  geom_hline(yintercept = varphi.data(abs(data_scores_cor$`Score (FABIA)`)) , linetype = "dashed", color = "goldenrod4") +
  geom_vline(xintercept = -varphi.data(abs(data_scores_cor$`Score (MOFA)`)) , linetype = "dashed", color = "goldenrod4") +
  geom_vline(xintercept = varphi.data(abs(data_scores_cor$`Score (MOFA)`)), linetype = "dashed", color = "goldenrod4") +
  labs(x = "Score (MOFA)", y = "Score (FABIA)") +
  geom_rect(aes(xmin = -varphi.data(abs(data_scores_cor$`Score (MOFA)`)) , xmax = -Inf, ymin = -varphi.data(abs(data_scores_cor$`Score (FABIA)`)) , ymax = -Inf), fill = "goldenrod4", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_scores_cor$`Score (MOFA)`)) , xmax = Inf, ymin =varphi.data(abs(data_scores_cor$`Score (FABIA)`)) , ymax = Inf), fill = "goldenrod4", alpha = 0.021)+
  geom_text_repel(data=data_scores_cor,
                  aes(label=sample),show.legend=FALSE,color='black')+
  theme_minimal()

# Second plot
plot2 = ggplot(data_scores_cor, aes(x = `Score (MOFA)`, y = `Score (MFA)`)) + 
  geom_point(color = "dodgerblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -varphi.data(abs(data_scores_cor$`Score (MFA)`)) , linetype = "dashed", color = "goldenrod4") +
  geom_hline(yintercept = varphi.data(abs(data_scores_cor$`Score (MFA)`)) , linetype = "dashed", color = "goldenrod4") +
  geom_vline(xintercept = -varphi.data(abs(data_scores_cor$`Score (MOFA)`)) , linetype = "dashed", color = "goldenrod4") +
  geom_vline(xintercept = varphi.data(abs(data_scores_cor$`Score (MOFA)`)), linetype = "dashed", color = "goldenrod4") +
  labs(x = "Score (MOFA)", y = "Score (MFA)") +
  geom_rect(aes(xmin = -varphi.data(abs(data_scores_cor$`Score (MOFA)`)) , xmax = -Inf, ymin = -varphi.data(abs(data_scores_cor$`Score (MFA)`)) , ymax = -Inf), fill = "goldenrod4", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_scores_cor$`Score (MOFA)`)) , xmax = Inf, ymin =varphi.data(abs(data_scores_cor$`Score (MFA)`)) , ymax = Inf), fill = "goldenrod4", alpha = 0.021)+
  geom_text_repel(data=data_scores_cor,
                  aes(label=sample),show.legend=FALSE,color='black')+
  theme_minimal()

# Third plot
plot3 = ggplot(data_scores_cor, aes(x = `Score (MOFA)`, y = `Score (GFA)`)) + 
  geom_point(color = "dodgerblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -varphi.data(abs(data_scores_cor$`Score (GFA)`)) , linetype = "dashed", color = "goldenrod4") +
  geom_hline(yintercept = varphi.data(abs(data_scores_cor$`Score (GFA)`)) , linetype = "dashed", color = "goldenrod4") +
  geom_vline(xintercept = -varphi.data(abs(data_scores_cor$`Score (MOFA)`)) , linetype = "dashed", color = "goldenrod4") +
  geom_vline(xintercept = varphi.data(abs(data_scores_cor$`Score (MOFA)`)), linetype = "dashed", color = "goldenrod4") +
  labs(x = "Score (MOFA)", y = "Score (GFA)") +
  geom_rect(aes(xmin = -varphi.data(abs(data_scores_cor$`Score (MOFA)`)) , xmax = -Inf, ymin = -varphi.data(abs(data_scores_cor$`Score (GFA)`)) , ymax = -Inf), fill = "goldenrod4", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_scores_cor$`Score (MOFA)`)) , xmax = Inf, ymin =varphi.data(abs(data_scores_cor$`Score (GFA)`)) , ymax = Inf), fill = "goldenrod4", alpha = 0.021)+
  geom_text_repel(data=data_scores_cor,
                  aes(label=sample),show.legend=FALSE,color='black')+
  theme_minimal()

# Four plot
plot4 = ggplot(data_scores_cor, aes(x = `Score (FABIA)`, y = `Score (MFA)`)) + 
  geom_point(color = "dodgerblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -varphi.data(abs(data_scores_cor$`Score (MFA)`)) , linetype = "dashed", color = "goldenrod4") +
  geom_hline(yintercept = varphi.data(abs(data_scores_cor$`Score (MFA)`)) , linetype = "dashed", color = "goldenrod4") +
  geom_vline(xintercept = -varphi.data(abs(data_scores_cor$`Score (FABIA)`)) , linetype = "dashed", color = "goldenrod4") +
  geom_vline(xintercept = varphi.data(abs(data_scores_cor$`Score (FABIA)`)), linetype = "dashed", color = "goldenrod4") +
  labs(x = "Score (FABIA)", y = "Score (MFA)") +
  geom_rect(aes(xmin = -varphi.data(abs(data_scores_cor$`Score (FABIA)`)) , xmax = -Inf, ymin = -varphi.data(abs(data_scores_cor$`Score (MFA)`)) , ymax = -Inf), fill = "goldenrod4", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_scores_cor$`Score (FABIA)`)) , xmax = Inf, ymin =varphi.data(abs(data_scores_cor$`Score (MFA)`)) , ymax = Inf), fill = "goldenrod4", alpha = 0.021)+
  geom_text_repel(data=data_scores_cor,
                  aes(label=sample),show.legend=FALSE,color='black')+
  theme_minimal()

# Five plot
plot5 = ggplot(data_scores_cor, aes(x = `Score (FABIA)`, y = `Score (GFA)`)) + 
  geom_point(color = "#0078BD") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#F6F6F6") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#F6F6F6") +
  geom_hline(yintercept = -varphi.data(abs(data_scores_cor$`Score (GFA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_hline(yintercept = varphi.data(abs(data_scores_cor$`Score (GFA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_vline(xintercept = -varphi.data(abs(data_scores_cor$`Score (FABIA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_vline(xintercept = varphi.data(abs(data_scores_cor$`Score (FABIA)`)), linetype = "dashed", color = "#DCDDE0") +
  labs(x = "Score (FABIA)", y = "Score (GFA)") +
  geom_rect(aes(xmin = -varphi.data(abs(data_scores_cor$`Score (FABIA)`)) , xmax = -Inf, ymin = -varphi.data(abs(data_scores_cor$`Score (GFA)`)) , ymax = -Inf), fill = "#DCDDE0", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_scores_cor$`Score (FABIA)`)) , xmax = Inf, ymin =varphi.data(abs(data_scores_cor$`Score (GFA)`)) , ymax = Inf), fill = "#DCDDE0", alpha = 0.021)+
  geom_text_repel(data=data_scores_cor,
                  aes(label=sample),show.legend=FALSE,color='black')+
  theme_minimal()

# Six plot
plot6 = ggplot(data_scores_cor, aes(x = `Score (MFA)`, y = `Score (GFA)`)) + 
  geom_point(color = "#0078BD") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#F6F6F6") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#F6F6F6") +
  geom_hline(yintercept = -varphi.data(abs(data_scores_cor$`Score (GFA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_hline(yintercept = varphi.data(abs(data_scores_cor$`Score (GFA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_vline(xintercept = -varphi.data(abs(data_scores_cor$`Score (MFA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_vline(xintercept = varphi.data(abs(data_scores_cor$`Score (MFA)`)), linetype = "dashed", color = "#DCDDE0") +
  labs(x = "Score (MFA)", y = "Score (GFA)") +
  geom_rect(aes(xmin = -varphi.data(abs(data_scores_cor$`Score (MFA)`)) , xmax = -Inf, ymin = -varphi.data(abs(data_scores_cor$`Score (GFA)`)) , ymax = -Inf), fill = "#DCDDE0", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_scores_cor$`Score (MFA)`)) , xmax = Inf, ymin =varphi.data(abs(data_scores_cor$`Score (GFA)`)) , ymax = Inf), fill = "#DCDDE0", alpha = 0.021)+
  geom_text_repel(data=data_scores_cor,
                  aes(label=sample),show.legend=FALSE,color='black')+
  theme_minimal()

# Arrange plots side by side
grid.arrange(plot1, plot2, plot3, plot4, plot5,plot6, ncol = 3, nrow = 2)


new_var_names <- c("Score (FABIA)", "Score (MOFA)", "Score (MFA)", "Score (GFA)")
colnames(data_scores1)[c(2,11,12,13)] <- new_var_names
# First plot
plot1 = ggplot(data_scores1, aes(x = `Score (MOFA)`, y = `Score (FABIA)`)) + 
  geom_point(color = "dodgerblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -varphi.data(abs(data_scores_cor$`Score (FABIA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_hline(yintercept = varphi.data(abs(data_scores_cor$`Score (FABIA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_vline(xintercept = -varphi.data(abs(data_scores_cor$`Score (MOFA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_vline(xintercept = varphi.data(abs(data_scores_cor$`Score (MOFA)`)), linetype = "dashed", color = "#DCDDE0") +
  labs(x = "Score (MOFA)", y = "Score (FABIA)") +
  geom_rect(aes(xmin = -varphi.data(abs(data_scores_cor$`Score (MOFA)`)) , xmax = -Inf, ymin = -varphi.data(abs(data_scores_cor$`Score (FABIA)`)) , ymax = -Inf), fill = "#DCDDE0", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_scores_cor$`Score (MOFA)`)) , xmax = Inf, ymin =varphi.data(abs(data_scores_cor$`Score (FABIA)`)) , ymax = Inf), fill = "#DCDDE0", alpha = 0.021)+
  geom_text_repel(data=data_scores_cor,
                  aes(label=sample),show.legend=FALSE,color='black')+
  theme_minimal()

# Second plot
plot2 = ggplot(data_scores1, aes(x = `Score (MOFA)`, y = `Score (MFA)`)) + 
  geom_point(color = "dodgerblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -varphi.data(abs(data_scores_cor$`Score (MFA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_hline(yintercept = varphi.data(abs(data_scores_cor$`Score (MFA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_vline(xintercept = -varphi.data(abs(data_scores_cor$`Score (MOFA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_vline(xintercept = varphi.data(abs(data_scores_cor$`Score (MOFA)`)), linetype = "dashed", color = "#DCDDE0") +
  labs(x = "Score (MOFA)", y = "Score (MFA)") +
  geom_rect(aes(xmin = -varphi.data(abs(data_scores_cor$`Score (MOFA)`)) , xmax = -Inf, ymin = -varphi.data(abs(data_scores_cor$`Score (MFA)`)) , ymax = -Inf), fill = "#DCDDE0", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_scores_cor$`Score (MOFA)`)) , xmax = Inf, ymin =varphi.data(abs(data_scores_cor$`Score (MFA)`)) , ymax = Inf), fill = "#DCDDE0", alpha = 0.021)+
  geom_text_repel(data=data_scores_cor,
                  aes(label=sample),show.legend=FALSE,color='black')+
  theme_minimal()

# Third plot
plot3 = ggplot(data_scores1, aes(x = `Score (MOFA)`, y = `Score (GFA)`)) + 
  geom_point(color = "dodgerblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -varphi.data(abs(data_scores_cor$`Score (GFA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_hline(yintercept = varphi.data(abs(data_scores_cor$`Score (GFA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_vline(xintercept = -varphi.data(abs(data_scores_cor$`Score (MOFA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_vline(xintercept = varphi.data(abs(data_scores_cor$`Score (MOFA)`)), linetype = "dashed", color = "#DCDDE0") +
  labs(x = "Score (MOFA)", y = "Score (GFA)") +
  geom_rect(aes(xmin = -varphi.data(abs(data_scores_cor$`Score (MOFA)`)) , xmax = -Inf, ymin = -varphi.data(abs(data_scores_cor$`Score (GFA)`)) , ymax = -Inf), fill = "#DCDDE0", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_scores_cor$`Score (MOFA)`)) , xmax = Inf, ymin =varphi.data(abs(data_scores_cor$`Score (GFA)`)) , ymax = Inf), fill = "#DCDDE0", alpha = 0.021)+
  geom_text_repel(data=data_scores_cor,
                  aes(label=sample),show.legend=FALSE,color='black')+
  theme_minimal()

# Four plot
plot4 = ggplot(data_scores1, aes(x = `Score (FABIA)`, y = `Score (MFA)`)) + 
  geom_point(color = "dodgerblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -varphi.data(abs(data_scores_cor$`Score (MFA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_hline(yintercept = varphi.data(abs(data_scores_cor$`Score (MFA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_vline(xintercept = -varphi.data(abs(data_scores_cor$`Score (FABIA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_vline(xintercept = varphi.data(abs(data_scores_cor$`Score (FABIA)`)), linetype = "dashed", color = "#DCDDE0") +
  labs(x = "Score (FABIA)", y = "Score (MFA)") +
  geom_rect(aes(xmin = -varphi.data(abs(data_scores_cor$`Score (FABIA)`)) , xmax = -Inf, ymin = -varphi.data(abs(data_scores_cor$`Score (MFA)`)) , ymax = -Inf), fill = "#DCDDE0", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_scores_cor$`Score (FABIA)`)) , xmax = Inf, ymin =varphi.data(abs(data_scores_cor$`Score (MFA)`)) , ymax = Inf), fill = "#DCDDE0", alpha = 0.021)+
  geom_text_repel(data=data_scores_cor,
                  aes(label=sample),show.legend=FALSE,color='black')+
  theme_minimal()

# Five plot
plot5 = ggplot(data_scores1, aes(x = `Score (FABIA)`, y = `Score (GFA)`)) + 
  geom_point(color = "#0078BD") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#F6F6F6") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#F6F6F6") +
  geom_hline(yintercept = -varphi.data(abs(data_scores_cor$`Score (GFA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_hline(yintercept = varphi.data(abs(data_scores_cor$`Score (GFA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_vline(xintercept = -varphi.data(abs(data_scores_cor$`Score (FABIA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_vline(xintercept = varphi.data(abs(data_scores_cor$`Score (FABIA)`)), linetype = "dashed", color = "#DCDDE0") +
  labs(x = "Score (FABIA)", y = "Score (GFA)") +
  geom_rect(aes(xmin = -varphi.data(abs(data_scores_cor$`Score (FABIA)`)) , xmax = -Inf, ymin = -varphi.data(abs(data_scores_cor$`Score (GFA)`)) , ymax = -Inf), fill = "#DCDDE0", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_scores_cor$`Score (FABIA)`)) , xmax = Inf, ymin =varphi.data(abs(data_scores_cor$`Score (GFA)`)) , ymax = Inf), fill = "#DCDDE0", alpha = 0.021)+
  geom_text_repel(data=data_scores_cor,
                  aes(label=sample),show.legend=FALSE,color='black')+
  theme_minimal()

# Six plot
plot6 = ggplot(data_scores1, aes(x = `Score (MFA)`, y = `Score (GFA)`)) + 
  geom_point(color = "#0078BD") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#F6F6F6") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#F6F6F6") +
  geom_hline(yintercept = -varphi.data(abs(data_scores_cor$`Score (GFA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_hline(yintercept = varphi.data(abs(data_scores_cor$`Score (GFA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_vline(xintercept = -varphi.data(abs(data_scores_cor$`Score (MFA)`)) , linetype = "dashed", color = "#DCDDE0") +
  geom_vline(xintercept = varphi.data(abs(data_scores_cor$`Score (MFA)`)), linetype = "dashed", color = "#DCDDE0") +
  labs(x = "Score (MFA)", y = "Score (GFA)") +
  geom_rect(aes(xmin = -varphi.data(abs(data_scores_cor$`Score (MFA)`)) , xmax = -Inf, ymin = -varphi.data(abs(data_scores_cor$`Score (GFA)`)) , ymax = -Inf), fill = "#DCDDE0", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_scores_cor$`Score (MFA)`)) , xmax = Inf, ymin =varphi.data(abs(data_scores_cor$`Score (GFA)`)) , ymax = Inf), fill = "#DCDDE0", alpha = 0.021)+
  geom_text_repel(data=data_scores_cor,
                  aes(label=sample),show.legend=FALSE,color='black')+
  theme_minimal()

# Arrange plots side by side
grid.arrange(plot1, plot2, plot3, plot4, plot5,plot6, ncol = 3, nrow = 2)

# ------------------------------ JACCARD INDEX ------------------------------------
# Jaccard index
jaccard.similarity.index <- function(data){
  # Extract index columns
  index_columns <- grep("_index$", names(data), value = TRUE)
  # Initialize an empty dataframe with a single row of NA values
  jaccard_df <- data.frame(
    Column1 = NA_character_,
    Column2 = NA_character_,
    similar_count = NA_real_,
    dissimilar_count = NA_real_,
    sum_column1 = NA_real_,
    sum_column2 = NA_real_,
    jaccard_index = NA_real_,
    stringsAsFactors = FALSE
  )
  jaccard_df <- jaccard_df[0, ]  # Remove the initial row
  # Calculate Jaccard similarity for each pair of index columns
  for (i in 1:(length(index_columns) - 1)) {
    for (j in (i + 1):length(index_columns)) {
      column1 <- index_columns[i]
      column2 <- index_columns[j]
      if (length(data[[column1]]) != length(data[[column2]])) {
        stop("Columns must have the same length.")
      }
      # Calculate Jaccard similarity index
      similar_count <- sum(data[[column1]] == 1 & data[[column2]] == 1)
      dissimilar_count <- sum(data[[column1]] != data[[column2]])
      # Modify this part to correctly calculate 
      sum_column1 <- sum(data[[column1]] == 1)
      sum_column2 <- sum(data[[column2]] == 1)
      #total_sum
      jaccard_index <- similar_count / (sum_column1 + sum_column2 - similar_count)
      ji = ifelse(jaccard_index == -1, 1, jaccard_index)
      # Append the results to the dataframe
      ji <- round(ji, 2)
      jaccard_df <- rbind(jaccard_df, c(column1, column2, ji, similar_count, dissimilar_count, sum_column1, sum_column2))
    }
  }
  # Rename the columns in the final result
  colnames(jaccard_df) <- c("Column1", "Column2", "ji", "similar_count", "dissimilar_count", "sum_column1", "sum_column2")
  return(jaccard_df)
}

# Install necessary package
if (!require(ggVennDiagram)) install.packages("ggVennDiagram")

# Load required libraries
library(ggVennDiagram)


# Define your list of binary vectors
x = list(
  A = data_scores_cor_index$MFA_score_index,  # FABIA set
  B = data_scores_cor_index$GFA_score_index # MOFA set
)

# Convert the binary vectors into positions for the Venn diagram
x <- list(
  A = which(x$A == 1), 
  B = which(x$B == 1)
)

# Generate the 4D Venn diagram
ggVennDiagram(x, relative_width = 0.1) + 
  scale_fill_gradient(low = "#BBBFBF", high = "#41709C")

# ---------------------------------- Overall ------------------------------------

# Define your list of binary vectors
x = list(
  A = data_scores_cor_index$FABIA_score_index,  # FABIA set
  B = data_scores_cor_index$MOFA_score_index, # MOFA set
  C = data_scores_cor_index$MFA_score_index,  # MFA set,
  D = data_scores_cor_index$GFA_score_index  # GFA set
)

# Convert the binary vectors into positions for the Venn diagram
x <- list(
  A = which(x$A == 1), 
  B = which(x$B == 1),
  C = which(x$C == 1),
  D = which(x$D == 1)
)

# Generate the 4D Venn diagram
ggVennDiagram(x, relative_width = 0.1) + 
  scale_fill_gradient(low = "#BBBFBF", high = "#41709C")

# --------------------------------------- FEATURES ------------------------------------------

# Merging Weights/loading data for FABIA, MOFA, MFA, GFA
merged_fabia_mofa_loading <- merge(loading_fabia2, loading_mofa2, by = 'feature', all = TRUE) # Merge loading_fabia2 and loading_mofa2
merged_fabia_mofa_mfa_loading <- merge(merged_fabia_mofa_loading, loading_mfa2, by = 'feature', all = TRUE) # Merge merged_fabia_mofa_loading with loading_mfa2
data_weights <- merge(merged_fabia_mofa_mfa_loading, loading_gfa, by = 'feature', all = TRUE) # Merge merged_fabia_mofa_mfa with loading_gfa

data_weights <- merge(merge(merge(loading_fabia2, loading_mofa2, by = 'feature', all = TRUE), loading_mfa2, by = 'feature', all = TRUE), loading_gfa, by = 'feature', all = TRUE) # Merge merged_fabia_mofa_mfa with factor_gfa

# Weights Datasets
data_weights 

# Visuaization of scores using correlation matrix plot

data_weights_cor <- data_weights
ggpairs(data_weights_cor[c(2,4,6,7)])

# Define new variable names
new_var_names <- c("Loadings (FABIA)", "Loadings (MOFA)", "Loadings (MFA)", "Loadings (GFA)")
colnames(data_weights_cor)[c(2,4,6,7)] <- new_var_names
pairs.panels(data_weights_cor[c(2,4,6,7)], hist.col="lightgoldenrod", breaks = 16, pch = 8, cex = 0.5, ellipses=TRUE, lm=TRUE)

# ------------------------------------- VISUALIZATION OF WEIGHTS -------------------------------------------

data_weights_sub <- data_weights[c(1,2,3,4,6,7)]
new_var_names_sub <- c("feature", "FABIA", "view", "MOFA", "MFA", "GFA")
colnames(data_weights_sub) <- new_var_names_sub

# scaling
scale_data_weights_sub <- data_weights_sub[c(1,2,9,11,12,13)]
scale_data_weights_sub$`Loadings (FABIA)` = scale_matrix(scale_data_weights_sub$`Loadings (FABIA)`)
scale_data_weights_sub$`Loadings (MOFA)` = scale_matrix(scale_data_weights_sub$`Loadings (MOFA)`)
scale_data_weights_sub$`Loadings (MFA)` = scale_matrix(scale_data_weights_sub$`Loadings (MFA)`)
scale_data_weights_sub$`Loadings (GFA)` = scale_matrix(scale_data_weights_sub$`Loadings (GFA)`)

# Filter the data for only mRNA
data_weights.methylation <- data_weights_sub %>%
  filter(view == "Methylation")
data_weights.drugs <- data_weights_sub %>%
  filter(view == "Drugs")

# Correlation plot
pairs.panels(data_weights.methylation[c(2,4,5,6)], hist.col="#4D7CA8", breaks = 16, 
             cex.labels = 1.5, cex.axis = 1.2,pch = 8, cex = 0.55, ellipses=FALSE, lm=TRUE)

pairs.panels(data_weights.drugs[c(2,4,5,6)], hist.col="#DEA363", breaks = 16, 
             cex.labels = 1.5, cex.axis = 1.2,pch = 8, cex = 0.55, ellipses=FALSE, lm=TRUE)


plot1 = ggplot(data_weights.methylation, aes(x = MOFA, y = MFA)) + 
  geom_point(color = "blue4") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey95") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey95") +
  geom_hline(yintercept = -varphi.data(abs(data_weights.methylation$MFA)) , linetype = "dashed", color = "grey100") +
  geom_hline(yintercept = varphi.data(abs(data_weights.methylation$MFA)) , linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = -varphi.data(abs(data_weights.methylation$MOFA)), linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = varphi.data(abs(data_weights.methylation$MOFA)), linetype = "dashed", color = "grey100") +
  labs(x = "Loading (MOFA)", y = "Loading (MFA)") +
  # Add shading beyond the vertical line
  geom_rect(aes(xmin = -varphi.data(abs(data_weights.methylation$MOFA)), xmax = -Inf, ymin =-varphi.data(abs(data_weights.methylation$MFA))  , ymax = -Inf), fill = "grey100", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_weights.methylation$MOFA)), xmax = Inf, ymin =varphi.data(abs(data_weights.methylation$MFA))  , ymax = Inf), fill = "grey100", alpha = 0.021)+
  theme_minimal()

plot2 = ggplot(data_weights.drugs, aes(x = MOFA, y = MFA)) + 
  geom_point(color = "grey1") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey95") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey95") +
  geom_hline(yintercept = -varphi.data(abs(data_weights.drugs$MFA)) , linetype = "dashed", color = "grey100") +
  geom_hline(yintercept = varphi.data(abs(data_weights.drugs$MFA)) , linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = -varphi.data(abs(data_weights.drugs$MOFA)), linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = varphi.data(abs(data_weights.drugs$MOFA)), linetype = "dashed", color = "grey100") +
  labs(x = "Loading (MOFA)", y = "Loading (MFA)") +
  # Add shading beyond the vertical line
  geom_rect(aes(xmin = -varphi.data(abs(data_weights.drugs$MOFA)), xmax = -Inf, ymin =-varphi.data(abs(data_weights.drugs$MFA))  , ymax = -Inf), fill = "grey100", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_weights.drugs$MOFA)), xmax = Inf, ymin =varphi.data(abs(data_weights.drugs$MFA))  , ymax = Inf), fill = "grey100", alpha = 0.021)+
  theme_minimal()


# Arrange plots side by side
grid.arrange(plot1, plot2, ncol = 2)

# MOFA VS FABIA

plot1 = ggplot(data_weights.methylation, aes(x = MOFA, y = FABIA)) + 
  geom_point(color = "blue4") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey95") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey95") +
  geom_hline(yintercept = -varphi.data(abs(data_weights.methylation$FABIA)) , linetype = "dashed", color = "grey100") +
  geom_hline(yintercept = varphi.data(abs(data_weights.methylation$FABIA)) , linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = -varphi.data(abs(data_weights.methylation$MOFA)), linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = varphi.data(abs(data_weights.methylation$MOFA)), linetype = "dashed", color = "grey100") +
  labs(x = "Loading (MOFA)", y = "Loading (FABIA)") +
  # Add shading beyond the vertical line
  geom_rect(aes(xmin = -varphi.data(abs(data_weights.methylation$MOFA)), xmax = -Inf, ymin =-varphi.data(abs(data_weights.methylation$FABIA))  , ymax = -Inf), fill = "grey100", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_weights.methylation$MOFA)), xmax = Inf, ymin =varphi.data(abs(data_weights.methylation$FABIA))  , ymax = Inf), fill = "grey100", alpha = 0.021)+
  theme_minimal()

plot2 = ggplot(data_weights.drugs, aes(x = MOFA, y = FABIA)) + 
  geom_point(color = "grey1") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey95") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey95") +
  geom_hline(yintercept = -varphi.data(abs(data_weights.drugs$FABIA)) , linetype = "dashed", color = "grey100") +
  geom_hline(yintercept = varphi.data(abs(data_weights.drugs$FABIA)) , linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = -varphi.data(abs(data_weights.drugs$MOFA)), linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = varphi.data(abs(data_weights.drugs$MOFA)), linetype = "dashed", color = "grey100") +
  labs(x = "Loading (MOFA)", y = "Loading (FABIA)") +
  # Add shading beyond the vertical line
  geom_rect(aes(xmin = -varphi.data(abs(data_weights.drugs$MOFA)), xmax = -Inf, ymin =-varphi.data(abs(data_weights.drugs$FABIA))  , ymax = -Inf), fill = "grey100", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_weights.drugs$MOFA)), xmax = Inf, ymin =varphi.data(abs(data_weights.drugs$FABIA))  , ymax = Inf), fill = "grey100", alpha = 0.021)+
  theme_minimal()

# Arrange plots side by side
grid.arrange(plot1, plot2, ncol = 2)

plot1 = ggplot(data_weights.methylation, aes(x = MOFA, y = GFA)) + 
  geom_point(color = "blue4") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey95") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey95") +
  geom_hline(yintercept = -varphi.data(abs(data_weights.methylation$GFA)) , linetype = "dashed", color = "grey100") +
  geom_hline(yintercept = varphi.data(abs(data_weights.methylation$GFA)) , linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = -varphi.data(abs(data_weights.methylation$MOFA)), linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = varphi.data(abs(data_weights.methylation$MOFA)), linetype = "dashed", color = "grey100") +
  labs(x = "Loading (MOFA)", y = "Loading (GFA)") +
  # Add shading beyond the vertical line
  geom_rect(aes(xmin = -varphi.data(abs(data_weights.methylation$MOFA)), xmax = -Inf, ymin =-varphi.data(abs(data_weights.methylation$GFA))  , ymax = -Inf), fill = "grey100", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_weights.methylation$MOFA)), xmax = Inf, ymin =varphi.data(abs(data_weights.methylation$GFA))  , ymax = Inf), fill = "grey100", alpha = 0.021)+
  theme_minimal()

plot2 = ggplot(data_weights.drugs, aes(x = MOFA, y = GFA)) + 
  geom_point(color = "grey1") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey95") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey95") +
  geom_hline(yintercept = -varphi.data(abs(data_weights.drugs$GFA)) , linetype = "dashed", color = "grey100") +
  geom_hline(yintercept = varphi.data(abs(data_weights.drugs$GFA)) , linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = -varphi.data(abs(data_weights.drugs$MOFA)), linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = varphi.data(abs(data_weights.drugs$MOFA)), linetype = "dashed", color = "grey100") +
  labs(x = "Loading (MOFA)", y = "Loading (GFA)") +
  # Add shading beyond the vertical line
  geom_rect(aes(xmin = -varphi.data(abs(data_weights.drugs$MOFA)), xmax = -Inf, ymin =-varphi.data(abs(data_weights.drugs$GFA))  , ymax = -Inf), fill = "grey100", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_weights.drugs$MOFA)), xmax = Inf, ymin =varphi.data(abs(data_weights.drugs$GFA))  , ymax = Inf), fill = "grey100", alpha = 0.021)+
  theme_minimal()


# Arrange plots side by side
grid.arrange(plot1, plot2, ncol = 2)

# FABIA VS MFA
plot1 = ggplot(data_weights.methylation, aes(x = FABIA, y = MFA)) + 
  geom_point(color = "blue4") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey95") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey95") +
  geom_hline(yintercept = -varphi.data(abs(data_weights.methylation$MFA)) , linetype = "dashed", color = "grey100") +
  geom_hline(yintercept = varphi.data(abs(data_weights.methylation$MFA)) , linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = -varphi.data(abs(data_weights.methylation$FABIA)), linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = varphi.data(abs(data_weights.methylation$FABIA)), linetype = "dashed", color = "grey100") +
  labs(x = "Loading (FABIA)", y = "Loading (MFA)") +
  # Add shading beyond the vertical line
  geom_rect(aes(xmin = -varphi.data(abs(data_weights.methylation$FABIA)), xmax = -Inf, ymin =-varphi.data(abs(data_weights.methylation$MFA))  , ymax = -Inf), fill = "grey100", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_weights.methylation$FABIA)), xmax = Inf, ymin =varphi.data(abs(data_weights.methylation$MFA))  , ymax = Inf), fill = "grey100", alpha = 0.021)+
  theme_minimal()

plot2 = ggplot(data_weights.drugs, aes(x = FABIA, y = MFA)) + 
  geom_point(color = "grey1") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey95") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey95") +
  geom_hline(yintercept = -varphi.data(abs(data_weights.drugs$MFA)) , linetype = "dashed", color = "grey100") +
  geom_hline(yintercept = varphi.data(abs(data_weights.drugs$MFA)) , linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = -varphi.data(abs(data_weights.drugs$FABIA)), linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = varphi.data(abs(data_weights.drugs$FABIA)), linetype = "dashed", color = "grey100") +
  labs(x = "Loading (FABIA)", y = "Loading (MFA)") +
  # Add shading beyond the vertical line
  geom_rect(aes(xmin = -varphi.data(abs(data_weights.drugs$FABIA)), xmax = -Inf, ymin =-varphi.data(abs(data_weights.drugs$MFA))  , ymax = -Inf), fill = "grey100", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_weights.drugs$FABIA)), xmax = Inf, ymin =varphi.data(abs(data_weights.drugs$MFA))  , ymax = Inf), fill = "grey100", alpha = 0.021)+
  theme_minimal()


# Arrange plots side by side
grid.arrange(plot1, plot2, ncol = 2)

# FABIA VS GFA
plot1 = ggplot(data_weights.methylation, aes(x = FABIA, y = MFA)) + 
  geom_point(color = "blue4") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey95") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey95") +
  geom_hline(yintercept = -varphi.data(abs(data_weights.methylation$MFA)) , linetype = "dashed", color = "grey100") +
  geom_hline(yintercept = varphi.data(abs(data_weights.methylation$MFA)) , linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = -varphi.data(abs(data_weights.methylation$FABIA)), linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = varphi.data(abs(data_weights.methylation$FABIA)), linetype = "dashed", color = "grey100") +
  labs(x = "Loading (FABIA)", y = "Loading (GFA)") +
  # Add shading beyond the vertical line
  geom_rect(aes(xmin = -varphi.data(abs(data_weights.methylation$FABIA)), xmax = -Inf, ymin =-varphi.data(abs(data_weights.methylation$MFA))  , ymax = -Inf), fill = "grey100", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_weights.methylation$FABIA)), xmax = Inf, ymin =varphi.data(abs(data_weights.methylation$MFA))  , ymax = Inf), fill = "grey100", alpha = 0.021)+
  theme_minimal()

plot2 = ggplot(data_weights.drugs, aes(x = FABIA, y = GFA)) + 
  geom_point(color = "grey1") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey95") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey95") +
  geom_hline(yintercept = -varphi.data(abs(data_weights.drugs$GFA)) , linetype = "dashed", color = "grey100") +
  geom_hline(yintercept = varphi.data(abs(data_weights.drugs$GFA)) , linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = -varphi.data(abs(data_weights.drugs$FABIA)), linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = varphi.data(abs(data_weights.drugs$FABIA)), linetype = "dashed", color = "grey100") +
  labs(x = "Loading (FABIA)", y = "Loading (GFA)") +
  # Add shading beyond the vertical line
  geom_rect(aes(xmin = -varphi.data(abs(data_weights.drugs$FABIA)), xmax = -Inf, ymin =-varphi.data(abs(data_weights.drugs$GFA))  , ymax = -Inf), fill = "grey100", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_weights.drugs$FABIA)), xmax = Inf, ymin =varphi.data(abs(data_weights.drugs$GFA))  , ymax = Inf), fill = "grey100", alpha = 0.021)+
  theme_minimal()


# Arrange plots side by side
grid.arrange(plot1, plot2, ncol = 2)

# MFA VS GFA
plot1 = ggplot(data_weights.methylation, aes(x = MFA, y = GFA)) + 
  geom_point(color = "blue4") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey95") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey95") +
  geom_hline(yintercept = -varphi.data(abs(data_weights.methylation$GFA)) , linetype = "dashed", color = "grey100") +
  geom_hline(yintercept = varphi.data(abs(data_weights.methylation$GFA)) , linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = -varphi.data(abs(data_weights.methylation$MFA)), linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = varphi.data(abs(data_weights.methylation$MFA)), linetype = "dashed", color = "grey100") +
  labs(x = "Loading (MFA)", y = "Loading (GFA)") +
  # Add shading beyond the vertical line
  geom_rect(aes(xmin = -varphi.data(abs(data_weights.methylation$MFA)), xmax = -Inf, ymin =-varphi.data(abs(data_weights.methylation$MFA))  , ymax = -Inf), fill = "grey100", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_weights.methylation$MFA)), xmax = Inf, ymin =varphi.data(abs(data_weights.methylation$MFA))  , ymax = Inf), fill = "grey100", alpha = 0.021)+
  theme_minimal()

plot2 = ggplot(data_weights.drugs, aes(x = MFA, y = GFA)) + 
  geom_point(color = "grey1") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey95") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey95") +
  geom_hline(yintercept = -varphi.data(abs(data_weights.drugs$GFA)) , linetype = "dashed", color = "grey100") +
  geom_hline(yintercept = varphi.data(abs(data_weights.drugs$GFA)) , linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = -varphi.data(abs(data_weights.drugs$MFA)), linetype = "dashed", color = "grey100") +
  geom_vline(xintercept = varphi.data(abs(data_weights.drugs$MFA)), linetype = "dashed", color = "grey100") +
  labs(x = "Loading (MFA)", y = "Loading (GFA)") +
  # Add shading beyond the vertical line
  geom_rect(aes(xmin = -varphi.data(abs(data_weights.drugs$MFA)), xmax = -Inf, ymin =-varphi.data(abs(data_weights.drugs$GFA))  , ymax = -Inf), fill = "grey100", alpha = 0.021)+
  geom_rect(aes(xmin = varphi.data(abs(data_weights.drugs$MFA)), xmax = Inf, ymin =varphi.data(abs(data_weights.drugs$GFA))  , ymax = Inf), fill = "grey100", alpha = 0.021)+
  theme_minimal()


# Arrange plots side by side
grid.arrange(plot1, plot2, ncol = 2)

# ------------------ JACCARD INDEX ---------------------------
# Jaccard index for weights
ji_data_weights_cor_methyl <- data_weights.methylation
data_weights_cor_methyl_index <- calculate_threshold_and_index(ji_data_weights_cor_methyl)
weights_index_columns_methyl.cll <- grep("_index$", names(data_weights_cor_methyl_index), value = TRUE)

ji_data_weights_cor_drugs <- data_weights.drugs
data_weights_cor_drugs_index <- calculate_threshold_and_index(ji_data_weights_cor_drugs)
weights_index_columns_drugs.cll <- grep("_index$", names(data_weights_cor_drugs_index), value = TRUE)

jaccard_weights_data_methyl.cll = jaccard.similarity.index(data_weights_cor_methyl_index)
jaccard_weights_data_drugs.cll = jaccard.similarity.index(data_weights_cor_drugs_index)

#
# Install necessary package
if (!require(ggVennDiagram)) install.packages("ggVennDiagram")

# Load required libraries
library(ggVennDiagram)


# Define your list of binary vectors
x = list(
  A = data_weights_cor_drugs_index$FABIA_index,  # FABIA set
  B = data_weights_cor_drugs_index$MOFA_index # MOFA set
)

# Convert the binary vectors into positions for the Venn diagram
x <- list(
  A = which(x$A == 1), 
  B = which(x$B == 1)
)

# Generate the 4D Venn diagram
ggVennDiagram(x, relative_width = 0.1) + 
  scale_fill_gradient(low = "#BBBFBF", high = "#41709C")

# ---------------------------------- Overall ------------------------------------

# Define your list of binary vectors
x = list(
  A = data_weights_cor_drugs_index$FABIA_index,  # FABIA set
  B = data_weights_cor_drugs_index$MOFA_index, # MOFA set
  C = data_weights_cor_drugs_index$MFA_index,  # MFA set,
  D = data_weights_cor_drugs_index$GFA_index  # GFA set
)

# Convert the binary vectors into positions for the Venn diagram
x <- list(
  A = which(x$A == 1), 
  B = which(x$B == 1),
  C = which(x$C == 1),
  D = which(x$D == 1)
)

# Generate the 4D Venn diagram
ggVennDiagram(x, relative_width = 0.1) + 
  scale_fill_gradient(low = "#BBBFBF", high = "#41709C")




















# Cosine similarity
#install.packages("lsa")
library(lsa)
cosine(mRNA_data_17$loading_FABIA, mRNA_data_17$loading_GFA)
cosine(mRNA_data_17$loading_FABIA, mRNA_data_17$loading_MFA)
cosine(mRNA_data_17$loading_FABIA, mRNA_data_17$loading_MOFA)
cosine(mRNA_data_17$loading_MFA, mRNA_data_17$loading_MOFA)
cosine(mRNA_data_17$loading_GFA, mRNA_data_17$loading_MOFA)

cosine(Protein_data_17$loading_FABIA, Protein_data_17$loading_GFA)
cosine(Protein_data_17$loading_FABIA, Protein_data_17$loading_MFA)
cosine(Protein_data_17$loading_FABIA, Protein_data_17$loading_MOFA)
cosine(Protein_data_17$loading_MFA, Protein_data_17$loading_MOFA)
cosine(Protein_data_17$loading_GFA, Protein_data_17$loading_MOFA)

#######################################################################################
# Visualization
mRNA_data_17$ID = row_number(mRNA_data_17)
ggplot(mRNA_data_17, aes(x=ID)) + 
  geom_point(aes(y = loading_FABIA), color = "darkred") + 
  geom_point(aes(y = loading_MFA), color="steelblue")  +
  geom_point(aes(y = (-1)*loading_MOFA), color="blue")+
  geom_point(aes(y = loading_GFA), color="yellow3")+
  labs(x = "Features", y = "Loading") +
  theme_bw()

# ------------------------------------------ END VISUALIZATION --------------------------------------------
# Visuaization of weights using correlation matrix plot
library(GGally); library(psych)

data_weights_cor <- data_weights
ggpairs(data_weights_cor[c(2,4,6,7)])


# Visualize scores using scatter plot
# Plot
data_scores_cor$IGHV <- as.factor(data_scores_cor$IGHV)
ggplot(data_scores_cor, aes(x = sample, y = `Score (GFA)`)) +
  geom_point(aes(fill = IGHV), size = 4, shape = 21) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  ggtitle("GFA: F 1") +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = c("0" = "yellow","1" = "steelblue1", "NA" = "grey")) +  # Manually assign colors to Condition
  labs(x = "Samples", y = "Score")

# Top samples
# Sort the values and get the top 50 values
sorted_values <- sort(abs(data_scores_cor$`Score (FABIA)`), decreasing = TRUE)
top_50_values <- head(sorted_values, n = 15)

# Identify the corresponding feature names for the top 50 values
top_50_samples <- data_scores_cor$sample[match(top_50_values, abs(data_scores_cor$`Score (FABIA)`))]

print(top_50_samples)

# Visualize weights/loading using scatter plot

# Filter the data for only View=drugs
data_weights_cor_drugs <- data_weights_cor %>%
  filter(view == "Drugs")

# Top drugs
# Sort the values and get the top 50 values
sorted_values <- sort(abs(data_weights_cor_drugs$`Loading (FABIA)`), decreasing = TRUE)
top_50_values <- head(sorted_values, n = 15)

N =100
# Sort the values and get the top 50 values
sorted_values <- sort(abs(data_weights_cor_drugs$`Loading (FABIA)`), decreasing = TRUE)
top_50_values <- head(sorted_values, n = N)
top_50_features_FABIA <- data_weights_cor_drugs$feature[match(top_50_values, abs(data_weights_cor_drugs$`Loading (FABIA)`))]

sorted_values <- sort(abs(data_weights_cor_drugs$`Loading (MOFA)`), decreasing = TRUE)
top_50_values <- head(sorted_values, n = N)
top_50_features_MOFA <- data_weights_cor_drugs$feature[match(top_50_values, abs(data_weights_cor_drugs$`Loading (MOFA)`))]

sorted_values <- sort(abs(data_weights_cor_drugs$`Loading (MFA)`), decreasing = TRUE)
top_50_values <- head(sorted_values, n = N)
top_50_features_MFA <- data_weights_cor_drugs$feature[match(top_50_values, abs(data_weights_cor_drugs$`Loading (MFA)`))]

sorted_values <- sort(abs(data_weights_cor_drugs$`Loading (GFA)`), decreasing = TRUE)
top_50_values <- head(sorted_values, n = N)
top_50_features_GFA <- data_weights_cor_drugs$feature[match(top_50_values, abs(data_weights_cor_drugs$`Loading (GFA)`))]

# Intersection
find_intersection <- function(vec1, vec2, vec3, vec4) {
  # Find intersection of vec1 and vec2
  intersection_vec <- intersect(vec1, vec2)
  
  # Find intersection of intersection_vec and vec3
  intersection_vec <- intersect(intersection_vec, vec3)
  
  # Find intersection of intersection_vec and vec4
  intersection_vec <- intersect(intersection_vec, vec4)
  
  return(intersection_vec)
}

result_intersection <- find_intersection(top_50_features_FABIA,top_50_features_MOFA, top_50_features_GFA, top_50_features_MFA)
print(result_intersection)

# Identify the corresponding feature names for the top 50 values
top_50_features <- data_weights_cor_drugs$feature[match(top_50_values, abs(data_weights_cor_drugs$`Loading (FABIA)`))]

print(top_50_features)

# Filter the data for only View=methylation
data_weights_cor_meth <- data_weights_cor %>%
  filter(view == "Methylation")

# Top methylation
N =200
# Sort the values and get the top 50 values
sorted_values <- sort(abs(data_weights_cor_meth$loading_FABIA1), decreasing = TRUE)
top_50_values <- head(sorted_values, n = N)
top_50_features_FABIA <- data_weights_cor_meth$feature[match(top_50_values, abs(data_weights_cor_meth$loading_FABIA1))]

sorted_values <- sort(abs(data_weights_cor_meth$loading_MOFA1), decreasing = TRUE)
top_50_values <- head(sorted_values, n = N)
top_50_features_MOFA <- data_weights_cor_meth$feature[match(top_50_values, abs(data_weights_cor_meth$loading_MOFA1))]

sorted_values <- sort(abs(data_weights_cor_meth$loading_MFA1), decreasing = TRUE)
top_50_values <- head(sorted_values, n = N)
top_50_features_MFA <- data_weights_cor_meth$feature[match(top_50_values, abs(data_weights_cor_meth$loading_MFA1))]

sorted_values <- sort(abs(data_weights_cor_meth$loading_GFA1), decreasing = TRUE)
top_50_values <- head(sorted_values, n = N)
top_50_features_GFA <- data_weights_cor_meth$feature[match(top_50_values, abs(data_weights_cor_meth$loading_GFA1))]

# Intersection
find_intersection <- function(vec1, vec2, vec3){#, vec4) {
  # Find intersection of vec1 and vec2
  intersection_vec <- intersect(vec1, vec2)
  
  # Find intersection of intersection_vec and vec3
  intersection_vec <- intersect(intersection_vec, vec3)
  
  # Find intersection of intersection_vec and vec4
  #intersection_vec <- intersect(intersection_vec, vec4)
  
  return(intersection_vec)
}

result_intersection <- find_intersection(top_50_features_MOFA, top_50_features_MFA, top_50_features_FABIA)
print(result_intersection)

# Intersection
find_intersection <- function(vec1, vec2, vec3, vec4) {
  # Find intersection of vec1 and vec2
  intersection_vec <- intersect(vec1, vec2)
  
  # Find intersection of intersection_vec and vec3
  intersection_vec <- intersect(intersection_vec, vec3)
  
  # Find intersection of intersection_vec and vec4
  intersection_vec <- intersect(intersection_vec, vec4)
  
  return(intersection_vec)
}

result_intersection <- find_intersection(top_50_features_FABIA, top_50_features_MOFA, top_50_features_MFA, top_50_features_GFA)
print(result_intersection)




# Correlation plots
# Define new variable names
new_var_names <- c("Loading (FABIA)", "Loading (MOFA)", "Loading (MFA)", "Loading (GFA)")
colnames(data_weights_cor_meth)[c(2,4,6,7)] <- new_var_names
pairs.panels(data_weights_cor_meth[c(2,4,6,7)], hist.col="deepskyblue4", breaks = 16, pch = 8, cex = 0.5, ellipses=TRUE, lm=TRUE)

# Plot using a single palette
plot1 = ggplot(data_weights_cor_drugs, aes(x = `Loading (MFA)`, y = `Loading (GFA)`)) + 
  geom_point(color = "dodgerblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Loading (MFA)", y = "Loading (GFA)") +
  ggtitle("Drugs") +
  theme(
    axis.ticks.x = element_blank(),
    plot.title = element_text(hjust = 1))+
  theme_minimal()

# Second plot
plot2 = ggplot(data_weights_cor_meth, aes(x = `Loading (MFA)`, y = `Loading (GFA)`)) + 
  geom_point(color = "deepskyblue4") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Loading (MFA)", y = "Loading (MFA)") +
  ggtitle("Methylation") +
  theme_minimal()

# Arrange plots side by side
grid.arrange(plot1, plot2, ncol = 2)

# Functions 
# Threshold calculator
calculate_threshold_and_index <- function(data) {
  for (i in seq_along(data)) {
    if (is.numeric(data[[i]])) {
      threshold <- varphi.data(data[[i]])
      
      # Create corresponding index column name
      index_column <- paste0(names(data)[i], "_index")
      
      # Create index column based on the threshold
      data[[index_column]] <- ifelse(data[[i]] > threshold, 1, 0)
    }
  }
  return(data)
}

# Jaccard Similarity Index
varphi.data <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    varphi <- quantile(numeric_values, probs = 0.90)
    return(varphi)
  } else {
    return(NA)  # Return NA if not all values are numeric
  }
}

# Function to re-scale a vector to (-1,1) scale
rescale_minusone_to_one <- function(vector){
  scaled_vector <- vector / max(abs(vector), 1e-6) # Using 1e-6 to avoid division by zero
  return(scaled_vector)
}

# subset data to only scores
data_scores_cor_subset <- data_scores_cor[,c(2,11,13,14)]

# Obtain thresholds
data_scores_cor_subset.i <- data_scores_cor_subset %>%
  mutate_if(is.numeric, ~ rescale_minusone_to_one(.))

data_scores_cor_subset.a <- data_scores_cor_subset.i %>%
  mutate_if(is.numeric, ~ abs(.))

varphi.data(abs(data_scores_cor_subset.a$`Score (FABIA)`))
varphi.data(abs(data_scores_cor_subset.a$`Score (MOFA)`))
varphi.data(abs(data_scores_cor_subset.a$`Score (MFA)`))
varphi.data(abs(data_scores_cor_subset.a$`Score (GFA)`))

# calculate_threshold_and_index
data_scores_cor_subset.b <- calculate_threshold_and_index(data_scores_cor_subset.a)
data_scores_cor_index.b <- grep("_index$", names(data_scores_cor_subset.b), value = TRUE)

# Subset the dataframe based on selected columns
data_scores_cor_subset.c <- data_scores_cor_subset.b[, c(data_scores_cor_index.b)]
jaccard_scores = jaccard.index(data_scores_cor_subset.c)

# 2. Loading/Weight
# Jaccard Similarity Index
varphi.data <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    varphi <- quantile(numeric_values, probs = 0.60)
    return(varphi)
  } else {
    return(NA)  # Return NA if not all values are numeric
  }
}

# subset data to only loading
data_weights_cor_drugs_subset <- data_weights_cor_drugs[,c(2,4,6,7)]

# Obtain thresholds
data_weights_cor_drugs_subset.i <- data_weights_cor_drugs_subset %>%
  mutate_if(is.numeric, ~ rescale_minusone_to_one(.))

data_weights_cor_drugs_subset.a <- data_weights_cor_drugs_subset.i %>%
  mutate_if(is.numeric, ~ abs(.))

varphi.data(abs(data_weights_cor_drugs_subset.a$`Loading (FABIA)`))
varphi.data(abs(data_weights_cor_drugs_subset.a$`Loading (MOFA)`))
varphi.data(abs(data_weights_cor_drugs_subset.a$`Loading (MFA)`))
varphi.data(abs(data_weights_cor_drugs_subset.a$`Loading (GFA)`))

# calculate_threshold_and_index
data_weights_cor_drugs_subset.b <- calculate_threshold_and_index(data_weights_cor_drugs_subset.a)
data_weights_cor_drugs_index.b <- grep("_index$", names(data_weights_cor_drugs_subset.b), value = TRUE)

# Subset the dataframe based on selected columns
data_weights_cor_drugs_subset.c <- data_weights_cor_drugs_subset.b[, c(data_weights_cor_drugs_index.b)]
jaccard_scores_drugs = jaccard.index(data_weights_cor_drugs_subset.c)

# (ii) Methylation
varphi.data <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    varphi <- quantile(numeric_values, probs = 0.95)
    return(varphi)
  } else {
    return(NA)  # Return NA if not all values are numeric
  }
}
# subset data to only loading
data_weights_cor_meth_subset <- data_weights_cor_meth[,c(2,4,6,7)]

# Obtain thresholds
data_weights_cor_meth_subset.i <- data_weights_cor_meth_subset %>%
  mutate_if(is.numeric, ~ rescale_minusone_to_one(.))

data_weights_cor_meth_subset.a <- data_weights_cor_meth_subset.i %>%
  mutate_if(is.numeric, ~ abs(.))

varphi.data(abs(data_weights_cor_meth_subset.a$`Loading (FABIA)`)); varphi.data(abs(data_weights_cor_meth_subset.a$`Loading (MOFA)`)); varphi.data(abs(data_weights_cor_meth_subset.a$`Loading (MFA)`)); varphi.data(abs(data_weights_cor_meth_subset.a$`Loading (GFA)`))

# calculate_threshold_and_index
data_weights_cor_meth_subset.b <- calculate_threshold_and_index(data_weights_cor_meth_subset.a)
data_weights_cor_meth_index.b <- grep("_index$", names(data_weights_cor_meth_subset.b), value = TRUE)

# Subset the dataframe based on selected columns
data_weights_cor_meth_subset.c <- data_weights_cor_meth_subset.b[, c(data_weights_cor_meth_index.b)]
jaccard_scores_meth = jaccard.index(data_weights_cor_meth_subset.c)

# ----------------------------------- MICROBIOME MULTI-MODAL DATA --------------------------------------
# 1. Microbiome multi-modal data - Data already filtered and normalized

# Load libraries
library(data.table); library(purrr); library(ggplot2); library(ggpubr); library(MOFA2); library(tidyr)

# Load data
dt <- fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/microbiome/data.txt.gz") 
metadata <-fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/microbiome/metadata.txt.gz")

# Create the MOFA object and train the model
bacteria_dt_subset <- dt[dt$view =='Bacteria']
bacteria_dt_subset <- pivot_wider(bacteria_dt_subset, names_from = sample, values_from = value) # wider version
bacteria_dt_subset_df <- as.data.frame(bacteria_dt_subset)
rownames(bacteria_dt_subset_df) <- bacteria_dt_subset_df$feature; bacteria_dt_subset2<-bacteria_dt_subset_df[,-c(1,2)]

virus_dt_subset <- dt[dt$view =='Viruses']
virus_dt_subset <- pivot_wider(virus_dt_subset, names_from = sample, values_from = value) # wide version
virus_dt_subset_df <- as.data.frame(virus_dt_subset)
rownames(virus_dt_subset_df) <- virus_dt_subset_df$feature; virus_dt_subset2<-virus_dt_subset_df[,-c(1,2)]

#virus_dt_subset_cleaned <- virus_dt_subset %>%
#  mutate(across(everything(), ~replace(., is.na(.), 0)))

fungi_dt_subset <- dt[dt$view =='Fungi']
fungi_dt_subset <- pivot_wider(fungi_dt_subset, names_from = sample, values_from = value) # wide version
fungi_dt_subset_df <- as.data.frame(fungi_dt_subset)
rownames(fungi_dt_subset_df) <- fungi_dt_subset_df$feature; fungi_dt_subset2<-fungi_dt_subset_df[,-c(1,2)]

# Data
dt_subset <- rbind(bacteria_dt_subset_df, fungi_dt_subset_df)
df_clean <- dt_subset %>%
  select(where(~ !any(is.na(.))))

dt_subset_cleaned <- dt_subset %>%
  mutate(across(everything(), ~replace(., is.na(.), 0)))

dt_subset <- list(Bacteria = as.matrix(bacteria_dt_subset2), Fungi = as.matrix(fungi_dt_subset2))
MOFAobject <- create_mofa(dt_subset)
MOFAobject

# Plot data overview
plot_data_overview(MOFAobject)

# Define MOFA options
## Data options
data_opts <- get_default_data_options(MOFAobject)
data_opts

## Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 3
model_opts

## Training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
train_opts

## Train the MOFA model
MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts
)
outfile_object_mmd = paste0(getwd(),"model_object_mmd.hdf5")
MOFAobject.mmd.trained <- run_mofa(MOFAobject, outfile_object_mmd)

# Overview of the trained MOFA model
## Slots - The MOFA object consists of multiple slots where relevant data and information is stored. For descriptions, you can read the documentation using ?MOFA. The most important slots are:
# data: input data used to train the model (features are centered at zero mean)
# samples_metadata: sample metadata information
# expectations: expectations of the posterior distributions for the Weights and the Factors
slotNames(MOFAobject)

## Add sample metadata to the model
samples_metadata(MOFAobject) <- metadata

## Correlation between factors # A good sanity check is to verify that the Factors are largely uncorrelated. 
plot_factor_cor(MOFAobject.mmd.trained)

## Plot variance decomposition
### Variance decomposition by Factor
plot_variance_explained(MOFAobject.mmd.trained, max_r2=15)

### Total variance explained per view
plot_variance_explained(MOFAobject.mmd.trained, plot_total = T, factors = "all")[[2]]

# Factor scores
factor_mofa <- get_factors(MOFAobject.mmd.trained, factors = c(1,2), as.data.frame = T)
factor_mofa <- factor_mofa %>%
  spread(factor, value)
factor_mofa2 = factor_mofa[, !(colnames(factor_mofa) %in% c('group'))]
colnames(factor_mofa2) =  c('sample', 'score_MOFA1', 'score_MOFA2')

# Weights/Loading 
loading_mofa<- get_weights(MOFAobject.mmd.trained, factors = c(1,2), as.data.frame = T)
loading_mofa2 <- loading_mofa %>%
  spread(factor, value)
colnames(loading_mofa2) =  c('feature', 'view', 'loading_MOFA1', 'loading_MOFA2')

# FABIA

# merge the two datasets
mmd_fabia <- rbind(bacteria_dt_subset2, fungi_dt_subset2)

# FABIA MODEL
set.seed(123)
FABIAobject.mmd.trained <- fabia(mmd_fabia,
                                 p = 2, alpha = 0.01, cyc = 1000, spl = 0.5, spz = 0.5, random = 1.0,  
                                 center = 2, norm = 2, lap = 1.0, nL = 1 
)

# Factor scores
factor_fabia = data.frame(t(FABIAobject.mmd.trained@Z)) 
factor_fabia$sample = rownames(factor_fabia)
colnames(factor_fabia) <- c("score_FABIA1", "score_FABIA2", "sample")
factor_fabia2 <- merge(x=factor_fabia, y=metadata, by="sample", all = TRUE)

# Weights/Loading 
loading_fabia = data.frame(FABIAobject.mmd.trained@L)
loading_fabia$feature = rownames(loading_fabia)
colnames(loading_fabia) <- c("loading_FABIA1", "loading_FABIA2", "feature")
loading_fabia2 = loading_fabia

# MFA
set.seed(123)
mmd_mfa <- rbind(bacteria_dt_subset2, fungi_dt_subset2) #FABIAobject.mmd.trained@X#
mmd_mfa <- as.data.frame(mmd_mfa)
#CLL_mfa$feature <- rownames(as.data.frame(CLL_mfa))
mmd_mfa_tr <- data.frame(t(mmd_mfa))
MFAobject.mmd.trained <- FactoMineR::MFA(mmd_mfa_tr, group = c(180, 18), type = c("s", "s"),
                                         name.group = c("Bacteria", "Fungi"),
                                         graph = FALSE)

# Factor scores
factor_mfa = data.frame(MFAobject.mmd.trained$ind$coord)
colnames(factor_mfa) <- c("score_MFA1", "score_MFA2", "score_MFA3", "score_MFA4", "score_MFA5")
factor_mfa$sample = rownames(factor_mfa)
factor_mfa2 = factor_mfa[,-c(3:5)]

# Weights/Loading 
loading_mfa = data.frame(MFAobject.mmd.trained$quanti.var$coord)
colnames(loading_mfa) <- c("loading_MFA1", "loading_MFA2", "loading_MFA3", "loading_MFA4", "loading_MFA5")
loading_mfa$feature = rownames(loading_mfa)
loading_mfa2 = loading_mfa[, -c(3:5)]

# GFA
fabia.scaledData.mmd <- FABIAobject.mmd.trained@X
gfa_mmd.data <- as.data.frame(fabia.scaledData.mmd)
merged_GFA_mmd.data = list(t(gfa_mmd.data))
model_option <- getDefaultOpts()
model_option$iter.max <- 1000
model_option$iter.burnin <- 10
GFAobject.mmd.trained <- gfa(t(merged_GFA_mmd.data), K= 2, opts=model_option)

# Factor scores
factor_gfa = as.data.frame(GFAobject.mmd.trained$X)
colnames(factor_gfa) <- c("score_GFA1", "score_GFA2")
factor_gfa$sample = rownames(factor_gfa)

# Weights/Loading 
loading_gfa = as.data.frame(GFAobject.mmd.trained$W)
colnames(loading_gfa) <- c("loading_GFA1", "loading_GFA2")
loading_gfa$feature = rownames(loading_gfa)

# Merging Factor scores data for FABIA, MOFA, MFA, GFA
merged_fabia_mofa <- merge(factor_fabia2, factor_mofa2, by = 'sample', all = TRUE) # Merge factor_fabia2 and factor_mofa2
merged_fabia_mofa_mfa <- merge(merged_fabia_mofa, factor_mfa2, by = 'sample', all = TRUE) # Merge merged_fabia_mofa with factor_mfa2
data_scores <- merge(merged_fabia_mofa_mfa, factor_gfa, by = 'sample', all = TRUE) # Merge merged_fabia_mofa_mfa with factor_gfa

# Merging Weights/loading data for FABIA, MOFA, MFA, GFA
merged_fabia_mofa_loading <- merge(loading_fabia2, loading_mofa2, by = 'feature', all = TRUE) # Merge loading_fabia2 and loading_mofa2
merged_fabia_mofa_mfa_loading <- merge(merged_fabia_mofa_loading, loading_mfa2, by = 'feature', all = TRUE) # Merge merged_fabia_mofa_loading with loading_mfa2
data_weights <- merge(merged_fabia_mofa_mfa_loading, loading_gfa, by = 'feature', all = TRUE) # Merge merged_fabia_mofa_mfa with loading_gfa

# Visuaization of scores using correlation matrix plot
library(GGally); library(psych)

data_scores_cor_F1 <- data_scores
ggpairs(data_scores_cor_F1[c("score_FABIA2", "score_MOFA1", "score_MFA1", "score_GFA1")])

# Define new variable names
new_var_names <- c("Score (FABIA 2)", "Score (MOFA 1)", "Score (MFA 1)", "Score (GFA 1)")
colnames(data_scores_cor_F1)[c("score_FABIA1", "score_MOFA2", "score_MFA2", "score_GFA2")] <- new_var_names
pairs.panels(data_scores_cor_F1[c("score_FABIA1", "score_MOFA2", "score_MFA2", "score_GFA2")], hist.col="lightgoldenrod", breaks = 16, pch = 8, cex = 0.5, ellipses=TRUE, lm=TRUE)

# Visualize scores using scatter plot
# Plot
data_scores_cor_F1$Category <- as.factor(data_scores_cor_F1$Category)
ggplot(data_scores_cor_F1, aes(x = score_MFA1, y = -1*score_MFA2)) +
  geom_point(aes(fill = Category), size = 4, shape = 21) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  ggtitle("GFA: F 1") +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = c("Healthy, antibiotics" = "yellow","Healthy, no antibiotics" = "steelblue1", "Non septic ICU" = "grey", "Sepsis" = "green")) +  # Manually assign colors to Condition
  labs(x = "Score 1", y = "Score 2")

# Visualize weights/loading using scatter plot

# Visuaization of weights using correlation matrix plot
library(GGally); library(psych)

data_weights_cor_mmd <- data_weights
ggpairs(data_weights_cor[c("loading_FABIA2", "loading_MOFA1", "loading_MFA1", "loading_GFA1")])


# Filter the data for only View=drugs
data_weights_cor_bacteria <- data_weights_cor_mmd %>%
  filter(view == "Bacteria")
# Filter the data for only View=methylation
data_weights_cor_fungi <- data_weights_cor_mmd %>%
  filter(view == "Fungi")

# Correlation plots
# Define new variable names
new_var_names <- c("Loading (FABIA)", "Loading (MOFA)", "Loading (MFA)", "Loading (GFA)")
colnames(data_weights_cor_meth)[c(2,4,6,7)] <- new_var_names
pairs.panels(data_weights_cor_bacteria[c("loading_FABIA1", "loading_MOFA2", "loading_MFA2", "loading_GFA2")], hist.col="deepskyblue4", breaks = 16, pch = 8, cex = 0.5, ellipses=TRUE, lm=TRUE)

# Plot using a single palette
plot1 = ggplot(data_weights_cor_drugs, aes(x = `Loading (MFA)`, y = `Loading (GFA)`)) + 
  geom_point(color = "dodgerblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Loading (MFA)", y = "Loading (GFA)") +
  ggtitle("Drugs") +
  theme(
    axis.ticks.x = element_blank(),
    plot.title = element_text(hjust = 1))+
  theme_minimal()

# Second plot
plot2 = ggplot(data_weights_cor_meth, aes(x = `Loading (MFA)`, y = `Loading (GFA)`)) + 
  geom_point(color = "deepskyblue4") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Loading (MFA)", y = "Loading (MFA)") +
  ggtitle("Methylation") +
  theme_minimal()

# Arrange plots side by side
grid.arrange(plot1, plot2, ncol = 2)

# Functions 
# Threshold calculator
calculate_threshold_and_index <- function(data) {
  for (i in seq_along(data)) {
    if (is.numeric(data[[i]])) {
      threshold <- varphi.data(data[[i]])
      
      # Create corresponding index column name
      index_column <- paste0(names(data)[i], "_index")
      
      # Create index column based on the threshold
      data[[index_column]] <- ifelse(data[[i]] > threshold, 1, 0)
    }
  }
  return(data)
}

# Jaccard Similarity Index
varphi.data <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    varphi <- quantile(numeric_values, probs = 0.60)
    return(varphi)
  } else {
    return(NA)  # Return NA if not all values are numeric
  }
}

# Function to re-scale a vector to (-1,1) scale
rescale_minusone_to_one <- function(vector){
  scaled_vector <- vector / max(abs(vector), 1e-6) # Using 1e-6 to avoid division by zero
  return(scaled_vector)
}

data_scores_cor_f2_subset.i <- data_scores_cor_F1[,c("score_FABIA1", "score_MOFA2", "score_MFA2", "score_GFA2")]
# subset data to only scores
data_scores_cor_f2_subset.ii <- df_normalized[,c("score_FABIA1", "score_MOFA2", "score_MFA2", "score_GFA2")]

# Obtain thresholds
#data_scores_cor_f1_subset.i <- data_scores_cor_f1_subset %>%
#  mutate_if(is.numeric, ~ rescale_minusone_to_one(.))

data_scores_cor_f2_subset.a <- data_scores_cor_f2_subset.ii %>%
  mutate_if(is.numeric, ~ abs(.))

varphi.data(abs(data_scores_cor_f2_subset.a$score_FABIA1))
varphi.data(abs(data_scores_cor_f2_subset.a$score_MOFA2))
varphi.data(abs(data_scores_cor_f2_subset.a$score_MFA2))
varphi.data(abs(data_scores_cor_f2_subset.a$score_GFA2))

# calculate_threshold_and_index
data_scores_cor_f2_subset.b <- calculate_threshold_and_index(data_scores_cor_f2_subset.a)
data_scores_cor_f2_index.b <- grep("_index$", names(data_scores_cor_f2_subset.b), value = TRUE)

# Subset the dataframe based on selected columns
data_scores_cor_f2_subset.c <- data_scores_cor_f2_subset.b[, c(data_scores_cor_f2_index.b)]
jaccard_scores = jaccard.index(data_scores_cor_f2_subset.c)

# 2. Loading/Weight
data_weights_cor_bacteria_F1 <- data_weights_cor_bacteria[,c("loading_FABIA1", "loading_MOFA1", "loading_MFA1", "loading_GFA1")]
data_weights_cor_fungi_F1 <- data_weights_cor_fungi[,c("loading_FABIA1", "loading_MOFA1", "loading_MFA1", "loading_GFA1")]
data_weights_cor_bacteria_F2 <- data_weights_cor_bacteria[,c("loading_FABIA1", "loading_MOFA2", "loading_MFA2", "loading_GFA2")]
data_weights_cor_fungi_F2 <- data_weights_cor_fungi[,c("loading_FABIA1", "loading_MOFA2", "loading_MFA2", "loading_GFA2")]

# Jaccard Similarity Index
varphi.data <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    varphi <- quantile(numeric_values, probs = 0.60)
    return(varphi)
  } else {
    return(NA)  # Return NA if not all values are numeric
  }
}
data_weights_cor_bacteria_F1$loading_MFA1 = as.numeric(data_weights_cor_bacteria_F1$loading_MFA1)
# Obtain thresholds
data_weights_cor_bacteria_subset.i <- data_weights_cor_fungi_F1 %>%
  mutate_if(is.numeric, ~ rescale_minusone_to_one(.))

data_weights_cor_bacteria_subset.a <- data_weights_cor_bacteria_subset.i %>%
  mutate_if(is.numeric, ~ abs(.))

varphi.data(abs(data_weights_cor_bacteria_subset.a$loading_FABIA2))
varphi.data(abs(data_weights_cor_bacteria_subset.a$loading_MOFA1))
varphi.data(abs(data_weights_cor_bacteria_subset.a$loading_MFA1))
varphi.data(abs(data_weights_cor_bacteria_subset.a$loading_GFA1))

# calculate_threshold_and_index
data_weights_cor_drugs_subset.b <- calculate_threshold_and_index(data_weights_cor_drugs_subset.a)
data_weights_cor_drugs_index.b <- grep("_index$", names(data_weights_cor_drugs_subset.b), value = TRUE)

# Subset the dataframe based on selected columns
data_weights_cor_drugs_subset.c <- data_weights_cor_drugs_subset.b[, c(data_weights_cor_drugs_index.b)]
jaccard_scores_drugs = jaccard.index(data_weights_cor_drugs_subset.c)

# (ii) Methylation
varphi.data <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    varphi <- quantile(numeric_values, probs = 0.95)
    return(varphi)
  } else {
    return(NA)  # Return NA if not all values are numeric
  }
}
# subset data to only loading
data_weights_cor_meth_subset <- data_weights_cor_meth[,c(2,4,6,7)]

# Obtain thresholds
data_weights_cor_meth_subset.i <- data_weights_cor_meth_subset %>%
  mutate_if(is.numeric, ~ rescale_minusone_to_one(.))

data_weights_cor_meth_subset.a <- data_weights_cor_meth_subset.i %>%
  mutate_if(is.numeric, ~ abs(.))

varphi.data(abs(data_weights_cor_meth_subset.a$`Loading (FABIA)`)); varphi.data(abs(data_weights_cor_meth_subset.a$`Loading (MOFA)`)); varphi.data(abs(data_weights_cor_meth_subset.a$`Loading (MFA)`)); varphi.data(abs(data_weights_cor_meth_subset.a$`Loading (GFA)`))

# calculate_threshold_and_index
data_weights_cor_meth_subset.b <- calculate_threshold_and_index(data_weights_cor_meth_subset.a)
data_weights_cor_meth_index.b <- grep("_index$", names(data_weights_cor_meth_subset.b), value = TRUE)

# Subset the dataframe based on selected columns
data_weights_cor_meth_subset.c <- data_weights_cor_meth_subset.b[, c(data_weights_cor_meth_index.b)]
jaccard_scores_meth = jaccard.index(data_weights_cor_meth_subset.c)

# ------------------------------------------------ FUNCTIONS -----------------------------------------------------------------

# Jaccard index
jaccard.similarity.index <- function(data){
  # Extract index columns
  index_columns <- grep("_index$", names(data), value = TRUE)
  # Initialize an empty dataframe with a single row of NA values
  jaccard_df <- data.frame(
    Column1 = NA_character_,
    Column2 = NA_character_,
    similar_count = NA_real_,
    dissimilar_count = NA_real_,
    sum_column1 = NA_real_,
    sum_column2 = NA_real_,
    jaccard_index = NA_real_,
    stringsAsFactors = FALSE
  )
  jaccard_df <- jaccard_df[0, ]  # Remove the initial row
  # Calculate Jaccard similarity for each pair of index columns
  for (i in 1:(length(index_columns) - 1)) {
    for (j in (i + 1):length(index_columns)) {
      column1 <- index_columns[i]
      column2 <- index_columns[j]
      if (length(data[[column1]]) != length(data[[column2]])) {
        stop("Columns must have the same length.")
      }
      # Calculate Jaccard similarity index
      similar_count <- sum(data[[column1]] == 1 & data[[column2]] == 1)
      dissimilar_count <- sum(data[[column1]] != data[[column2]])
      sum_column1 <- sum(data[[column1]] == 1)
      sum_column2 <- sum(data[[column2]] == 1)
      jaccard_index <- similar_count / (sum_column1 + sum_column2 - similar_count)
      ji = ifelse(jaccard_index == -1, 1, jaccard_index)
      # Append the results to the dataframe
      ji <- round(ji, 2)
      jaccard_df <- rbind(jaccard_df, c(column1, column2, ji, similar_count, dissimilar_count, sum_column1, sum_column2))
    }
  }
  # Rename the columns in the final result
  colnames(jaccard_df) <- c("Column1", "Column2", "ji", "similar_count", "dissimilar_count", "sum_column1", "sum_column2")
  return(jaccard_df)
}
# ---------------------------------------------- CREATING SCORES DATAFRAME ---------------------------------------------------
# merge with metadata
df_sampleInfo_Hip_D17b$samples = rownames(df_sampleInfo_Hip_D17b)#$Internal.ID
data_score_methods.17 = merge(data_score_df.17, df_sampleInfo_Hip_D17b, by = 'samples', all = TRUE)

# Jaccard index for scores
ji_data_scores_cor <- data_scores_cor[c(1,2,11,12,13)]
new_var_names <- c('sample', 'FABIA_score', 'MOFA_score', 'MFA_score', 'GFA_score')
colnames(ji_data_scores_cor) <- new_var_names
data_scores_cor_index <- calculate_threshold_and_index(ji_data_scores_cor)
score_index_columns.cll <- grep("_index$", names(data_scores_cor_index), value = TRUE)

# Jaccard Index
jaccard_score_data.cll = jaccard.similarity.index(data_scores_cor_index)

# Venn
# install.packages("ggVennDiagram")
library(ggVennDiagram)
x = data_scores_cor_index
# List of items
x <- list(x$FABIA_score_index, x$MOFA_score_index)

# 2D Venn diagram
ggVennDiagram(x) 
