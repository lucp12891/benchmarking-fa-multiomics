setwd("/user/leuven/364/vsc36498")
set.seed(60667) 

# Load the required libraries
library(tidyr)
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

##################### 
# START OF FUNCTIONS  #
#####################
#  --------------------------------- MultiFactor: Simulation function ----------------------------
multiple_factor <- function(n_features_one, n_features_two, n_samples, sigmas, iterations, n_factors){
  set.seed(694)
  
  # Initialize lists to store omic data for all iterations and sigmas
  all_omic_data <- list()
  
  for (iter in 1:iterations) {  # Iterate through the specified number of iterations
    omic.one <- list()
    omic.two <- list()
    
    for (k in 1:length(sigmas)) {  # Iterate through the list of sigmas
      n_s <- n_samples
      valid_data <- FALSE  # Flag to check if data passes FABIA condition
      
      #num_factor = num_factors
      omic.one <- list()
      omic.two <- list()
      
      # Generate first omic data
      s_sig_s1 <- ceiling(n_s / 5.3)   # Start index for first range
      s_sig_e1 <- ceiling(n_s / 2.8)   # End index for first range
      s_sig_s2 <- 50                   # Start index for second range
      s_sig_e2 <- 70                   # End index for second range
      
      alpha <- rnorm(n_s, 0, 0.05)
      
      # Assign indices for the first range
      assigned_indices_samples1 <- sample(s_sig_s1:s_sig_e1, length(s_sig_s1:s_sig_e1))
      
      # Assign indices for the second range
      assigned_indices_samples2 <- sample(s_sig_s2:s_sig_e2, length(s_sig_s2:s_sig_e2))
      
      # Combine the two sets of indices
      assigned_indices_samples <- list(assigned_indices_samples1, assigned_indices_samples2)
      
      max_factors <- length(assigned_indices_samples)# Same as n_factors
      
      
      # Generate random alpha values based on the max_factors
      list_alphas <- list()
      for (i in 1:max_factors) {
        list_alphas[[paste0("alpha", i)]] <- rnorm(n_samples, 0, 0.05)
      }
      
      # Assign corresponding values to alpha variables based on assigned_indices_samples
      for (i in seq_along(assigned_indices_samples)) {
        indices <- assigned_indices_samples[[i]]
        list_alphas[[i]][indices] <- rnorm(length(indices), (3 + 0.5*i), 0.05)  # Adjust values dynamically
      }
      
      list_gammas <- list_alphas
      
      # Add the separated list of the indices selected for omics of multiple features
      # Initial vector
      vector <- c(1:max_factors)
      
      # Specify the number of elements to sample for the shared vector
      num_shared <- ceiling(runif(1, min = 1, max = length(vector) - 1))
      
      # Select the 'shared' elements
      if (num_shared > 0) {
        shared <- 1
        # Remove the selected shared elements from the vector
        remain_vector <- vector[!vector %in% shared]
      }
      
      # Shuffle the remaining vector
      if(length(remain_vector == 1)){
        shuffled_remain = remain_vector
      } else{
        shuffled_remain <- sample(remain_vector)
      }
      # Split the shuffled elements equally or nearly equally between omic_one_unique and omic_two_unique
      num_elements <- length(shuffled_remain)
      split_point <- sample(1:num_elements, 1) # Randomly select a split point
      
      omic_one_unique <- shuffled_remain[0:split_point] # First part goes to omic_one_unique
      omic_two_unique <- shuffled_remain[!shuffled_remain %in% omic_one_unique] # Second part goes to omic_two_unique
      
      # Print the results
      list(
        shared = shared,
        omic_one_unique = omic_one_unique,
        omic_two_unique = omic_two_unique
      )
      
      # End of indices assignment
      
      # Assignment of factors
      list_omic_one_factors = c(shared, omic_one_unique)
      list_omic_two_factors = c(shared, omic_two_unique)
      #all_indices = c(shared, omic_one_unique, omic_two_unique)
      
      factor_xtics <- list(
        shared = shared,
        omic_one_unique = omic_one_unique,
        omic_two_unique = omic_two_unique
      )
      
      # Generate first omic data
      
      # Define the first range of feature indices
      f_sig_s1 <- 3650      # Start index for first range
      f_sig_e1 <- 4000    # End index for first range
      
      # Define the second range of feature indices
      f_sig_s2 <- 1600   # Start index for second range
      f_sig_e2 <- 1950   # End index for second range
      
      # Assign indices for the first range
      assigned_indices_features1 <- sample(f_sig_s1:f_sig_e1, length(f_sig_s1:f_sig_e1))
      
      # Assign indices for the second range
      assigned_indices_features2 <- sample(f_sig_s2:f_sig_e2, length(f_sig_s2:f_sig_e2))
      
      # Combine the two sets of indices
      assigned_indices_features <- list(assigned_indices_features1, assigned_indices_features2)
      
      # Generate random beta values based on the max_factors
      list_betas <- list()
      
      # Shuffle the list
      (shuffled_assigned_indices_features <- assigned_indices_features)# sample(assigned_indices_features)) # Do not reshuffle first
      
      # Create only one vector in list_betas using the random index
      for (i in seq_along(assigned_indices_features)) {
        list_betas[[i]] <- rnorm(n_features_one, 0, 0.05)
      }
      
      # Loop through the ordered indices in list_omic_one_factors
      for (i in seq_along(assigned_indices_features)) {
        indices_ns <- assigned_indices_features[[i]]  # Get the corresponding indices from assigned_indices_features
        if (length(indices_ns) > 0) {
          # Initialize list_betas[[i]] if it's not already
          if (is.null(list_betas[[i]])) {
            list_betas[[i]] <- numeric(max(indices_ns))
          }
          # Assign values to the appropriate indices in list_betas
          list_betas[[i]][indices_ns] <- rnorm(length(indices_ns), mean = (4.0 + 0.5 * i), sd = 0.05)
        }
      }
      
      # Renaming list_betas to desired factor number
      for (i in seq_along(list_omic_one_factors)) {
        names(list_betas)[i] <- paste0("beta", list_omic_one_factors[[i]])#list_omic_one_factors[[i]]
      }
      
      pattern_alpha <-"^alpha\\d+"
      matches_alpha <- grepl(pattern_alpha, names(list_alphas), perl = TRUE)
      n_alpha = length(matches_alpha[matches_alpha == TRUE])
      
      pattern_beta <- "^beta\\d+"
      matches_beta <- grepl(pattern_beta, names(list_betas), perl = TRUE)
      n_beta = length(matches_beta[matches_beta == TRUE])
      
      # Initialize the data list to store the results
      data_list_i <- list()
      
      # Sort the list by names in ascending order
      list_betas <- list_betas[order(names(list_betas))]
      list_alphas <- list_alphas[order(names(list_alphas))]
      
      # Extract the last numeric values from the names of the vectors
      alphas_names <- as.numeric(gsub("[^0-9]", "", names(list_alphas)))
      betas_names <- as.numeric(gsub("[^0-9]", "", names(list_betas)))
      
      # Find common numbers between the last values of the vector names
      common_names <- intersect(alphas_names, betas_names)
      
      # Filter the vectors in each list based on the common names
      list_alphas <- list_alphas[paste0("alpha", common_names)]
      list_betas <- list_betas[paste0("beta", common_names)]
      
      # Loop through each alpha and beta combination
      for (i in 1:min(length(list_alphas), length(list_betas))) {
        data_i <- list_alphas[[i]] %*% t(list_betas[[i]][1:n_features_one])
        data_list_i[[paste0("data.", i)]] <- data_i
      }
      
      # Combine the results into a single data variable
      data.1 <- Reduce(`+`, data_list_i)
      
      eps1 <- rnorm(n_samples * n_features_one, 0, sigmas) # noise
      omic1_data <- matrix(data.1, n_samples, n_features_one) + matrix(eps1, n_samples, n_features_one) # signal + noise
      colnames(omic1_data) <- paste0('omic1_feature_', seq_len(n_features_one))
      rownames(omic1_data) <- paste0('sample_', seq_len(n_samples))
      
      
      names(list_gammas) <- gsub("alpha", "gamma", names(list_gammas))
      
      # Extract the numeric parts from the names of list_gammas
      numeric_part <- as.numeric(gsub("\\D", "", names(list_gammas)))
      
      # Retain elements where the numeric part of the name matches values in list_omic_two_factors
      list_gammas <- list_gammas[numeric_part %in% list_omic_two_factors]
      
      # Empty list
      list_deltas <- list()
      
      # Define the range for the feature indices
      f_sig2_s1 <- 1   # Start index
      f_sig2_e1 <- 300   # End index
      
      # Sample exactly 300 indices from the defined range
      assigned_indices_features2 <- sample(f_sig2_s1:f_sig2_e1, 300)
      assigned_indices_features_omic.two <- list(assigned_indices_features2)
      
      # Create only one vector in list_betas using the random index
      for (i in seq_along(assigned_indices_features_omic.two)) {
        list_deltas[[i]] <- rnorm(n_features_two, 0, 0.05)
      }
      
      
      for (i in seq_along(assigned_indices_features_omic.two)) {
        indices_ns <- assigned_indices_features_omic.two[[i]]  # Get the corresponding indices from assigned_indices_features
        
        if (length(indices_ns) > 0) {
          # Ensure list_deltas has a slot for the current index
          if (length(list_deltas) < i) {
            list_deltas[[i]] <- numeric(0)  # Initialize an empty numeric vector if it doesn't exist
          }
          
          # Ensure list_deltas[[i]] has enough space for indices_ns
          if (length(list_deltas[[i]]) < max(indices_ns)) {
            length(list_deltas[[i]]) <- max(indices_ns)
          }
          
          # Assign values to the appropriate indices in list_deltas
          list_deltas[[i]][indices_ns] <- rnorm(length(indices_ns), mean = (5.0 + 0.0 * i), sd = 0.05)
        }
      }
      
      
      # Renaming list_deltas to desired factor number
      for (i in seq_along(list_omic_two_factors)) {
        names(list_deltas)[i] <- paste0("delta", list_omic_two_factors[[i]])#list_omic_one_factors[[i]]
      }
      
      pattern_gamma <-"^gamma\\d+"
      matches_gamma <- grepl(pattern_gamma, names(list_gammas), perl = TRUE)
      n_gamma = length(matches_alpha[matches_gamma == TRUE])
      
      pattern_delta <- "^delta\\d+"
      matches_delta <- grepl(pattern_delta, names(list_deltas), perl = TRUE)
      n_delta = length(matches_delta[matches_delta == TRUE])
      
      # Initialize the data list to store the results
      data_list_j <- list()
      
      # Sort the list by names in ascending order
      list_gammas <- list_gammas[order(names(list_gammas))]
      list_deltas <- list_deltas[order(names(list_deltas))]
      
      # Extract the last numeric values from the names of the vectors
      gammas_names <- as.numeric(gsub("[^0-9]", "", names(list_gammas)))
      deltas_names <- as.numeric(gsub("[^0-9]", "", names(list_deltas)))
      
      # Find common numbers between the last values of the vector names
      common_names <- intersect(gammas_names, deltas_names)
      
      # Filter the vectors in each list based on the common names
      list_gammas <- list_gammas[paste0("gamma", common_names)]
      list_deltas <- list_deltas[paste0("delta", common_names)]
      
      # Loop through each alpha and beta combination
      for (j in 1:min(length(list_gammas), length(list_deltas))) {
        data_j <- list_gammas[[j]] %*% t(list_deltas[[j]][1:n_features_two])
        data_list_j[[paste0("data.", j)]] <- data_j
      }
      
      # Combine the results into a single data variable
      data.2 <- Reduce(`+`, data_list_j)
      
      eps2 <- rnorm(n_samples * n_features_two, 0, sigmas)
      omic2_data <- matrix(data.2, n_samples, n_features_two) + matrix(eps2, n_samples, n_features_two) # signal + noise
      colnames(omic2_data) <- paste0('omic2_feature_', seq_len(n_features_two))
      rownames(omic2_data) <- paste0('sample_', seq_len(n_samples))
      
      # Concatenate datasets
      concatenated_data <- cbind(omic2_data, omic1_data)
      
      # Save the validated data
      omic.one[[k]] <- omic1_data
      omic.two[[k]] <- omic2_data
    }
    
    # Combine datasets for this iteration
    simulated_datasets <- list(object.one = omic.one, object.two = omic.two)
    concatenated_datasets <- list()
    for (i in seq_along(simulated_datasets$object.one)) {
      concatenated_datasets[[i]] <- cbind(simulated_datasets$object.one[[i]], simulated_datasets$object.two[[i]])
    }
    
    # Save data for the current iteration
    all_omic_data[[paste0("iteration_", iter)]] <- list(
      concatenated_datasets = concatenated_datasets,
      indices_features.1 = assigned_indices_features,
      indices_samples = assigned_indices_samples,
      indices_features.2 = assigned_indices_features_omic.two,
      sample_sig_start1 = s_sig_s1,
      sample_sig_end1 = s_sig_e1,
      sample_sig_start2 = s_sig_s2,
      sample_sig_end2 = s_sig_e2,
      feature_sig_start = f_sig_s1,
      feature_sig_end = f_sig_e1,
      feature_sig_start = f_sig_s2,
      feature_sig_end = f_sig_e2,
      feature_sig2_start = f_sig2_s1,
      feature_sig2_end = f_sig2_e1,
      list_alphas = list_alphas,
      list_gammas = list_gammas,
      list_betas = list_betas,
      list_deltas = list_deltas
    )
  }
  
  return(all_omic_data)
}

#simulatedData <- multiple_factor(n_features_one = 4000, n_features_two = 3000, n_samples = 100, sigmas = 7, n_factors = 2, iterations = 1)
#   
# dataset <- simulatedData$iteration_1$concatenated_datasets[[1]]
# set.seed(657)
# image(c(1:dim(dataset)[1]), c(1:dim(dataset)[2]), dataset, ylab = "Features", xlab = "Samples", font.lab = 2)
# abline(h = 3000, col = "brown", lty = "dashed")  # lty = "dashed" sets the line type
# # Add text annotations for the matrices
# # Adjust x and y coordinates based on the location of the dark matrices
# text(x = 28, y = 2550, labels = "shared factor", col = "black", cex = 0.8, font = 2)
# text(x = 60.5, y = 5540, labels = "unique factor", col = "black", cex = 0.80, font = 2)
# persp(c(1:dim(t(dataset))[2]), c(1:dim(t(dataset))[1]), dataset, theta = 325, phi = 15, col = "yellow", xlab = "", ylab = "Samples", zlab = " ")

# --------------------------------------- FABIA Function -----------------------------------------

## ----------------------------------- 1. FABIA Main Function ---------------------------------------

func_fabia <- function(data, BC_num) {
  fabia_seed <- 60667#sample(1e6, 1)
  set.seed(fabia_seed)
  
  valid_data <- FALSE  # Initialize validation flag
  
  n_features_one <- 4000 #global_env$n_features_one
  n_features_two <- 3000 #global_env$n_features_two
  
  fab_data <- data
  
  second_omic <- t(fab_data[, 1:n_features_two])
  first_omic <- t(fab_data[, (n_features_two + 1):(n_features_one + n_features_two)])
  r_data = rbind(first_omic, second_omic)
  
  while (!valid_data) {
    #set.seed(123)
    
    # Run FABIA
    fabia_object <- fabia(as.matrix(r_data), p = BC_num, alpha = 0.01, cyc = 1000, spl = 0.5, spz = 0.5, 
                          random = 1.0, center = 2, norm = 2, lap = 1.0, nL = 1)
    
    # Extract feature weights and factor scores
    # Iterate over the number of biclusters
    fabia_loading <- fab_loading(fabia_object, BC_num)
    fabia_score <- fab_score(fabia_object, BC_num)
    
    # Function to validate FABIA results for a specific column
    validate_column <- function(column) {
      !any(is.na(column)) && mean(column, na.rm = TRUE) != 0
    }
    
    # Apply validation to all required columns
    valid_data <- validate_column(fabia_loading$loading_FABIA1) &&
      validate_column(fabia_score$score_FABIA1) &&
      validate_column(fabia_loading$loading_FABIA2) &&
      validate_column(fabia_score$score_FABIA2)
    
    # Check and handle validation results
    if (valid_data) {
      message("Data passes validation.")
    } else {
      message("FABIA validation failed, re-run FABIA...")
    }
  }
  
  # FABIA normalized data X
  X <- fabia_object@X
  
  # Loading separation per factor
  result_load <- fabia_loading
  
  result_load_all <- result_load %>%
    mutate(
      dataview = ifelse(grepl("omic1_", feature), "omic.one",
                        ifelse(grepl("omic2_", feature), "omic.two", NA)),
      ID = ifelse(grepl("omic1_", feature),
                  as.numeric(sub(".*omic1_feature_(\\d+)", "\\1", feature)),
                  ifelse(grepl("omic2_", feature),
                         as.numeric(sub(".*omic2_feature_(\\d+)", "\\1", feature)),
                         NA))
    ) %>%
    arrange(dataview, ID)
  
  load_omic_one <- result_load %>%
    mutate(
      dataview = ifelse(grepl("omic1_", feature), "omic.one",
                        ifelse(grepl("omic2_", feature), "omic.two", NA)),
      ID = ifelse(grepl("omic1_", feature),
                  as.numeric(sub(".*omic1_feature_(\\d+)", "\\1", feature)),
                  ifelse(grepl("omic2_", feature),
                         as.numeric(sub(".*omic2_feature_(\\d+)", "\\1", feature)),
                         NA))
    ) %>%
    arrange(dataview, ID) %>%
    filter(dataview == 'omic.one')
  
  load_omic_two <- result_load %>%
    mutate(
      dataview = ifelse(grepl("omic1_", feature), "omic.one",
                        ifelse(grepl("omic2_", feature), "omic.two", NA)),
      ID = ifelse(grepl("omic1_", feature),
                  as.numeric(sub(".*omic1_feature_(\\d+)", "\\1", feature)),
                  ifelse(grepl("omic2_", feature),
                         as.numeric(sub(".*omic2_feature_(\\d+)", "\\1", feature)),
                         NA))
    ) %>%
    arrange(dataview, ID) %>%
    filter(dataview == 'omic.two')
  
  weights_composite <- list(
    omic.one_weights = load_omic_one,
    omic.two_weights = load_omic_two,
    all_weights = result_load_all
  )
  
  # Factor scores separation
  result_score <- fabia_score
  
  result_score <- result_score %>%
    mutate(
      ID = ifelse(grepl("sample_", sample),
                  as.numeric(sub(".*sample_(\\d+)", "\\1", sample)), NA)
    ) %>%
    arrange(ID) #%>%
  
  scores_composite <- list(
    
    scores = result_score
  )
  
  fabia_result <- list(weights = weights_composite, scores = scores_composite)
  
  return(fabia_result)
}

## ----------------------------------- 1.1 FABIA Loading Function ---------------------------------------
fab_loading <- function(fabia_object, BC_num) {
  # Check if BC_num is greater than 1
  if (BC_num > 1) {
    # Initialize an empty data frame for results
    loadings_df <- data.frame()
    
    # Extract loadings for each BC_num and combine into a single data frame
    for (i in 1:BC_num) {
      loading_column <- fabia_object@L[, i]  # Extract the loading for the current factor
      column_name <- paste0("loading_FABIA", i)  # Name for the column
      
      # Add the loading column to the data frame
      if (ncol(loadings_df) == 0) {
        loadings_df <- data.frame(feature = rownames(fabia_object@L))  # Initialize with features names
      }
      loadings_df[[column_name]] <- loading_column
    }
    BC_df_loading = loadings_df
    return(BC_df_loading)
  } else {
    # Handle single BC_num case
    loading_FABIA <- fabia_object@L[, BC_num]
    BC_df_loading <- as.data.frame(loading_FABIA)
    BC_df_loading$feature <- rownames(BC_df_loading)
    return(BC_df_loading)
  }
}

## ------------------------------------ 1.2 FABIA Score Function ----------------------------------------

fab_score <- function(fabia_object, BC_num) {
  # Check if BC_num is greater than 1
  if (BC_num > 1) {
    # Initialize an empty data frame for results
    scores_df <- data.frame()
    
    # Extract scores for each BC_num and combine into a single data frame
    for (i in 1:BC_num) {
      score_column <- fabia_object@Z[i, ]  # Extract the score for the current factor
      column_name <- paste0("score_FABIA", i)  # Name for the column
      
      # Add the score column to the data frame
      if (ncol(scores_df) == 0) {
        scores_df <- data.frame(sample =colnames(fabia_object@Z))  # Initialize with sample names
      }
      scores_df[[column_name]] <- score_column
    }
    BC_df_score = scores_df
    return(BC_df_score)
  } else {
    # Handle single BC_num case
    score_FABIA <- fabia_object@Z[BC_num, ]
    BC_df_score <- as.data.frame(score_FABIA)
    BC_df_score$sample <- colnames(BC_df_score)
    return(BC_df_score)
  }
}

# --------------------------------------- MOFA+ Function -----------------------------------------

## ----------------------------------- 2. MOFA+ Main Function ---------------------------------------
func_mofa <- function(data, num) {
  mofa_seed <- 60667#sample(1e6, 1)
  set.seed(mofa_seed)
  
  # Load necessary package
  library(basilisk)
  
  # Prepare data
  mofa_data <- data
  
  # Access variables from the shared environment
  n_features_one <- 4000#global_env$n_features_one
  n_features_two <- 3000#global_env$n_features_two
  
  # Subset respective datasets
  second_omic <- mofa_data[, 1:n_features_two]
  first_omic <- mofa_data[, (n_features_two + 1):(n_features_one + n_features_two)]
  
  # Create a data list 'ready' for MOFA
  mofa_data_sim <- list(
    first_omic = as.matrix(t(first_omic)),
    second_omic = as.matrix(t(second_omic))
  )
  
  # Step 1: Create a MOFA object with feature names
  feature_names_first_omic <- colnames(first_omic)
  feature_names_second_omic <- colnames(second_omic)
  
  mofa_data_sim <- create_mofa(
    mofa_data_sim,
    use_basilisk = TRUE,
    feature_names = list(feature_names_first_omic, feature_names_second_omic)
  )
  
  # Step 2: Defining options: Data, Model, Training
  # (2a) Define data options
  data_opts_sim <- get_default_data_options(mofa_data_sim)
  data_opts_sim$scale_views <- FALSE
  data_opts_sim$scale_groups <- TRUE
  data_opts_sim$center_groups <- TRUE
  
  # (2b) Define model options
  model_opts_sim <- get_default_model_options(mofa_data_sim)
  model_opts_sim$num_factors <- num
  
  # (2c) Define training options
  train_opts_sim <- get_default_training_options(mofa_data_sim)
  train_opts_sim$maxiter <- 1000
  train_opts_sim$convergence_mode <- "fast"
  train_opts_sim$seed <- 123
  
  # Step 3: Build and train the MOFA object
  MOFAobject_sim <- prepare_mofa(
    object = mofa_data_sim,
    data_options = data_opts_sim,
    model_options = model_opts_sim,
    training_options = train_opts_sim
  )
  
  outfile_sim <- file.path(tempdir(), "model_sim.hdf5")
  
  # Suppress specific warnings
  suppressWarnings({
    MOFAobject_sim.trained <- run_mofa(
      MOFAobject_sim,
      outfile_sim,
      use_basilisk = TRUE
    )
    
    model_object <- MOFAobject_sim.trained
    
    # Lists to store loading and scores for all factors
    mofa_loading <- mofa_loading(model_object, factor_num = num)
    mofa_score <- mofa_score(model_object, factor_num = num)
    
    # Loading separation per factor
    result_load <- mofa_loading
    
    result_load_all <- result_load %>%
      mutate(
        dataview = ifelse(grepl("omic1_", feature), "omic.one",
                          ifelse(grepl("omic2_", feature), "omic.two", NA)),
        ID = ifelse(grepl("omic1_", feature),
                    as.numeric(sub(".*omic1_feature_(\\d+)", "\\1", feature)),
                    ifelse(grepl("omic2_", feature),
                           as.numeric(sub(".*omic2_feature_(\\d+)", "\\1", feature)),
                           NA))
      ) %>%
      arrange(dataview, ID)
    
    load_omic_one <- result_load %>%
      mutate(
        dataview = ifelse(grepl("omic1_", feature), "omic.one",
                          ifelse(grepl("omic2_", feature), "omic.two", NA)),
        ID = ifelse(grepl("omic1_", feature),
                    as.numeric(sub(".*omic1_feature_(\\d+)", "\\1", feature)),
                    ifelse(grepl("omic2_", feature),
                           as.numeric(sub(".*omic2_feature_(\\d+)", "\\1", feature)),
                           NA))
      ) %>%
      arrange(dataview, ID) %>%
      filter(dataview == 'omic.one')
    
    load_omic_two <- result_load %>%
      mutate(
        dataview = ifelse(grepl("omic1_", feature), "omic.one",
                          ifelse(grepl("omic2_", feature), "omic.two", NA)),
        ID = ifelse(grepl("omic1_", feature),
                    as.numeric(sub(".*omic1_feature_(\\d+)", "\\1", feature)),
                    ifelse(grepl("omic2_", feature),
                           as.numeric(sub(".*omic2_feature_(\\d+)", "\\1", feature)),
                           NA))
      ) %>%
      arrange(dataview, ID) %>%
      filter(dataview == 'omic.two')
    
    weights_composite <- list(
      omic.one_weights = load_omic_one,
      omic.two_weights = load_omic_two,
      all_weights = result_load_all
    )
    
    # Factor scores separation
    result_score <- mofa_score
    
    result_score <- result_score %>%
      mutate(
        ID = ifelse(grepl("sample_", sample),
                    as.numeric(sub(".*sample_(\\d+)", "\\1", sample)), NA)
      ) %>%
      arrange(ID) #%>%
    
    scores_composite <- list(
      
      scores = result_score
    )
    
    mofa_result <- list(weights = weights_composite, scores = scores_composite)
    
    return(mofa_result)
  })
}

## ----------------------------------- 2.1 MOFA+ Loading Function ---------------------------------------
mofa_loading <- function(model_object, factor_num) {
  # Initialize an empty data frame for results
  loadings_df <- data.frame(sample = NULL)
  
  # Loop through each factor and extract scores
  for (i in 1:factor_num) {
    # Extract factor scores
    factor_loadings <- get_weights(model_object, factors = i, as.data.frame = T)
    # Keep sample and value columns
    factor_loadings <- factor_loadings[, c("feature", "value")]
    # Rename columns dynamically
    colnames(factor_loadings) <- c("feature", paste0("loading_MOFA", i))
    # Merge the current factor scores with the main scores data frame
    if (nrow(loadings_df) == 0) {
      loadings_df <- factor_loadings
    } else {
      loadings_df <- merge(loadings_df, factor_loadings, by = "feature")
    }
  }
  
  return(loadings_df)
}

## ------------------------------------ 2.2 MOFA+ Score Function ----------------------------------------

mofa_score <- function(model_object, factor_num) {
  # Initialize an empty data frame for results
  scores_df <- data.frame(sample = NULL)
  
  # Loop through each factor and extract scores
  for (i in 1:factor_num) {
    # Extract factor scores
    factor_scores <- get_factors(model_object, factors = i, as.data.frame = TRUE)
    # Keep sample and value columns
    factor_scores <- factor_scores[, c("sample", "value")]
    # Rename columns dynamically
    colnames(factor_scores) <- c("sample", paste0("score_MOFA", i))
    # Merge the current factor scores with the main scores data frame
    if (nrow(scores_df) == 0) {
      scores_df <- factor_scores
    } else {
      scores_df <- merge(scores_df, factor_scores, by = "sample")
    }
  }
  
  return(scores_df)
}

# --------------------------------------- MFA Function -----------------------------------------

## ----------------------------------- 3. MFA Main Function ---------------------------------------

# func_mfa <- function(data, num) {
#   
#   #set.seed(123)
#   BC_num = num
#   load_data <- data # dataset
#   
#   # Access variables from the shared environment
#   n_features_one <- 4000#global_env$n_features_one
#   n_features_two <- 3000#global_env$n_features_two
#   
#   first_omic <- load_data[, 1:n_features_one]
#   second_omic <- load_data[, (n_features_one + 1):(n_features_one + n_features_two)]
#   cdata = cbind(first_omic, second_omic)
#   
#   mfa_data <- cdata
#   # Create a list specifying which columns are quantitative or qualitative
#   
#   mfa_object = MFA(as.matrix(mfa_data), 
#                    group = c(n_features_one, n_features_two), 
#                    type = c("s","s"), 
#                    name.group = c("first.omic", "second.omic"),
#                    graph = FALSE)
#   
#   
#   mfa_loading <- mfa_loading(mfa_object, BC_num)
#   mfa_score <- mfa_score(mfa_object, BC_num)
#   
#   # Loading separation per factor
#   
#   result_load <- mfa_loading
#   
#   result_load_all <- result_load %>%
#     mutate(
#       dataview = ifelse(grepl("omic1_", feature), "omic.one",
#                         ifelse(grepl("omic2_", feature), "omic.two", NA)),
#       ID = ifelse(grepl("omic1_", feature),
#                   as.numeric(sub(".*omic1_feature_(\\d+)", "\\1", feature)),
#                   ifelse(grepl("omic2_", feature),
#                          as.numeric(sub(".*omic2_feature_(\\d+)", "\\1", feature)),
#                          NA))
#     ) %>%
#     arrange(dataview, ID) 
#   
#   load_omic_one <- result_load %>%
#     mutate(
#       dataview = ifelse(grepl("omic1_", feature), "omic.one",
#                         ifelse(grepl("omic2_", feature), "omic.two", NA)),
#       ID = ifelse(grepl("omic1_", feature),
#                   as.numeric(sub(".*omic1_feature_(\\d+)", "\\1", feature)),
#                   ifelse(grepl("omic2_", feature),
#                          as.numeric(sub(".*omic2_feature_(\\d+)", "\\1", feature)),
#                          NA))
#     ) %>%
#     arrange(dataview, ID) %>%
#     filter(dataview == 'omic.one')
#   
#   # Calculate the standard deviation of gfa_load$loading_GFA1 where signal_b == FALSE
#   indices_features.1a = simulated_data[[paste0("iteration_", i)]]$indices_features.1[1]
#   
#   indices_features.1a_vector <- unlist(indices_features.1a)  # Convert to atomic vector
#   indices_features.1a_sorted <- sort(indices_features.1a_vector)
#   
#   indices_features.1b = simulated_data[[paste0("iteration_", i)]]$indices_features.1[2]
#   
#   indices_features.1b_vector <- unlist(indices_features.1b)  # Convert to atomic vector
#   indices_features.1b_sorted <- sort(indices_features.1b_vector)
#   
#   in_range.b <- random.data.a$ID %in% indices_features.1b_sorted
#   in_range.a <- random.data.a$ID %in% indices_features.1a_sorted
#   
#   simulated_features_b <- data.frame(feature = load_omic_one$feature, signal_a = in_range.a, signal_b = in_range.b)
#   mfa_omic_one_load <- merge(x=load_omic_one,y=simulated_features_b, by="feature", all = TRUE)
#   
#   var_false <- var(mfa_omic_one_load$loading_MFA1[mfa_omic_one_load$signal_a == FALSE & mfa_omic_one_load$signal_b == FALSE], na.rm = TRUE)
#   mean_false <- mean(mfa_omic_one_load$loading_MFA1[mfa_omic_one_load$signal_a == FALSE & mfa_omic_one_load$signal_b == FALSE], na.rm = TRUE)
#   
#   # Suppress elevated factor
#   mfa_omic_one_load$loading_MFA1 <- ifelse(
#     mfa_omic_one_load$signal_b == TRUE,
#     rnorm(length(mfa_omic_one_load$loading_MFA1[mfa_omic_one_load$signal_b == TRUE]), mean_false, var_false), # Action if TRUE
#     mfa_omic_one_load$loading_MFA1# Action if FALSE
#   )
#   
#   # FACTOR 2
#   # Calculate the standard deviation of gfa_load$loading_GFA1 where signal_b == FALSE
#   indices_features.1a = simulated_data[[paste0("iteration_", i)]]$indices_features.1[1]
#   
#   indices_features.1a_vector <- unlist(indices_features.1a)  # Convert to atomic vector
#   indices_features.1a_sorted <- sort(indices_features.1a_vector)
#   
#   indices_features.1b = simulated_data[[paste0("iteration_", i)]]$indices_features.1[2]
#   
#   indices_features.1b_vector <- unlist(indices_features.1b)  # Convert to atomic vector
#   indices_features.1b_sorted <- sort(indices_features.1b_vector)
#   
#   in_range.b <- random.data.a$ID %in% indices_features.1b_sorted
#   in_range.a <- random.data.a$ID %in% indices_features.1a_sorted
#   
#   simulated_features_a <- data.frame(feature = load_omic_one$feature, signal_a = in_range.a, signal_b = in_range.b)
#   mfa_omic_one_load <- merge(x=mfa_omic_one_load,y=simulated_features_a, by="feature", all = TRUE)
#   
#   var_false <- var(mfa_omic_one_load$loading_MFA2[mfa_omic_one_load$signal_a == FALSE & mfa_omic_one_load$signal_b == FALSE], na.rm = TRUE)
#   mean_false <- mean(mfa_omic_one_load$loading_MFA2[mfa_omic_one_load$signal_a == FALSE & mfa_omic_one_load$signal_b == FALSE], na.rm = TRUE)
#   
#   # Suppress elevated factor
#   mfa_omic_one_load$loading_MFA2 <- ifelse(
#     mfa_omic_one_load$signal_a == TRUE,
#     rnorm(length(mfa_omic_one_load$loading_MFA2[mfa_omic_one_load$signal_a == TRUE]), mean_false, var_false), # Action if TRUE
#     mfa_omic_one_load$loading_MFA2# Action if FALSE
#   )
#   
#   load_omic_two <- result_load %>%
#     mutate(
#       dataview = ifelse(grepl("omic1_", feature), "omic.one",
#                         ifelse(grepl("omic2_", feature), "omic.two", NA)),
#       ID = ifelse(grepl("omic1_", feature),
#                   as.numeric(sub(".*omic1_feature_(\\d+)", "\\1", feature)),
#                   ifelse(grepl("omic2_", feature),
#                          as.numeric(sub(".*omic2_feature_(\\d+)", "\\1", feature)),
#                          NA))
#     ) %>%
#     arrange(dataview, ID) %>%
#     filter(dataview == 'omic.two')
#   # - - - - - - - - - - - - - - - - - - - MFA SCORES - - - - - - - - - - - - - - - - - #
#   
#   # Factor scores separation
#   result_score <- mfa_score
#   
#   # Factor 1
#   
#   scores_factor_one <- result_score %>%
#     mutate(
#       ID = ifelse(grepl("sample_", sample),
#                   as.numeric(sub(".*sample_(\\d+)", "\\1", sample)), NA)) %>%
#     arrange(ID) %>%
#     select(ID, sample, score_MFA1)  # Use `select()` to choose specific columns
#   
#   # Calculate the standard deviation of gfa_load$loading_GFA1 where signal_b == FALSE
#   indices_samples.a = simulated_data[[paste0("iteration_", i)]]$indices_samples[1]
#   indices_samples.b = simulated_data[[paste0("iteration_", i)]]$indices_samples[2]
#   
#   indices_samples.a_vector <- unlist(indices_samples.a)  # Convert to atomic vector
#   indices_samples.a_sorted <- sort(indices_samples.a_vector)
#   indices_samples.b_vector <- unlist(indices_samples.b)  # Convert to atomic vector
#   indices_samples.b_sorted <- sort(indices_samples.b_vector)
#   
#   in_range_sample_a <- scores_factor_one$ID %in% indices_samples.a_sorted
#   in_range_sample_b <- scores_factor_one$ID %in% indices_samples.b_sorted
#   
#   
#   simulated_samples_a <- data.frame(sample = scores_factor_one$sample, signal_a = in_range_sample_a)
#   simulated_samples_b <- data.frame(sample = scores_factor_one$sample, signal_b = in_range_sample_b)
#   
#   mfa_omic_one_score <- merge(merge(x=scores_factor_one,y=simulated_samples_a, by="sample", all = TRUE), simulated_samples_b,  by="sample", all = TRUE)
#   
#   var_false <- var(
#     mfa_omic_one_score$score_MFA1[mfa_omic_one_score$signal_a == FALSE & mfa_omic_one_score$signal_b == FALSE], 
#     na.rm = TRUE
#   )
#   
#   mean_false <- mean(
#     mfa_omic_one_score$score_MFA1[mfa_omic_one_score$signal_a == FALSE & mfa_omic_one_score$signal_b == FALSE], 
#     na.rm = TRUE
#   )
#   
#   # Suppress elevated factor
#   mfa_omic_one_score$score_MFA1 <- ifelse(
#     mfa_omic_one_score$signal_b == TRUE,
#     rnorm(length(mfa_omic_one_score$score_MFA1[mfa_omic_one_score$signal_b == TRUE]), mean_false, var_false), # Action if TRUE
#     mfa_omic_one_score$score_MFA1# Action if FALSE
#   )
#   
#   # Factor 2
#   
#   scores_factor_two <- result_score %>%
#     mutate(
#       ID = ifelse(grepl("sample_", sample),
#                   as.numeric(sub(".*sample_(\\d+)", "\\1", sample)), NA)) %>%
#     arrange(ID) %>%
#     select(ID, sample, score_MFA2)  # Use `select()` to choose specific columns
#   
#   mfa_omic_two_score <- merge(merge(x=scores_factor_two,y=simulated_samples_a, by="sample", all = TRUE),simulated_samples_b,  by="sample", all = TRUE)
#   
#   var_false <- var(
#     mfa_omic_two_score$score_MFA2[mfa_omic_two_score$signal_a == FALSE & mfa_omic_two_score$signal_b == FALSE], 
#     na.rm = TRUE
#   )
#   
#   mean_false <- mean(
#     mfa_omic_two_score$score_MFA2[mfa_omic_two_score$signal_a == FALSE & mfa_omic_two_score$signal_b == FALSE], 
#     na.rm = TRUE
#   )
#   
#   # Suppress elevated factor
#   mfa_omic_two_score$score_MFA2 <- ifelse(
#     mfa_omic_two_score$signal_a == TRUE,
#     rnorm(length(mfa_omic_two_score$score_MFA2[mfa_omic_two_score$signal_a == TRUE]), mean_false, var_false), # Action if TRUE
#     mfa_omic_two_score$score_MFA2# Action if FALSE
#   )
#   
#   # Select only specific columns from Y
#   mfa_omic_two_score_selected <- mfa_omic_two_score %>%
#     select(sample, score_MFA2)  # Replace 'column1', 'column2' with desired column names
#   
#   # Perform the merge, keeping all columns from X
#   factor_scores <- merge(
#     mfa_omic_one_score, 
#     mfa_omic_two_score_selected, 
#     by = "sample", 
#     all.x = TRUE
#   )
#   
#   
#   weights_composite <- list(omic.one_weights = mfa_omic_one_load, omic.two_weights = load_omic_two, original_weights = result_load_all)
#   scores_composite <- list(factor_scores = factor_scores, original_scores = result_score)
#   
#   mfa_result <- list(weights = weights_composite, scores = scores_composite)
#   
#   return(mfa_result)
# }
func_mfa <- function(data, num, sim_object, i) {
  mfa_seed <- 60667#sample(1e6, 1)
  set.seed(mfa_seed)
  
  # Saved key variables to global environment
  features_omic1 <- global_env$features_omic1
  features_omic2 <- global_env$features_omic2
  indices_features.OMIC1.A <- global_env$indices_features.OMIC1.A 
  indices_features.OMIC1.B <- global_env$indices_features.OMIC1.B  
  indices_features.OMIC2.A <- global_env$indices_features.OMIC2.A
  indices_samples.1A <- global_env$indices_samples.1A  
  indices_samples.2B <- global_env$indices_samples.2B 
  true_load_omic1_FACTOR1 <- global_env$true_load_omic1_FACTOR1
  true_load_omic1_FACTOR2 <- global_env$true_load_omic1_FACTOR2
  true_load_omic2_FACTOR1 <- global_env$true_load_omic2_FACTOR1
  true_load_omic2_FACTOR2 <- global_env$true_load_omic2_FACTOR2
  
  # Set parameters
  BC_num <- num
  load_data <- data
  
  # Feature split sizes
  n_features_one <- 4000
  n_features_two <- 3000
  
  # Split omics
  first_omic <- load_data[, 1:n_features_one]
  second_omic <- load_data[, (n_features_one + 1):(n_features_one + n_features_two)]
  cdata <- cbind(first_omic, second_omic)
  
  # Run MFA
  mfa_object <- MFA(as.matrix(cdata), 
                    group = c(n_features_one, n_features_two), 
                    type = c("s", "s"), 
                    name.group = c("first.omic", "second.omic"),
                    graph = FALSE)
  
  # Extract loadings and scores
  result_load <- mfa_loading(mfa_object, BC_num)
  result_score <- mfa_score(mfa_object, BC_num)
  
  # Annotate features with dataview and ID
  result_load_all <- result_load %>%
    mutate(
      dataview = ifelse(grepl("omic1_", feature), "omic.one",
                        ifelse(grepl("omic2_", feature), "omic.two", NA)),
      ID = ifelse(grepl("omic1_", feature),
                  as.numeric(sub(".*omic1_feature_(\\d+)", "\\1", feature)),
                  ifelse(grepl("omic2_", feature),
                         as.numeric(sub(".*omic2_feature_(\\d+)", "\\1", feature)),
                         NA))
    ) %>%
    arrange(dataview, ID)
  
  # Subset omic.one with all MFA loadings intact
  load_omic_one <- result_load_all %>% filter(dataview == "omic.one")
  
  # Load signal annotations
  indices_features.1a <- unlist(sim_object[[paste0("iteration_", i)]]$indices_features.1[1])
  indices_features.1b <- unlist(sim_object[[paste0("iteration_", i)]]$indices_features.1[2])
  
  in_range.a <- load_omic_one$ID %in% indices_features.1a
  in_range.b <- load_omic_one$ID %in% indices_features.1b
  
  sim_feat_labels <- data.frame(feature = load_omic_one$feature,
                                signal_a = in_range.a,
                                signal_b = in_range.b)
  
  mfa_omic_one_load <- merge(load_omic_one, sim_feat_labels, by = "feature", all = TRUE)
  
  # Suppress MFA1 loadings in signal_b
  var_false <- var(mfa_omic_one_load$loading_MFA1[!mfa_omic_one_load$signal_a & !mfa_omic_one_load$signal_b], na.rm = TRUE)
  mean_false <- mean(mfa_omic_one_load$loading_MFA1[!mfa_omic_one_load$signal_a & !mfa_omic_one_load$signal_b], na.rm = TRUE)
  
  mfa_omic_one_load$loading_MFA1 <- ifelse(
    mfa_omic_one_load$signal_b,
    rnorm(sum(mfa_omic_one_load$signal_b, na.rm = TRUE), mean_false, var_false),
    mfa_omic_one_load$loading_MFA1
  )
  
  # Suppress MFA2 loadings in signal_a
  var_false <- var(mfa_omic_one_load$loading_MFA2[!mfa_omic_one_load$signal_a & !mfa_omic_one_load$signal_b], na.rm = TRUE)
  mean_false <- mean(mfa_omic_one_load$loading_MFA2[!mfa_omic_one_load$signal_a & !mfa_omic_one_load$signal_b], na.rm = TRUE)
  
  mfa_omic_one_load$loading_MFA2 <- ifelse(
    mfa_omic_one_load$signal_a,
    rnorm(sum(mfa_omic_one_load$signal_a, na.rm = TRUE), mean_false, var_false),
    mfa_omic_one_load$loading_MFA2
  )
  
  # Load omic.two
  load_omic_two <- result_load_all %>% filter(dataview == "omic.two")
  
  # Load signal annotations
  indices_features.2a <- unlist(sim_object[[paste0("iteration_", i)]]$indices_features.2[1])
  
  in_range.a <- load_omic_two$ID %in% indices_features.2a
  
  sim_feat_labels2 <- data.frame(feature = load_omic_two$feature,
                                 signal_a = in_range.a)
  
  mfa_omic_two_load <- merge(load_omic_two, sim_feat_labels2, by = "feature", all = TRUE)
  
  # -- Process scores --
  # Factor 1
  scores_factor_one <- result_score %>%
    mutate(ID = as.numeric(sub(".*sample_(\\d+)", "\\1", sample))) %>%
    arrange(ID) %>%
    select(ID, sample, score_MFA1)
  
  indices_samples.a <- unlist(sim_object[[paste0("iteration_", i)]]$indices_samples[1])
  indices_samples.b <- unlist(sim_object[[paste0("iteration_", i)]]$indices_samples[2])
  
  in_range_sample_a <- scores_factor_one$ID %in% indices_samples.a
  in_range_sample_b <- scores_factor_one$ID %in% indices_samples.b
  
  sim_samp_labels <- data.frame(sample = scores_factor_one$sample,
                                signal_a = in_range_sample_a,
                                signal_b = in_range_sample_b)
  
  mfa_omic_one_score <- merge(scores_factor_one, sim_samp_labels, by = "sample", all = TRUE)
  
  var_false <- var(mfa_omic_one_score$score_MFA1[!mfa_omic_one_score$signal_a & !mfa_omic_one_score$signal_b], na.rm = TRUE)
  mean_false <- mean(mfa_omic_one_score$score_MFA1[!mfa_omic_one_score$signal_a & !mfa_omic_one_score$signal_b], na.rm = TRUE)
  
  mfa_omic_one_score$score_MFA1 <- ifelse(
    mfa_omic_one_score$signal_b,
    rnorm(sum(mfa_omic_one_score$signal_b, na.rm = TRUE), mean_false, var_false),
    mfa_omic_one_score$score_MFA1
  )
  
  # Factor 2
  scores_factor_two <- result_score %>%
    mutate(ID = as.numeric(sub(".*sample_(\\d+)", "\\1", sample))) %>%
    arrange(ID) %>%
    select(ID, sample, score_MFA2)
  
  mfa_omic_two_score <- merge(scores_factor_two, sim_samp_labels, by = "sample", all = TRUE)
  
  var_false <- var(mfa_omic_two_score$score_MFA2[!mfa_omic_two_score$signal_a & !mfa_omic_two_score$signal_b], na.rm = TRUE)
  mean_false <- mean(mfa_omic_two_score$score_MFA2[!mfa_omic_two_score$signal_a & !mfa_omic_two_score$signal_b], na.rm = TRUE)
  
  mfa_omic_two_score$score_MFA2 <- ifelse(
    mfa_omic_two_score$signal_a,
    rnorm(sum(mfa_omic_two_score$signal_a, na.rm = TRUE), mean_false, var_false),
    mfa_omic_two_score$score_MFA2
  )
  
  # Combine scores
  factor_scores <- merge(
    mfa_omic_one_score,
    mfa_omic_two_score %>% select(sample, score_MFA2),
    by = "sample", all.x = TRUE
  )
  
  # Final composite
  weights_composite <- list(
    omic.one_weights = mfa_omic_one_load,
    omic.two_weights = mfa_omic_two_load,
    original_weights = result_load_all
  )
  scores_composite <- list(
    factor_scores = factor_scores,
    original_scores = result_score
  )
  
  return(list(weights = weights_composite, scores = scores_composite))
}

## ----------------------------------- 3.1 MFA Loading Function ---------------------------------------

mfa_loading <- function(mfa_object, BC_num) {
  # Check if BC_num is greater than 1
  if (BC_num > 1) {
    # Initialize an empty data frame for results
    loadings_mfa_df <- data.frame()
    
    # Extract loadings for each BC_num and combine into a single data frame
    for (i in 1:BC_num) {
      loading_column <- mfa_object$quanti.var$contrib[, i]  # Extract the loading for the current factor
      column_name <- paste0("loading_MFA", i)  # Name for the column
      
      # Add the score column to the data frame
      if (ncol(loadings_mfa_df) == 0) {
        loadings_mfa_df <- data.frame(feature = rownames(mfa_object$quanti.var$contrib))  # Initialize with feature names
      }
      loadings_mfa_df[[column_name]] <- loading_column
    }
    mfa_df_loading = loadings_mfa_df
    return(mfa_df_loading)
  } else {
    # Handle single BC_num case
    loading_MFA <- mfa_object$quanti.var$contrib[, BC_num]
    mfa_df_loading <- as.data.frame(loading_MFA)
    mfa_df_loading$feature <- rownames(mfa_df_loading)
    return(mfa_df_loading)
  }
}

## ------------------------------------ 3.2 MFA Score Function ----------------------------------------

mfa_score <- function(mfa_object, BC_num) {
  # Check if BC_num is greater than 1
  if (BC_num > 1) {
    # Initialize an empty data frame for results
    scores_mfa_df <- data.frame()
    
    # Extract scores for each BC_num and combine into a single data frame
    for (i in 1:BC_num) {
      score_column <- mfa_object$ind$contrib[, i]  # Extract the score for the current factor
      column_name <- paste0("score_MFA", i)  # Name for the column
      
      # Add the score column to the data frame
      if (ncol(scores_mfa_df) == 0) {
        scores_mfa_df <- data.frame(sample = rownames(mfa_object$ind$contrib))  # Initialize with sample names
      }
      scores_mfa_df[[column_name]] <- score_column
    }
    mfa_df_score = scores_mfa_df
    return(mfa_df_score)
  } else {
    # Handle single BC_num case
    score_MFA <- mfa_object$ind$contrib[, BC_num]
    mfa_df_score <- as.data.frame(score_MFA)
    mfa_df_score$sample <- rownames(mfa_df_score)
    return(mfa_df_score)
  }
}

# --------------------------------------- GFA Function -----------------------------------------

## ----------------------------------- 3. GFA Main Function ---------------------------------------
func_gfa <- function(data, num, sim_object, i) {
  gfa_seed <- 60667#sample(1e6, 1)
  set.seed(gfa_seed)
  
  # Saved key variables to global environment
  features_omic1 <- global_env$features_omic1
  features_omic2 <- global_env$features_omic2
  indices_features.OMIC1.A <- global_env$indices_features.OMIC1.A 
  indices_features.OMIC1.B <- global_env$indices_features.OMIC1.B  
  indices_features.OMIC2.A <- global_env$indices_features.OMIC2.A
  indices_samples.1A <- global_env$indices_samples.1A  
  indices_samples.2B <- global_env$indices_samples.2B 
  true_load_omic1_FACTOR1 <- global_env$true_load_omic1_FACTOR1
  true_load_omic1_FACTOR2 <- global_env$true_load_omic1_FACTOR2
  true_load_omic2_FACTOR1 <- global_env$true_load_omic2_FACTOR1
  true_load_omic2_FACTOR2 <- global_env$true_load_omic2_FACTOR2
  
  # Feature configuration
  n_features_one <- 4000
  n_features_two <- 3000
  
  first_omic <- data[, 1:n_features_one]
  second_omic <- data[, (n_features_one + 1):(n_features_one + n_features_two)]
  cdata <- as.matrix(cbind(first_omic, second_omic))
  
  gfa_dt <- as.data.frame(cdata)
  merged_GFA_data <- list(t(gfa_dt))
  
  model_option <- getDefaultOpts()
  model_option$iter.max <- 1000
  model_option$iter.burnin <- 10
  
  # Run GFA model (you must define this properly in your actual pipeline)
  gfa_object <- gfa(merged_GFA_data, K = 10, opts = model_option)
  
  # Extract loadings and scores
  result_load <- gfa_loading(gfa_object, num)
  result_score <- gfa_score(gfa_object, num)
  
  # True loading features
  df_factor_loading <- data.frame(
    feature = c(features_omic1, features_omic2),
    #beta1   = c(true_load_omic1_FACTOR1, rep(NA, length(true_load_omic2_FACTOR1))),
    #delta1  = c(rep(NA, length(true_load_omic1_FACTOR1)), true_load_omic2_FACTOR1),
    true_loading_F1 = c(true_load_omic1_FACTOR1, true_load_omic2_FACTOR1),
    true_loading_F2 = c(true_load_omic1_FACTOR2, true_load_omic2_FACTOR2)
  )
  
  # Annotate features
  result_load_all <- result_load %>%
    mutate(
      dataview = ifelse(grepl("omic1_", feature), "omic.one",
                        ifelse(grepl("omic2_", feature), "omic.two", NA)),
      ID = ifelse(grepl("omic1_", feature),
                  as.numeric(sub(".*omic1_feature_(\\d+)", "\\1", feature)),
                  ifelse(grepl("omic2_", feature),
                         as.numeric(sub(".*omic2_feature_(\\d+)", "\\1", feature)),
                         NA))
    ) %>%
    arrange(dataview, ID)
  
  load_omic_one <- result_load_all %>% filter(dataview == "omic.one")
  
  # Signal annotations for loadings
  indices_features.1a <- unlist(sim_object[[paste0("iteration_", i)]]$indices_features.1[1])
  indices_features.1b <- unlist(sim_object[[paste0("iteration_", i)]]$indices_features.1[2])
  
  in_range.a <- load_omic_one$ID %in% indices_features.1a
  in_range.b <- load_omic_one$ID %in% indices_features.1b
  
  sim_features <- data.frame(feature = load_omic_one$feature, signal_a = in_range.a, signal_b = in_range.b)
  gfa_omic_one_load <- merge(load_omic_one, sim_features, by = "feature", all = TRUE)
  
  # Suppress GFA1 (signal_b)
  var_false <- var(gfa_omic_one_load$loading_GFA1[!gfa_omic_one_load$signal_a & !gfa_omic_one_load$signal_b], na.rm = TRUE)
  mean_false <- mean(gfa_omic_one_load$loading_GFA1[!gfa_omic_one_load$signal_a & !gfa_omic_one_load$signal_b], na.rm = TRUE)
  
  gfa_omic_one_load$loading_GFA1 <- ifelse(
    gfa_omic_one_load$signal_b,
    rnorm(sum(gfa_omic_one_load$signal_b, na.rm = TRUE), mean_false, var_false),
    gfa_omic_one_load$loading_GFA1
  )
  
  # Suppress GFA2 (signal_a)
  var_false <- var(gfa_omic_one_load$loading_GFA2[!gfa_omic_one_load$signal_a & !gfa_omic_one_load$signal_b], na.rm = TRUE)
  mean_false <- mean(gfa_omic_one_load$loading_GFA2[!gfa_omic_one_load$signal_a & !gfa_omic_one_load$signal_b], na.rm = TRUE)
  
  gfa_omic_one_load$loading_GFA2 <- ifelse(
    gfa_omic_one_load$signal_a,
    rnorm(sum(gfa_omic_one_load$signal_a, na.rm = TRUE), mean_false, var_false),
    gfa_omic_one_load$loading_GFA2
  )
  
  # Omic two
  load_omic_two <- result_load_all %>% filter(dataview == "omic.two")
  
  # Load signal annotations
  indices_features.2a <- unlist(sim_object[[paste0("iteration_", i)]]$indices_features.2[1])
  
  in_range.a <- load_omic_two$ID %in% indices_features.2a
  
  sim_feat_labels2 <- data.frame(feature = load_omic_two$feature,
                                 signal_a = in_range.a)
  
  gfa_omic_two_load <- merge(load_omic_two, sim_feat_labels2, by = "feature", all = TRUE)
  
  # Process scores
  scores_factor_one <- result_score %>%
    mutate(ID = as.numeric(sub(".*sample_(\\d+)", "\\1", sample))) %>%
    arrange(ID) %>%
    select(ID, sample, score_GFA1)
  
  indices_samples.a <- unlist(sim_object[[paste0("iteration_", i)]]$indices_samples[1])
  indices_samples.b <- unlist(sim_object[[paste0("iteration_", i)]]$indices_samples[2])
  
  in_range_sample_a <- scores_factor_one$ID %in% indices_samples.a
  in_range_sample_b <- scores_factor_one$ID %in% indices_samples.b
  
  sim_samples <- data.frame(sample = scores_factor_one$sample, signal_a = in_range_sample_a, signal_b = in_range_sample_b)
  gfa_omic_one_score <- merge(scores_factor_one, sim_samples, by = "sample", all = TRUE)
  
  # Suppress GFA1 score
  var_false <- var(gfa_omic_one_score$score_GFA1[!gfa_omic_one_score$signal_a & !gfa_omic_one_score$signal_b], na.rm = TRUE)
  mean_false <- mean(gfa_omic_one_score$score_GFA1[!gfa_omic_one_score$signal_a & !gfa_omic_one_score$signal_b], na.rm = TRUE)
  
  gfa_omic_one_score$score_GFA1 <- ifelse(
    gfa_omic_one_score$signal_b,
    rnorm(sum(gfa_omic_one_score$signal_b, na.rm = TRUE), mean_false, var_false),
    gfa_omic_one_score$score_GFA1
  )
  
  # Suppress GFA2 score
  scores_factor_two <- result_score %>%
    mutate(ID = as.numeric(sub(".*sample_(\\d+)", "\\1", sample))) %>%
    arrange(ID) %>%
    select(ID, sample, score_GFA2)
  
  gfa_omic_two_score <- merge(scores_factor_two, sim_samples, by = "sample", all = TRUE)
  
  var_false <- var(gfa_omic_two_score$score_GFA2[!gfa_omic_two_score$signal_a & !gfa_omic_two_score$signal_b], na.rm = TRUE)
  mean_false <- mean(gfa_omic_two_score$score_GFA2[!gfa_omic_two_score$signal_a & !gfa_omic_two_score$signal_b], na.rm = TRUE)
  
  gfa_omic_two_score$score_GFA2 <- ifelse(
    gfa_omic_two_score$signal_a,
    rnorm(sum(gfa_omic_two_score$signal_a, na.rm = TRUE), mean_false, var_false),
    gfa_omic_two_score$score_GFA2
  )
  
  # Merge scores
  factor_scores <- merge(
    gfa_omic_one_score,
    gfa_omic_two_score %>% select(sample, score_GFA2),
    by = "sample", all.x = TRUE
  )
  
  weights_composite <- list(
    omic.one_weights = gfa_omic_one_load,
    omic.two_weights = gfa_omic_two_load,
    original_weights = result_load_all
  )
  scores_composite <- list(
    factor_scores = factor_scores,
    original_scores = result_score
  )
  
  return(list(weights = weights_composite, scores = scores_composite))
}

## ------------------------------------ 3.2 GFA Score Function ----------------------------------------
gfa_score <- function(gfa_object, BC_num) {
  # Initialize an empty data frame for results
  scores_gfa_df <- data.frame(sample = NULL)
  
  # Loop through each factor and extract scores
  for (i in 1:BC_num) {
    # Extract factor scores
    factor_scores <- as.data.frame(gfa_object$W)[i]
    factor_scores$feature <- rownames(factor_scores)
    # Rename columns dynamically
    colnames(factor_scores) <- c(paste0("score_GFA", i), "sample")
    # Merge the current factor scores with the main scores data frame
    if (nrow(scores_gfa_df) == 0) {
      scores_gfa_df <- factor_scores
    } else {
      scores_gfa_df <- merge(scores_gfa_df, factor_scores, by = "sample")
    }
  }
  
  return(scores_gfa_df)
}

## ------------------------------------ 3.2 GFA Score Function ----------------------------------------


gfa_loading <- function(gfa_object, BC_num) {
  # Initialize an empty data frame for results
  loadings_gfa_df <- data.frame(sample = NULL)
  
  # Loop through each factor and extract scores
  for (i in 1:BC_num) {
    # Extract factor scores
    factor_loadings <- as.data.frame(gfa_object$X)[i]
    factor_loadings$feature <- rownames(factor_loadings)
    # Rename columns dynamically
    colnames(factor_loadings) <- c(paste0("loading_GFA", i), "feature")
    # Merge the current factor scores with the main scores data frame
    if (nrow(loadings_gfa_df) == 0) {
      loadings_gfa_df <- factor_loadings
    } else {
      loadings_gfa_df <- merge(loadings_gfa_df, factor_loadings, by = "feature")
    }
  }
  
  return(loadings_gfa_df)
}


# -------------------------- JI:Define the function to compute the Jaccard index -----------------------------
jaccard.index.sim <- function(Column1, Column2) {
  
  if (length(Column1) != length(Column2)) {
    stop("Columns must have the same length.")
  }
  similar_count <- sum(Column1 == 1 & Column2 == 1)
  sum_column1 <- sum(Column1 == 1)
  sum_column2 <- sum(Column2 == 1)
  jaccard_index <- similar_count / (sum_column1 + sum_column2 - similar_count)
  ji = ifelse(jaccard_index == -1, 1, jaccard_index)
  return(ji)
}

# --------------------------- Metrics -------------------------------------------------------------------------

CompareMetrics <- function(data, signal_index, methods) {
  
  # Helper function to calculate metrics
  calculate_metrics <- function(actual, predicted) {
    # Generate confusion matrix
    cm <- table(Predicted = predicted, Actual = actual)
    
    # Pad the confusion matrix to ensure '0' and '1' rows/columns exist
    if (!("0" %in% rownames(cm))) cm <- rbind(cm, "0" = rep(0, ncol(cm)))
    if (!("1" %in% rownames(cm))) cm <- rbind(cm, "1" = rep(0, ncol(cm)))
    if (!("0" %in% colnames(cm))) cm <- cbind(cm, "0" = rep(0, nrow(cm)))
    if (!("1" %in% colnames(cm))) cm <- cbind(cm, "1" = rep(0, nrow(cm)))
    
    # Extract values with safe handling
    TP <- ifelse("1" %in% rownames(cm) && "1" %in% colnames(cm), cm["1", "1"], 0)
    TN <- ifelse("0" %in% rownames(cm) && "0" %in% colnames(cm), cm["0", "0"], 0)
    FP <- ifelse("1" %in% rownames(cm) && "0" %in% colnames(cm), cm["1", "0"], 0)
    FN <- ifelse("0" %in% rownames(cm) && "1" %in% colnames(cm), cm["0", "1"], 0)
    
    # Calculate metrics safely
    sensitivity <- ifelse((TP + FN) == 0, 0, TP / (TP + FN))
    specificity <- ifelse((TN + FP) == 0, 0, TN / (TN + FP))
    accuracy <- ifelse(sum(cm) == 0, 0, (TP + TN) / sum(cm))
    
    # AUC: Check if actual has two levels
    auc_value <- ifelse(length(unique(actual)) > 1, auc(roc(actual, predicted)), NA)
    
    # Return as data frame
    return(data.frame(
      TP = TP, FP = FP, TN = TN, FN = FN,
      Sensitivity = sensitivity,
      Specificity = specificity,
      Accuracy = accuracy,
      AUC = auc_value
    ))
  }
  
  # Initialize lists for results
  result_true_list <- list()
  result_comparison_list <- list()
  
  ### TRUE Comparisons: Compare methods to the ground truth (signal_index)
  for (method in methods) {
    result_true_list[[method]] <- calculate_metrics(data[[signal_index]], data[[method]])
  }
  
  ### Pairwise Comparisons Between Methods
  for (i in 1:(length(methods) - 1)) {
    for (j in (i + 1):length(methods)) {
      name <- paste(methods[i], "vs", methods[j], sep = "_")
      result_comparison_list[[name]] <- calculate_metrics(data[[methods[i]]], data[[methods[j]]])
    }
  }
  
  # Combine results into separate data frames
  per_measures_true <- do.call(rbind, lapply(names(result_true_list), function(x) {
    cbind(Method = x, result_true_list[[x]])
  }))
  
  per_measures_comparison <- do.call(rbind, lapply(names(result_comparison_list), function(x) {
    cbind(Method_Comparison = x, result_comparison_list[[x]])
  }))
  
  # Return results as a list of two datasets
  return(list(
    per_measures_true = per_measures_true,
    per_measures_comparison = per_measures_comparison
  ))
}

# -------------------------- Actual simulation: Define the main simulation function -----------------------------
current_seed <- 694#sample(1e6, 1)
set.seed(current_seed)

# Create a specific environment to hold the variables
global_env <- new.env() # shared_env
# Assign parameters to global environment
#global_env$n_features_one <- n_features_one
#global_env$n_features_two <- n_features_two

# Initialize result lists
jaccard_results <- list()
jaccard_comparison_results <- list()
per_measures_results <- list()
dataset_output <- list()
datasets <- list()
dataset_output_corrected <- list()


# Cutoff 1: Based on the ratio of the top value of the features weights/loading or samples scores and the standard deviation of noise in the simulated data
varphi.one <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    varphi <- max(abs(numeric_values)) / sqrt(sigma)
    return(varphi)
  } else {
    return(NA)  # Return NA if not all values are numeric
  }
}

# varphi.one(score_dta$score_FABIA2)
# varphi.one(score_dta$score_MOFA2)
# varphi.one(score_dta$score_MFA2)
# varphi.one(score_dta$score_GFA2)

# Cutoff 2: Based on the ratio of the top value of the features weights/loading or samples scores and the variance of noise in the simulated data
varphi.two <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    varphi <- max(abs(numeric_values)) / (sigma/2)
    return(varphi)
  } else {
    return(NA)  # Return NA if not all values are numeric
  }
}
# 
# varphi.two(score_dta$score_FABIA2)
# varphi.two(score_dta$score_MOFA2)
# varphi.two(score_dta$score_MFA2)
# varphi.two(score_dta$score_GFA2)
# 
# varphi.two(loading_dta$loading_FABIA1)
# varphi.two(loading_dta$loading_MOFA1)

# Cutoff 3: Based on the ratio of the top 100 features weights/loading or samples scores and the standard deviation of noise in the simulated data
varphi.three <- function(loading_or_score) {
  # Obtain the name of the vector
  vector_name <- deparse(substitute(loading_or_score))
  
  numeric_values <- as.numeric(unlist(loading_or_score))
  
  if (any(grepl("loading", vector_name)==TRUE)) {
    if (all(!is.na(numeric_values))) {
      x <- abs(loading_or_score)
      top.x <- mean(head(sort(x, decreasing = TRUE), 100))
      varphi <- top.x / sqrt(sigma)
      return(varphi)
    }
  }
  
  else if (any(grepl("score", vector_name)==TRUE)) {
    if (all(!is.na(numeric_values))) {
      x <- abs(loading_or_score)
      top.x <- mean(head(sort(x, decreasing = TRUE), 9))
      varphi <- top.x / sqrt(sigma)
      return(varphi)
    }
  }
  
  #stop("Neither 'loading' nor 'score' found in the method argument.")
}

# Cutoff 4: Based on the ratio of the top 100 features weights/loading or samples scores and the variance of noise in the simulated data
varphi.four <- function(loading_or_score) {
  # Obtain the name of the vector
  vector_name <- deparse(substitute(loading_or_score))
  
  numeric_values <- as.numeric(unlist(loading_or_score))
  
  if (any(grepl("loading", vector_name)==TRUE)) {
    if (all(!is.na(numeric_values))) {
      x <- abs(numeric_values)
      top.x <- mean(head(sort(x, decreasing = TRUE), 100))
      varphi <- top.x / (sigma/2)
      return(varphi)
    }
  }
  
  else if (any(grepl("score", vector_name)==TRUE)) {
    if (all(!is.na(numeric_values))) {
      x <- abs(numeric_values)
      top.x <- mean(head(sort(x, decreasing = TRUE), 9))
      varphi <- top.x / (sigma/2)
      return(varphi)
    }
  }
  
  stop("Neither 'loading' nor 'score' found in the method argument.")
}

# Cutoff 5: Based on the true proportion of the samples or features with signal in the matrix
varphi.five <- function(loading_or_score, omic = NULL, factor = NULL, type = c("loading", "score")) {
  type <- match.arg(type)
  numeric_values <- as.numeric(unlist(loading_or_score))
  indices_features.OMIC1.A <- global_env$indices_features.OMIC1.A 
  indices_features.OMIC1.B <- global_env$indices_features.OMIC1.B  
  indices_features.OMIC2.A <- global_env$indices_features.OMIC2.A
  indices_samples.1A <- global_env$indices_samples.1A  
  indices_samples.2B <- global_env$indices_samples.2B
  n_features_one <- global_env$n_features_one
  n_features_two <- global_env$n_features_two
  n_samples <- global_env$n_samples
  
  if (type == "loading") {
    if (omic == 1 && factor == 1) {
      x <- 1 - (length(indices_features.OMIC1.A) / n_features_one)
      buffer = 0.000001
      numeric_values = numeric_values + buffer # add buffer
      varphi <- quantile(abs(numeric_values), probs = x)[[1]]
    } else if (omic == 1 && factor == 2) {
      buffer = 0.000001
      numeric_values = numeric_values + buffer # add buffer
      x <- 1 - (length(indices_features.OMIC1.B) / n_features_one)
      varphi <- quantile(abs(numeric_values), probs = x)[[1]]
    } else if (omic == 2 && factor == 1) {
      buffer = 0.000001
      numeric_values = numeric_values + buffer # add buffer
      y <- 1 - (length(indices_features.OMIC2.A) / n_features_two)
      varphi <- quantile(abs(numeric_values), probs = y)[[1]]
    }else if (omic == 2 && factor == 2) {
      buffer = 0.000001
      numeric_values = numeric_values + buffer # add buffer
      #y <- 1 - (length(indices_features.OMIC2.A) / n_features_two)
      varphi <- max(abs(numeric_values))#, probs = y)[[1]]
    } else {
      stop("Invalid omic/factor combination for loading.")
    }
    return(varphi)
  }
  
  if (type == "score") {
    if (factor == 1) {
      y <- 1 - (length(indices_samples.1A) / n_samples)
      varphi <- quantile(abs(numeric_values), probs = y)[[1]]
    } else if (factor == 2) {
      y <- 1 - (length(indices_samples.2B) / n_samples)
      varphi <- quantile(abs(numeric_values), probs = y)[[1]]
    } else {
      stop("Invalid factor value for score.")
    }
    return(varphi)
  }
  
  stop("Invalid type: must be either 'loading' or 'score'")
}

# varphi.five_multisim <- function(loading_or_score, omic = NULL, factor = NULL, type = c("loading", "score")) {
#   type <- match.arg(type)
#   numeric_values <- as.numeric(unlist(loading_or_score))
#   
#   if (type == "loading") {
#     if (method == 'multiple.factor'){
#       if (omic == 1 && factor == 1) {
#         x <- (1 - (length(indices_features.OMIC1.A) / n_features_one))
#         varphi <- quantile(abs(numeric_values), probs = x)[[1]]
#       }else if (omic == 1 && factor == 2) {
#         x <- (1 - (length(indices_features.OMIC1.B) / n_features_one))
#         varphi <- quantile(abs(numeric_values), probs = x)[[1]]
#       }else if (omic == 2 && factor == 1) {
#         y <- (1 - (length(indices_features.OMIC2.A) / n_features_two))
#         varphi <- quantile(abs(numeric_values), probs = y)[[1]]
#       } 
#     } else {
#       stop("Invalid method specified.")
#     }
#     return(varphi)
#   }
#   
#   if (type == "score") {
#     if (method =='multiple.factor') {
#       if(factor == 1){
#         y <- (1 - (length(indices_samples.1A) / n_samples))
#         varphi <- quantile(abs(numeric_values), probs = y)[[1]]
#       } else if(factor == 2){
#         y <- (1 - (length(indices_samples.2B) / n_samples))
#         varphi <- quantile(abs(numeric_values), probs = y)[[1]]
#       }
#     } else {
#       stop("Invalid method specified.")
#     }
#     return(varphi)
#   }
#   
#   stop("Invalid type: must be either 'loading' or 'score'")
# }

# Cutoff 6: Half the 0 to 1 range of the normalized values of loading and scores
normalize_columns <- function(df, columns, clip_percentiles = TRUE) {
  normalize <- function(x, lower_bound = NULL, upper_bound = NULL) {
    if (!is.null(lower_bound) && !is.null(upper_bound)) {
      x <- pmin(pmax(x, lower_bound), upper_bound)  # Clip values to bounds
    }
    x_min <- min(x, na.rm = TRUE)
    x_max <- max(x, na.rm = TRUE)
    if (x_max - x_min == 0) {
      return(rep(0.5, length(x)))  # Assign 0.5 if all values are identical
    }
    (x - x_min) / (x_max - x_min)
  }
  
  df_normalized <- df
  
  for (col in columns) {
    if (clip_percentiles) {
      lower_bound <- quantile(df[[col]], 0.01, na.rm = TRUE)  # 1st percentile
      upper_bound <- quantile(df[[col]], 0.99, na.rm = TRUE)  # 99th percentile
    } else {
      lower_bound <- NULL
      upper_bound <- NULL
    }
    df_normalized[[col]] <- normalize(df[[col]], lower_bound, upper_bound)
  }
  
  return(df_normalized)
}

varphi.six <- function(loading_or_score) {
  # Convert to numeric
  numeric_values <- as.numeric(unlist(loading_or_score))
  
  if (all(!is.na(numeric_values))) {
    # Normalize values
    normalized_values <- normalize_columns(
      data.frame(values = numeric_values), 
      columns = "values"
    )$values
    
    # Define cutoff as half the range of normalized values
    varphi <- 0.5 * (max(normalized_values) - min(normalized_values))
    
    return(varphi)
  } else {
    return(NA)  # Return NA if not all values are numeric
  }
}

# Cutoff 7: 80% Quantile, more restrictive than lower quantiles.
varphi.seven <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    buffer = 0.000001
    numeric_values = numeric_values + buffer # add buffer
    varphi <- quantile(numeric_values, probs = 0.80)[[1]]
    return(varphi)
  } else {
    return(NA)  # Return NA if not all values are numeric
  }
}

# Cutoff 8: 85% Quantile, more restrictive than lower quantiles.
varphi.eight <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    buffer = 0.000001
    numeric_values = numeric_values + buffer # add buffer
    varphi <- quantile(numeric_values, probs = 0.85)[[1]]
    return(varphi)
  } else {
    return(NA)  # Return NA if not all values are numeric
  }
}
# Cutoff 9: 90% Quantile, more restrictive than lower quantiles.
varphi.nine <- function(loading_or_score) {
  numeric_values <- as.numeric(unlist(loading_or_score))
  if (all(!is.na(numeric_values))) {
    buffer = 0.000001
    numeric_values = numeric_values + buffer # add buffer
    varphi <- quantile(numeric_values, probs = 0.90)[[1]]
    return(varphi)
  } else {
    return(NA)  # Return NA if not all values are numeric
  }
}

# Cutoff 10: Perform k-means clustering (k = 2 for signal/noise)
varphi.ten <- function(loading_or_score) {
  # Obtain the name of the vector
  vector_name <- deparse(substitute(loading_or_score))
  
  numeric_values <- as.numeric(unlist(loading_or_score))
  
  if (any(grepl("loading", vector_name)==TRUE)) {
    if (all(!is.na(numeric_values))) {
      clusters <- kmeans(numeric_values, centers = 2)$cluster
      # Assign the higher mean cluster as "signal"
      varphi <- mean(numeric_values[clusters == which.max(tapply(numeric_values, clusters, mean))])
      return(varphi)
    }
  }
  
  else if (any(grepl("score", vector_name)==TRUE)) {
    if (all(!is.na(numeric_values))) {
      clusters <- kmeans(numeric_values, centers = 2)$cluster
      # Assign the higher mean cluster as "signal"
      varphi <- mean(numeric_values[clusters == which.max(tapply(numeric_values, clusters, mean))])
      return(varphi)
    }
  }
  
  #stop("Neither 'loading' nor 'score' found in the method argument.")
}

# Cutoff 11:Based on Gaussian Mixture Model (GMM) clustering (G = 2 for signal/noise)
varphi.eleven <- function(loading_or_score) {
  # Obtain the name of the vector
  vector_name <- deparse(substitute(loading_or_score))
  
  # Convert input to numeric vector
  numeric_values <- as.numeric(unlist(loading_or_score))
  
  # Check if the vector contains "loading" or "score" in its name
  if (any(grepl("loading", vector_name)) || any(grepl("score", vector_name))) {
    if (all(!is.na(numeric_values))) {
      # Install and load the mclust package
      #if (!requireNamespace("mclust", quietly = TRUE)) install.packages("mclust")
      #library(mclust)
      
      # Fit a Gaussian Mixture Model with 2 components
      gmm <- Mclust(numeric_values, G = 2)
      
      # Determine the cluster with the higher mean
      signal_cluster <- which.max(gmm$parameters$mean)
      
      # Calculate the cutoff as the mean of the higher mean cluster
      varphi <- mean(numeric_values[gmm$classification == signal_cluster])
      
      return(varphi)
    }
  }
}

# Cutoff 12: Based on the rolling mean and binary comparison
varphi.twelve <- function(loading_or_score) {
  # Obtain the name of the vector
  vector_name <- deparse(substitute(loading_or_score))
  
  # Convert input to numeric vector
  numeric_values <- as.numeric(unlist(loading_or_score))
  
  # Check if the vector contains "loading" or "score" in its name
  if (any(grepl("loading", vector_name)) || any(grepl("score", vector_name))) {
    if (all(!is.na(numeric_values))) {
      # Install and load zoo for rolling mean
      #if (!requireNamespace("zoo", quietly = TRUE)) install.packages("zoo")
      #library(zoo)
      
      # Calculate rolling mean with a window size of 10
      rolling_mean <- rollmean(numeric_values, k = 10, fill = NA, align = "center")
      
      # Replace NA with 0 for binary comparison
      rolling_mean[is.na(rolling_mean)] <- 0
      
      # Convert to binary: 1 if the value exceeds the rolling mean, 0 otherwise
      binary_vector <- ifelse(numeric_values > rolling_mean, 1, 0)
      
      # Filter numeric values classified as signal (binary 1)
      signal_values <- numeric_values[binary_vector == 1]
      
      # Calculate the cutoff as the mean of the "signal" values
      if (length(signal_values) > 0) {
        varphi <- mean(signal_values)
        return(varphi=varphi)
      } else {
        warning("No signal values found above the rolling mean.")
        return(NA)
      }
    }
  }
}

varphi_functions <- list(
  varphi.one = varphi.one, # top value / std dev
  varphi.two = varphi.two, # top value / variance
  varphi.three = varphi.three, # top x values / std dev
  varphi.four = varphi.four, # top x values / variance
  varphi.five = varphi.five, # prop of true signal
  varphi.six = varphi.six, # 0.5 for normalized values
  varphi.seven = varphi.seven, # 80% quantile
  varphi.eight = varphi.eight, # 85% quantile
  varphi.nine = varphi.nine, # 90% quantile
  varphi.ten = varphi.ten, # kmeans
  varphi.eleven = varphi.eleven, # GMM
  varphi.twelve = varphi.twelve # rolling mean
  
)

# extract_id_from_feature <- function(df, feature_col = "feature") {
#   df$ID <- as.numeric(sub(".*omic1_feature_(\\d+)", "\\1", df[[feature_col]]))
#   df
# }

extract_id_from_feature <- function(df, feature_col = "feature") {
  df$ID <- as.numeric(sub(".*omic[12]_feature_(\\d+)", "\\1", df[[feature_col]]))
  df
}

extract_id_from_sample <- function(df, sample_col = "sample") {
  df$ID <- as.numeric(sub(".*sample_(\\d+)", "\\1", df[[sample_col]]))
  df
}

run_with_seed <- function(expr, seed = NULL) {
  old_seed <- .Random.seed
  set.seed(seed)
  on.exit({ .Random.seed <<- old_seed }, add = TRUE)
  force(expr)
}

align_loading_columns_by_truth <- function(df, true_vector) {
  # Identify loading columns
  loading_cols <- names(df)[grepl("loading", names(df))]
  
  # Extract only loading columns
  loading_matrix <- df[, loading_cols, drop = FALSE]
  
  # Check dimensions
  if (nrow(loading_matrix) != length(true_vector)) {
    stop("Number of rows in loading columns must match length of true_vector.")
  }
  
  # Compute correlation with true vector
  cor_vals <- sapply(loading_matrix, function(col) cor(col, true_vector, use = "pairwise.complete.obs"))
  
  # Order columns by absolute correlation (descending)
  ordered_idx <- order(abs(cor_vals), decreasing = TRUE)
  ordered_cols <- loading_cols[ordered_idx]
  
  # Extract method name (e.g., "FABIA" from "loading_FABIA1")
  method_name <- sub(".*loading_([A-Za-z]+).*", "\\1", ordered_cols[1])
  
  # Rename columns based on match order
  new_names <- paste0("loading_", method_name, seq_along(ordered_cols))
  
  # Replace names in original df
  df_renamed <- df %>%
    dplyr::select(-all_of(loading_cols)) %>%
    dplyr::bind_cols(setNames(loading_matrix[, ordered_idx], new_names))
  
  return(df_renamed)
}

align_loading_columns_by_truth_overall <- function(weights_df, true_df, method_name) {
  # Merge by feature
  df_merged <- merge(true_df, weights_df, by = "feature")
  
  # Pull true values
  true_vector <- df_merged$true_loading_F1
  
  # Extract only numeric loading columns for this method
  loading_cols <- grep(paste0("^loading_", method_name, "[0-9]+$"), names(df_merged), value = TRUE)
  
  if (length(loading_cols) < 2) {
    stop(paste("Expected at least 2 loading columns for method", method_name))
  }
  
  # Ensure loading columns are numeric
  loading_matrix <- df_merged#[, loading_cols, drop = FALSE]
  for (col in loading_cols) {
    if (!is.numeric(loading_matrix[[col]])) {
      loading_matrix[[col]] <- as.numeric(as.character(loading_matrix[[col]]))
    }
  }
  
  # Compute correlations
  cor_vals <- sapply(loading_matrix[, loading_cols], function(col) cor(col, true_vector, use = "pairwise.complete.obs"))
  
  # Order by absolute correlation
  ordered_idx <- order(abs(cor_vals), decreasing = TRUE)
  ordered_cols <- loading_cols[ordered_idx]
  
  new_names <- paste0("matched_", method_name, seq_along(ordered_cols))
  old_names <- character(length(ordered_idx))
  for (j in seq_along(ordered_idx)) {
    old_names[j] <- paste0("loading_", method_name, ordered_idx[j])
  }
  
  # Rename the columns in df_merged directly
  names(df_merged)[match(old_names, names(df_merged))] <- new_names
  
  return(df_merged)
}

align_scores_columns_by_truth <- function(df, true_vector) {
  # Identify score columns
  score_cols <- names(df)[grepl("score", names(df))]
  
  # Extract only score columns
  score_matrix <- df[, score_cols, drop = FALSE]
  
  # Check dimensions
  if (nrow(score_matrix) != length(true_vector)) {
    stop("Number of rows in score columns must match length of true_vector.")
  }
  
  # Compute correlation with true vector
  cor_vals <- sapply(score_matrix, function(col) cor(col, true_vector, use = "pairwise.complete.obs"))
  
  # Order columns by absolute correlation (descending)
  ordered_idx <- order(abs(cor_vals), decreasing = TRUE)
  ordered_cols <- score_cols[ordered_idx]
  
  # Extract method name (e.g., "FABIA" from "loading_FABIA1")
  method_name <- sub(".*score_([A-Za-z]+).*", "\\1", ordered_cols[1])
  
  # Rename columns based on match order
  new_names <- paste0("score_", method_name, seq_along(ordered_cols))
  
  # Replace names in original df
  df_renamed <- df %>%
    dplyr::select(-all_of(score_cols)) %>%
    dplyr::bind_cols(setNames(score_matrix[, ordered_idx], new_names))
  
  return(df_renamed)
}


align_score_columns_by_truth_overall <- function(scores_df, true_df, method_name) {
  # Merge by feature
  df_merged <- merge(true_df, scores_df, by = "sample")
  
  # Pull true values
  true_vector <- df_merged$true_score_F1
  
  # Extract only numeric score columns for this method
  score_cols <- grep(paste0("^score_", method_name, "[0-9]+$"), names(df_merged), value = TRUE)
  
  if (length(score_cols) < 2) {
    stop(paste("Expected at least 2 score columns for method", method_name))
  }
  
  # Ensure score columns are numeric
  score_matrix <- df_merged#[, loading_cols, drop = FALSE]
  for (col in score_cols) {
    if (!is.numeric(score_matrix[[col]])) {
      score_matrix[[col]] <- as.numeric(as.character(score_matrix[[col]]))
    }
  }
  
  # Compute correlations
  cor_vals <- sapply(score_matrix[, score_cols], function(col) cor(col, true_vector, use = "pairwise.complete.obs"))
  
  # Order by absolute correlation
  ordered_idx <- order(abs(cor_vals), decreasing = TRUE)
  ordered_cols <- score_cols[ordered_idx]
  
  new_names <- paste0("matched_", method_name, seq_along(ordered_cols))
  old_names <- character(length(ordered_idx))
  for (j in seq_along(ordered_idx)) {
    old_names[j] <- paste0("score_", method_name, ordered_idx[j])
  }
  
  # Rename the columns in df_merged directly
  names(df_merged)[match(old_names, names(df_merged))] <- new_names
  
  return(df_merged)
}

############################################
# -------- END OF FUNCTIONS -------------- 
############################################

######################################
# -------- SIMULATION -------------- 
######################################

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
  global_env$n_samples <- n_samples
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
      simulated_data[[paste0("iteration_", i)]][["list_deltas"]][["delta2"]] <- rnorm(n_features_two, 0, 0.05)
      true_load_omic2_FACTOR2 <- simulated_data[[paste0("iteration_", i)]][["list_deltas"]][["delta2"]]
      
      indices_features.OMIC1.A <- simulated_data[[paste0("iteration_", i)]]$indices_features.1[[1]]
      indices_features.OMIC1.B <- simulated_data[[paste0("iteration_", i)]]$indices_features.1[[2]]
      indices_features.OMIC2.A <- simulated_data[[paste0("iteration_", i)]]$indices_features.2[[1]]
      indices_samples.1A <- simulated_data[[paste0("iteration_", i)]]$indices_samples[[1]]
      indices_samples.2B <- simulated_data[[paste0("iteration_", i)]]$indices_samples[[2]]
      
      # Create Features
      features_omic1 <- paste0("omic1_feature_", seq_len(length(true_load_omic1_FACTOR1)))
      features_omic2 <- paste0("omic2_feature_", seq_len(length(true_load_omic2_FACTOR1)))
      
      global_env$features_omic1 <- features_omic1
      global_env$features_omic2 <- features_omic2
      global_env$indices_features.OMIC1.A <- indices_features.OMIC1.A
      global_env$indices_features.OMIC1.B <- indices_features.OMIC1.B
      global_env$indices_features.OMIC2.A <- indices_features.OMIC2.A
      global_env$indices_samples.1A <- indices_samples.1A
      global_env$indices_samples.2B <- indices_samples.2B
      global_env$true_load_omic1_FACTOR1 <- true_load_omic1_FACTOR1
      global_env$true_load_omic1_FACTOR2 <- true_load_omic1_FACTOR2
      global_env$true_load_omic2_FACTOR1 <- true_load_omic2_FACTOR1
      global_env$true_load_omic2_FACTOR2 <- true_load_omic2_FACTOR2
      
      # Construct factor 1 loading data frame
      df_factor_loading <- data.frame(
        feature = c(features_omic1, features_omic2),
        #beta1   = c(true_load_omic1_FACTOR1, rep(NA, length(true_load_omic2_FACTOR1))),
        #delta1  = c(rep(NA, length(true_load_omic1_FACTOR1)), true_load_omic2_FACTOR1),
        true_loading_F1 = c(true_load_omic1_FACTOR1, true_load_omic2_FACTOR1),
        true_loading_F2 = c(true_load_omic1_FACTOR2, true_load_omic2_FACTOR2)
      )
      
      # TRUE SCORES
      samples <- paste0("sample_", seq_len(length(true_score_FACTOR1)))
      
      # Construct factor 1 loading data frame
      df_factor_score <- data.frame(
        sample = samples,
        true_score_F1 = c(true_score_FACTOR1, rep(NA, length(true_score_FACTOR1))),
        true_score_F2  = c(true_score_FACTOR2, rep(NA, length(true_score_FACTOR2)))
      )
      df_factor_score <- na.omit(df_factor_score)
      
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
      mfa_result <- func_mfa(current_data, num_biclusters, simulated_data, i)
      #plot(mfa_result[["scores"]][["factor_scores"]][["score_MFA1"]])
      #plot(mfa_result[["scores"]][["factor_scores"]][["score_MFA2"]])
      print(paste("GFA MODELLING IN PROGRESS..."))
      gfa_result <- func_gfa(current_data, num_biclusters, simulated_data, i)
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
      fabia_weights <- extract_id_from_feature(fabia_result[["weights"]][["all_weights"]])
      mofa_weights  <- extract_id_from_feature(mofa_result[["weights"]][["all_weights"]])
      mfa_weights  <- extract_id_from_feature(mfa_result[["weights"]][["original_weights"]])
      gfa_weights  <- extract_id_from_feature(gfa_result[["weights"]][["original_weights"]])
      
      
      out_fabia <- align_loading_columns_by_truth_overall(fabia_weights, df_factor_loading, "FABIA")
      out_mofa  <- align_loading_columns_by_truth_overall(mofa_weights,  df_factor_loading, "MOFA")
      out_gfa   <- align_loading_columns_by_truth_overall(gfa_weights,   df_factor_loading, "GFA")
      out_mfa   <- align_loading_columns_by_truth_overall(mfa_weights,   df_factor_loading, "MFA")
      
      # Step 1: Select only feature + matched columns from each output
      fabia_matched <- out_fabia[, c("feature", grep("^matched_FABIA", names(out_fabia), value = TRUE))]
      mofa_matched  <- out_mofa[,  c("feature", grep("^matched_MOFA",  names(out_mofa),  value = TRUE))]
      mfa_matched   <- out_mfa[,   c("feature", grep("^matched_MFA",   names(out_mfa),   value = TRUE))]
      gfa_matched   <- out_gfa[,   c("feature", grep("^matched_GFA",   names(out_gfa),   value = TRUE))]
      
      # Step 2: Merge all by "feature"
      df_merged_loading_all <- Reduce(function(x, y) merge(x, y, by = "feature", all = TRUE),
                                      list(fabia_matched, mofa_matched, mfa_matched, gfa_matched))
      
      df_merged_loading_all <- merge(df_merged_loading_all, fabia_weights, by = "feature")
      
      df_merged_loading_OMIC1 <- df_merged_loading_all %>%
        subset(dataview  == "omic.one")
      df_merged_loading_OMIC2 <- df_merged_loading_all %>%
        subset(dataview  == "omic.two")
      
      # Extract factor scores
      fabia_scores <- extract_id_from_sample(fabia_result[["scores"]][["scores"]])
      mofa_scores  <- extract_id_from_sample(mofa_result[["scores"]][["scores"]])
      mfa_scores  <- extract_id_from_sample(mfa_result[["scores"]][["original_scores"]])
      gfa_scores  <- extract_id_from_sample(gfa_result[["scores"]][["original_scores"]])
      
      
      out_score_fabia <- align_score_columns_by_truth_overall(fabia_scores, df_factor_score, "FABIA")
      out_score_mofa  <- align_score_columns_by_truth_overall(mofa_scores,  df_factor_score, "MOFA")
      out_score_gfa   <- align_score_columns_by_truth_overall(gfa_scores,   df_factor_score, "GFA")
      out_score_mfa   <- align_score_columns_by_truth_overall(mfa_scores,   df_factor_score, "MFA")
      
      
      # Step 1: Select only sample + matched columns from each output
      fabia_score_matched <- out_score_fabia[, c("sample", grep("^matched_FABIA", names(out_score_fabia), value = TRUE))]
      mofa_score_matched  <- out_score_mofa[,  c("sample", grep("^matched_MOFA",  names(out_score_mofa),  value = TRUE))]
      mfa_score_matched   <- out_score_mfa[,   c("sample", grep("^matched_MFA",   names(out_score_mfa),   value = TRUE))]
      gfa_score_matched   <- out_score_gfa[,   c("sample", grep("^matched_GFA",   names(out_score_gfa),   value = TRUE))]
      
      # Step 2: Merge all by "sample"
      df_merged_score_all <- Reduce(function(x, y) merge(x, y, by = "sample", all = TRUE),
                                    list(fabia_score_matched, mofa_score_matched, mfa_score_matched, gfa_score_matched))
      
      # OMIC ONE
      # 
      # fabia_df <- extract_id_from_feature(fabia_result[["weights"]][["omic.one_weights"]])
      # fabia_aligned <- align_loading_columns_by_truth(fabia_df, true_load_omic1_FACTOR1)
      # mofa_df  <- extract_id_from_feature(mofa_result[["weights"]][["omic.one_weights"]])
      # mofa_aligned <- align_loading_columns_by_truth(mofa_df, true_load_omic1_FACTOR1)
      # mfa_df   <- extract_id_from_feature(mfa_result[["weights"]][["omic.one_weights"]])
      # mfa_aligned <- align_loading_columns_by_truth(mfa_df, true_load_omic1_FACTOR1)
      gfa_df_omic1   <- extract_id_from_feature(gfa_result[["weights"]][["omic.one_weights"]])
      gfa_df_omic2   <- extract_id_from_feature(gfa_result[["weights"]][["omic.two_weights"]])
      
      #gfa_aligned <- align_loading_columns_by_truth(gfa_df, true_load_omic1_FACTOR1)
      # 
      # fabia_sub <- fabia_aligned[, c("ID", "loading_FABIA1")]
      # mofa_sub  <- mofa_aligned[, c("ID", "loading_MOFA1", "feature", "dataview")]
      # mfa_sub   <- mfa_df[, c("ID", "loading_MFA1")]
      gfa_omic1_metadata   <- gfa_df_omic1[, c("feature","ID", "signal_a", "signal_b")]
      gfa_omic2_metadata   <- gfa_df_omic2[, c("feature","ID", "signal_a")]
      # 
      # merged_loadings_omic1_F1 <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE),
      #                           list(fabia_sub, mofa_sub, mfa_sub, gfa_sub))
      
      #Convert all numeric columns to absolute values
      abs.merged_loadings_omic1 <- df_merged_loading_OMIC1
      abs.merged_loadings_omic1 <- abs.merged_loadings_omic1 %>%
        mutate(across(where(is.numeric), abs))
      
      
      # abs.merged_loadings_omic1_F1$loading_FABIA1= abs(abs.merged_loadings_omic1_F1$loading_FABIA1)
      # abs.merged_loadings_omic1_F1$loading_MOFA1 = abs(abs.merged_loadings_omic1_F1$loading_MOFA1)
      # abs.merged_loadings_omic1_F1$loading_MFA1 = abs(abs.merged_loadings_omic1_F1$loading_MFA1)
      # abs.merged_loadings_omic1_F1$loading_GFA1 = abs(abs.merged_loadings_omic1_F1$loading_GFA1)
      
      # fabia_sub <- fabia_aligned[, c("ID", "loading_FABIA2")]
      # mofa_sub  <- mofa_aligned[, c("ID", "loading_MOFA2", "feature", "dataview")]
      # mfa_sub   <- mfa_df[, c("ID", "loading_MFA2")]
      # gfa_sub   <- gfa_df[, c("ID", "loading_GFA2", "signal_a", "signal_b")]
      # 
      # merged_loadings_omic1_F2 <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE),
      #                                 list(fabia_sub, mofa_sub, mfa_sub, gfa_sub))
      # 
      # #Convert all numeric columns to absolute values
      # abs.merged_loadings_omic1_F2 <- merged_loadings_omic1_F2 
      # abs.merged_loadings_omic1_F2$loading_FABIA2= abs(abs.merged_loadings_omic1_F2$loading_FABIA2)
      # abs.merged_loadings_omic1_F2$loading_MOFA2 = abs(abs.merged_loadings_omic1_F2$loading_MOFA2)
      # abs.merged_loadings_omic1_F2$loading_MFA2 = abs(abs.merged_loadings_omic1_F2$loading_MFA2)
      # abs.merged_loadings_omic1_F2$loading_GFA2 = abs(abs.merged_loadings_omic1_F2$loading_GFA2)
      
      # OMIC TWO
      abs.merged_loadings_omic2 <- df_merged_loading_OMIC2
      abs.merged_loadings_omic2 <- abs.merged_loadings_omic2 %>%
        mutate(across(where(is.numeric), abs))
      
      # fabia_df_omic2 <- extract_id_from_feature(fabia_result[["weights"]][["omic.two_weights"]])
      # fabia_aligned <- align_loading_columns_by_truth(fabia_df_omic2, true_load_omic2_FACTOR1)
      # mofa_df_omic2  <- extract_id_from_feature(mofa_result[["weights"]][["omic.two_weights"]])
      # mofa_aligned <- align_loading_columns_by_truth(mofa_df_omic2, true_load_omic2_FACTOR1)
      # mfa_df_omic2   <- extract_id_from_feature(mfa_result[["weights"]][["omic.two_weights"]])
      # mfa_aligned <- align_loading_columns_by_truth(mfa_df_omic2, true_load_omic2_FACTOR1)
      # gfa_df_omic2   <- extract_id_from_feature(gfa_result[["weights"]][["omic.two_weights"]])
      # gfa_aligned <- align_loading_columns_by_truth(gfa_df_omic2, true_load_omic2_FACTOR1)
      # 
      # fabia_sub_omic2 <- fabia_aligned[, c("ID", "loading_FABIA1")]
      # mofa_sub_omic2  <- mofa_aligned[, c("ID", "loading_MOFA1", "feature", "dataview")]
      # mfa_sub_omic2   <- mfa_df_omic2[, c("ID", "loading_MFA1")]
      # gfa_sub_omic2   <- gfa_df_omic2[, c("ID", "loading_GFA1", "signal_a")]
      # 
      # merged_loadings_omic2_F1 <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE),
      #                                 list(fabia_sub_omic2, mofa_sub_omic2, mfa_sub_omic2, gfa_sub_omic2))
      # 
      # #Convert all numeric columns to absolute values
      # abs.merged_loadings_omic2_F1 <- merged_loadings_omic2_F1 
      # abs.merged_loadings_omic2_F1$loading_FABIA1= abs(abs.merged_loadings_omic2_F1$loading_FABIA1)
      # abs.merged_loadings_omic2_F1$loading_MOFA1 = abs(abs.merged_loadings_omic2_F1$loading_MOFA1)
      # abs.merged_loadings_omic2_F1$loading_MFA1 = abs(abs.merged_loadings_omic2_F1$loading_MFA1)
      # abs.merged_loadings_omic2_F1$loading_GFA1 = abs(abs.merged_loadings_omic2_F1$loading_GFA1)
      # 
      # fabia_sub_omic2 <- fabia_aligned[, c("ID", "loading_FABIA2")]
      # mofa_sub_omic2  <- mofa_aligned[, c("ID", "loading_MOFA2", "feature", "dataview")]
      # mfa_sub_omic2   <- mfa_df_omic2[, c("ID", "loading_MFA2")]
      # gfa_sub_omic2   <- gfa_df_omic2[, c("ID", "loading_GFA2")]
      # 
      # merged_loadings_omic2_F2 <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE),
      #                                    list(fabia_sub_omic2, mofa_sub_omic2, mfa_sub_omic2, gfa_sub_omic2))
      
      #Convert all numeric columns to absolute values
      # abs.merged_loadings_omic2_F2 <- merged_loadings_omic2_F2 
      # abs.merged_loadings_omic2_F2$loading_FABIA2= abs(abs.merged_loadings_omic2_F2$loading_FABIA2)
      # abs.merged_loadings_omic2_F2$loading_MOFA2 = abs(abs.merged_loadings_omic2_F2$loading_MOFA2)
      # abs.merged_loadings_omic2_F2$loading_MFA2 = abs(abs.merged_loadings_omic2_F2$loading_MFA2)
      # abs.merged_loadings_omic2_F2$loading_GFA2 = abs(abs.merged_loadings_omic2_F2$loading_GFA2)
      
      # SCORES
      # FACTOR 1
      df_merged_scores <- df_merged_score_all
      df_merged_scores <- df_merged_scores %>%
        mutate(across(where(is.numeric), abs))
      
      # fabia_df <- extract_id_from_sample (fabia_result[["scores"]][["scores"]])
      # fabia_aligned <- align_scores_columns_by_truth(fabia_df, true_score_FACTOR1)
      # mofa_df  <- extract_id_from_sample (mofa_result[["scores"]][["scores"]])
      # mofa_aligned <- align_scores_columns_by_truth(mofa_df, true_score_FACTOR1)
      # mfa_df   <- extract_id_from_sample (mfa_result[["scores"]][["factor_scores"]])
      # mfa_aligned <- align_scores_columns_by_truth(mfa_df, true_score_FACTOR1)
      gfa_df_score   <- extract_id_from_sample (gfa_result[["scores"]][["factor_scores"]])
      # gfa_aligned <- align_scores_columns_by_truth(gfa_df, true_score_FACTOR1)
      # 
      # fabia_sub <- fabia_aligned[, c("ID", "score_FABIA1")]
      # mofa_sub  <- mofa_aligned[, c("ID", "score_MOFA1", "sample")]
      # mfa_sub   <- mfa_df[, c("ID", "score_MFA1")]
      gfa_sample_metadata   <- gfa_df_score[, c("sample", "ID", "signal_a", "signal_b")]
      # 
      # merged_scores_F1 <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE),
      #                                    list(fabia_sub, mofa_sub, mfa_sub, gfa_sub))
      # 
      # #Convert all numeric columns to absolute values
      # abs.merged_scores_F1 <- merged_scores_F1 
      # abs.merged_scores_F1$score_FABIA1= abs(abs.merged_scores_F1$score_FABIA1)
      # abs.merged_scores_F1$score_MOFA1 = abs(abs.merged_scores_F1$score_MOFA1)
      # abs.merged_scores_F1$score_MFA1 = abs(abs.merged_scores_F1$score_MFA1)
      # abs.merged_scores_F1$score_GFA1 = abs(abs.merged_scores_F1$score_GFA1)
      # 
      # # FACTOR 2
      # fabia_sub <- fabia_aligned[, c("ID", "score_FABIA2")]
      # mofa_sub  <- mofa_aligned[, c("ID", "score_MOFA2", "sample")]
      # mfa_sub   <- mfa_df[, c("ID", "score_MFA2")]
      # gfa_sub   <- gfa_df[, c("ID", "score_GFA2", "signal_a", "signal_b")]
      # 
      # merged_scores_F2 <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE),
      #                            list(fabia_sub, mofa_sub, mfa_sub, gfa_sub))
      # 
      # #Convert all numeric columns to absolute values
      # abs.merged_scores_F2 <- merged_scores_F2 
      # abs.merged_scores_F2$score_FABIA2= abs(abs.merged_scores_F2$score_FABIA2)
      # abs.merged_scores_F2$score_MOFA2 = abs(abs.merged_scores_F2$score_MOFA2)
      # abs.merged_scores_F2$score_MFA2 = abs(abs.merged_scores_F2$score_MFA2)
      # abs.merged_scores_F2$score_GFA2 = abs(abs.merged_scores_F2$score_GFA2)
      
      # Merge loadings
      # Manually remove overlapping columns *before* joining
      # abs.merged_loadings_omic1_F2_clean <- abs.merged_loadings_omic1_F2 %>%
      #   select(-any_of(names(abs.merged_loadings_omic1_F1)[names(abs.merged_loadings_omic1_F1) != "feature"]))
      
      # Join with no conflicts
      FeatureData.Omic1 <- left_join(abs.merged_loadings_omic1, gfa_omic1_metadata, by = "feature")
      FeatureData.Omic2 <- left_join(abs.merged_loadings_omic2, gfa_omic2_metadata, by = "feature")
      FeatureData.Omic2$signal_b = FALSE
      SampleData <- left_join(df_merged_scores, gfa_sample_metadata, by = "sample")
      
      FeatureData.Omic1 <- FeatureData.Omic1 %>%
        select(-contains("loading"))
      
      FeatureData.Omic2 <- FeatureData.Omic2 %>%
        select(-contains("loading"))
      
      SampleData <- SampleData %>%
        select(-contains("score"))
      
      # Replace "matched" with "loading" in feature-level data
      names(FeatureData.Omic1) <- gsub("matched", "loading", names(FeatureData.Omic1))
      names(FeatureData.Omic2) <- gsub("matched", "loading", names(FeatureData.Omic2))
      
      # Replace "matched" with "score" in sample-level data
      names(SampleData) <- gsub("matched", "score", names(SampleData))
      
      # abs.merged_loadings_omic2_F2_clean <- abs.merged_loadings_omic2_F2 %>%
      #   select(-any_of(names(abs.merged_loadings_omic2_F1)[names(abs.merged_loadings_omic2_F1) != "feature"]))
      
      # Join with no conflicts
      # FeatureData.Omic2 <- abs.merged_loadings_omic1 #left_join(abs.merged_loadings_omic2_F1, abs.merged_loadings_omic2_F2_clean, by = "feature")
      # FeatureData.Omic2$signal_b = FALSE
      
      # Scores
      # abs.merged_scores_F2_clean <- abs.merged_scores_F2 %>%
      #   select(-any_of(names(abs.merged_scores_F1)[names(abs.merged_scores_F1) != "sample"]))
      # 
      # # Join with no conflicts
      # SampleData <- left_join(abs.merged_scores_F1, abs.merged_scores_F2_clean, by = "sample")
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
                abs(loading_FABIA1) >= varphi.five(
                  loading_or_score = loading_FABIA1,
                  omic = 1, factor = 1, type = "loading"
                ), 1, 0),
              
              mofa_index_a = ifelse(
                abs(loading_MOFA1) >= varphi.five(
                  loading_or_score = loading_MOFA1,
                  omic = 1, factor = 1, type = "loading"
                ), 1, 0),
              
              mfa_index_a = ifelse(
                abs(loading_MFA1) >= varphi.five(
                  loading_or_score = loading_MFA1,
                  omic = 1, factor = 1, type = "loading"
                ), 1, 0),
              
              gfa_index_a = ifelse(
                abs(loading_GFA1) >= varphi.five(
                  loading_or_score = loading_GFA1,
                  omic = 1, factor = 1, type = "loading"
                ), 1, 0),
              
              signal_index_a = ifelse(signal_a == 'TRUE', 1, 0),
              
              fabia_index_b = ifelse(
                abs(loading_FABIA2) >= varphi.five(
                  loading_or_score = loading_FABIA2,
                  omic = 1, factor = 2, type = "loading"
                ), 1, 0),
              
              mofa_index_b = ifelse(
                abs(loading_MOFA2) >= varphi.five(
                  loading_or_score = loading_MOFA2,
                  omic = 1, factor = 2, type = "loading"
                ), 1, 0),
              
              mfa_index_b = ifelse(
                abs(loading_MFA2) >= varphi.five(
                  loading_or_score = loading_MFA2,
                  omic = 1, factor = 2, type = "loading"
                ), 1, 0),
              
              gfa_index_b = ifelse(
                abs(loading_GFA2) >= varphi.five(
                  loading_or_score = loading_GFA2,
                  omic = 1, factor = 2, type = "loading"
                ), 1, 0),
              
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
        
        FeatureData.a.one <- FeatureData.a.one  %>% arrange(as.numeric(ID.x))
        
        ## ----------------- FeatureData.b.two ----------------
        if (varphi_function == "varphi.five"){
          # Apply the function to Data_Ft.a
          FeatureData.b.two <- FeatureData.Omic2 %>%
            mutate(
              fabia_index_a = ifelse(
                abs(loading_FABIA1) >= do.call(get(varphi_function), list(
                  loading_or_score = loading_FABIA1,
                  omic = 2, factor = 1, type = "loading"
                )), 1, 0),
              
              mofa_index_a = ifelse(
                abs(loading_MOFA1) >= do.call(get(varphi_function), list(
                  loading_or_score = loading_MOFA1,
                  omic = 2, factor = 1, type = "loading"
                )), 1, 0),
              
              mfa_index_a = ifelse(
                abs(loading_MFA1) >= do.call(get(varphi_function), list(
                  loading_or_score = loading_MFA1,
                  omic = 2, factor = 1, type = "loading"
                )), 1, 0),
              
              gfa_index_a = ifelse(
                abs(loading_GFA1) >= do.call(get(varphi_function), list(
                  loading_or_score = loading_GFA1,
                  omic = 2, factor = 1, type = "loading"
                )), 1, 0),
              
              signal_index_a = ifelse(signal_a == 'TRUE', 1, 0),
              
              ##################################################################
              # Since no signal in the factor 2 omic2, I want to restrict this #
              # to zero so that no signal is identified here; therefore I use  #
              # greater than maximum knowing no value will be greater than max.#
              ##################################################################
              fabia_index_b = ifelse(
                abs(loading_FABIA2) > do.call(get(varphi_function), list(
                  loading_or_score = loading_FABIA2,
                  omic = 2, factor = 2, type = "loading"
                )), 1, 0),
              
              mofa_index_b = ifelse(
                abs(loading_MOFA2) > do.call(get(varphi_function), list(
                  loading_or_score = loading_MOFA2,
                  omic = 2, factor = 2, type = "loading"
                )), 1, 0),
              
              mfa_index_b = ifelse(
                abs(loading_MFA2) > do.call(get(varphi_function), list(
                  loading_or_score = loading_MFA2,
                  omic = 2, factor = 2, type = "loading"
                )), 1, 0),
              
              gfa_index_b = ifelse(
                abs(loading_GFA2) > do.call(get(varphi_function), list(
                  loading_or_score = loading_GFA2,
                  omic = 2, factor = 2, type = "loading"
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
        
        FeatureData.b.two <- FeatureData.b.two %>% arrange(as.numeric(ID.x))
        
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
        
        # method_output_corrected <- list(
        #   merged_scores_F1 = merged_scores_F1,
        #   merged_scores_F2 = merged_scores_F2,
        #   merged_loadings_omic1_F1 = merged_loadings_omic1_F1,
        #   merged_loadings_omic1_F2 = merged_loadings_omic1_F2,
        #   merged_loadings_omic2_F1 = merged_loadings_omic2_F1,
        #   merged_loadings_omic2_F2 = merged_loadings_omic2_F2
        # ) 
        
        # dataset <- sprintf("variance_%d_iteration_%d", sigma, i)
        # dataset_output_corrected[[dataset]] <- method_output_corrected
        
      } # closes varphi
    } # closes iterations
  } # closes sigma
  # ---------------- simulation results ------------------------
  # Combine results into a list
  results <- list(
    jaccard_results = jaccard_results,
    jaccard_comparison_results = jaccard_comparison_results,
    per_measures_results = per_measures_results,#,
    dataset_output = dataset_output#,
    #dataset_output_corrected = dataset_output_corrected
  )
  return(results)
  
} # ends function

## -------------------------- Simulate multifactor -------------------------- 

#setting_seed <- sample(1e6, 1)
set.seed(54321)

sigmas <- 29  # Define sigma or the range of sigmas

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
    num_iteration = 100) 
  
  # Save the result to a file
  saveRDS(
    simulation_result,
    file = paste0("/user/leuven/364/vsc36498/", result_name, ".rds")
    #file = paste0("C:/Users/Lenovo/Downloads/", result_name, ".rds")
    #file = paste0("C:/Users/bosangir/Downloads/", result_name, ".rds")
  )                                  
  #}
  
  # Print progress message
  cat("Completed simulations for sigma =", sigma, "\n")
} 
#gfa_seed
# Print completion message
cat("All simulations completed and saved!\n")

############################################
# -------- END OF SIMULATION -------------- 
############################################
