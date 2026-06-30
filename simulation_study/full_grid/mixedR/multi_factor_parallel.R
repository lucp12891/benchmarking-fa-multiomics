setwd("/user/leuven/364/vsc36498")
set.seed(645) #644

# Load packages and suppress message
packages <- c("ggplot2", "dplyr", "tidyr", "viridis", "fabia", "FactoMineR", "readr", "readxl", "tidyverse",
                      "stringr", "basilisk", "MOFA2", "data.table", "GFA", "pROC", "mclust", "zoo")

for (pkg in packages) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

#  --------------------------------- MultiFactor: Simulation function ----------------------------
multiple_factor <- function(n_features_one, n_features_two, n_samples, sigmas, iterations, n_factors){
    set.seed(694)
    # # FABIA check function
    # func_fabia <- function(data, num = n_factors) {
    #   #set.seed(123)
    #   fabia_object <- fabia(as.matrix(data), p = n_factors, alpha = 0.01, cyc = 1000, spl = 0.5, spz = 0.5, 
    #                         random = 1.0, center = 2, norm = 2, lap = 1.0, nL = 1)
    #   fab_loading <- fab_loading(fabia_object, num)
    #   fab_score <- fab_score(fabia_object, num)
    #   X <- fabia_object@X
    #   fabia_result <- list(fab_loading = fab_loading, fab_score = fab_score, X = X)
    #   return(fabia_result)
    # }

  # Initialize lists to store omic data for all iterations and sigmas
  all_omic_data <- list()

  for (iter in 1:iterations) {  # Iterate through the specified number of iterations
    omic.one <- list()
    omic.two <- list()

    for (k in 1:length(sigmas)) {  # Iterate through the list of sigmas
      n_s <- n_samples
      valid_data <- FALSE  # Flag to check if data passes FABIA condition

      #while (!valid_data) {
        # Generate first omic data

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
        # omic.one <- list()
        # omic.one[[length(omic.one) + 1]] <- matrix(data.1, n_samples, n_features_one) + matrix(eps1, n_samples, n_features_one) # signal + noise
        # colnames(omic.one[[length(omic.one)]]) <- paste0('omic1_feature_', seq_len(n_features_one))
        # rownames(omic.one[[length(omic.one)]]) <- paste0('sample_', seq_len(n_samples))
        omic1_data <- matrix(data.1, n_samples, n_features_one) + matrix(eps1, n_samples, n_features_one) # signal + noise
        colnames(omic1_data) <- paste0('omic1_feature_', seq_len(n_features_one))
        rownames(omic1_data) <- paste0('sample_', seq_len(n_samples))
        
        #dataset <- omic.one[[1]]
        #image(c(1:dim(dataset)[1]), c(1:dim(dataset)[2]), dataset, ylab = "Features", xlab = "Samples")
        
        # Generate random gamma values based on the max_factors
        # Replace 'alpha' with 'gamma' in the names of list_gammas
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
        # omic.two <- list()
        # omic.two[[length(omic.two) + 1]] <- matrix(data.2, n_samples, n_features_two) + matrix(eps2, n_samples, n_features_two) # signal + noise
        # colnames(omic.two[[length(omic.two)]]) <- paste0('omic2_feature_', seq_len(n_features_two))
        # rownames(omic.two[[length(omic.two)]]) <- paste0('sample_', seq_len(n_samples))
        omic2_data <- matrix(data.2, n_samples, n_features_two) + matrix(eps2, n_samples, n_features_two) # signal + noise
        colnames(omic2_data) <- paste0('omic2_feature_', seq_len(n_features_two))
        rownames(omic2_data) <- paste0('sample_', seq_len(n_samples))
        
        
        #dataset <- data.2#dataset <- omic.two[[1]]
        #image(c(1:dim(dataset)[1]), c(1:dim(dataset)[2]), dataset, ylab = "Features", xlab = "Samples")
        
        # Concatenate datasets
        concatenated_data <- cbind(omic2_data, omic1_data)
        
        # # Run FABIA and check conditions
        # fabia_result <- func_fabia(concatenated_data, num = 2)
        # fab_loading1 <- fabia_result$fab_loading$loading_FABIA1
        # fab_score1 <- fabia_result$fab_score$score_FABIA1
        # fab_loading2 <- fabia_result$fab_loading$loading_FABIA2
        # fab_score2 <- fabia_result$fab_score$score_FABIA2
        # 
        # Validate FABIA results
        # if (!any(is.na(fab_loading1)) && mean(fab_loading1, na.rm = TRUE) != 0 &&
        #     !any(is.na(fab_score1)) && mean(fab_score1, na.rm = TRUE) != 0 &&
        #     !any(is.na(fab_loading2)) && mean(fab_loading2, na.rm = TRUE) != 0 &&
        #     !any(is.na(fab_score2)) && mean(fab_score2, na.rm = TRUE) != 0) {
        #   valid_data <- TRUE
        # } else {
        #   message("FABIA validation failed, regenerating concatenated_data...")
        # }
      #}
      
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
      feature_sig2_end = f_sig2_e1
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
    
    # List container for feature weights and factor score results
    #fab_loading <- list()
    #fab_score <- list()
    
    # Extract feature weights and factor scores
    # Iterate over the number of biclusters
    fabia_loading <- fab_loading(fabia_object, BC_num)
    fabia_score <- fab_score(fabia_object, BC_num)
    # 
    # # Validate FABIA results
    # if (!any(is.na(fab_loading$loading_FABIA1)) && mean(fab_loading$loading_FABIA1, na.rm = TRUE) != 0 &&
    #     !any(is.na(fab_score$score_FABIA1)) && mean(fab_score$score_FABIA1, na.rm = TRUE) != 0 && 
    #     !any(is.na(fab_loading$loading_FABIA2)) && mean(fab_loading$loading_FABIA2, na.rm = TRUE) != 0 &&
    #     !any(is.na(fab_score$score_FABIA2)) && mean(fab_score$score_FABIA2, na.rm = TRUE) != 0) {
    #   valid_data <- TRUE  # Data passes validation
    # } else {
    #   message("FABIA validation failed, re-run FABIA...")
    # }
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
    #select(sample, score_FABIA1)
  
  # scores_factor_two <- result_score %>%
  #   mutate(
  #     ID = ifelse(grepl("sample_", sample),
  #                 as.numeric(sub(".*sample_(\\d+)", "\\1", sample)), NA)
  #   ) %>%
  #   arrange(ID) %>%
  #   select(sample, score_FABIA2)
  
  scores_composite <- list(
    #factor_one_scores = scores_factor_one,
    #factor_two_scores = scores_factor_two,
    scores = result_score
  )
  
  fabia_result <- list(weights = weights_composite, scores = scores_composite)
  
  # Return results
  #fabia_result <- list(fab_loading = fab_loading, fab_score = fab_score, X = X)
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

#fab_score <- fab_score(fabia_object, BC_num)
#fab_loading <- fab_loading(fabia_object, BC_num)
#fabia_result <- func_fabia(data, BC_num)
# --------------------------------------- MOFA+ Function -----------------------------------------

## ----------------------------------- 2. MOFA+ Main Function ---------------------------------------
func_mofa <- function(data, num) {
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
      #select(sample, score_MOFA1)
    
    # scores_factor_two <- result_score %>%
    #   mutate(
    #     ID = ifelse(grepl("sample_", sample),
    #                 as.numeric(sub(".*sample_(\\d+)", "\\1", sample)), NA)
    #   ) %>%
    #   arrange(ID) %>%
    #   select(sample, score_MOFA2)
    
    scores_composite <- list(
      #factor_one_scores = scores_factor_one,
      #factor_two_scores = scores_factor_two,
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

# mofa_score = mofa_score(data, BC_num)
# mofa_loading = mofa_loading(data, BC_num)
# mofa_result = func_mofa(data, BC_num)

# --------------------------------------- MFA Function -----------------------------------------

## ----------------------------------- 3. MFA Main Function ---------------------------------------

func_mfa <- function(data, num) {
  
  #set.seed(123)
  BC_num = num
  load_data <- data # dataset
  
  # Access variables from the shared environment
  n_features_one <- 4000#global_env$n_features_one
  n_features_two <- 3000#global_env$n_features_two
  
  second_omic <- load_data[, 1:n_features_two]
  first_omic <- load_data[, (n_features_two + 1):(n_features_one + n_features_two)]
  cdata = cbind(first_omic, second_omic)
  
  mfa_data <- cdata
  
  # mfa analysis
  #mfa_data = data.frame(mfa_data)
  
  # Create a list specifying which columns are quantitative or qualitative
  
  mfa_object = MFA(as.matrix(mfa_data), 
                   group = c(n_features_one, n_features_two), 
                   type = c("s","s"), 
                   name.group = c("first.omic", "second.omic"),
                   graph = FALSE)
  
  #mfa_loading <- list()
  #mfa_score <- list()
  
  mfa_loading <- mfa_loading(mfa_object, BC_num)
  mfa_score <- mfa_score(mfa_object, BC_num)
  
  # Loading separation per factor
  
  result_load <- mfa_loading
  
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
  
  # Calculate the standard deviation of gfa_load$loading_GFA1 where signal_b == FALSE
  indices_features.1b = simulated_data[[paste0("iteration_", i)]]$indices_features.1[2]
  
  indices_features.1b_vector <- unlist(indices_features.1b)  # Convert to atomic vector
  indices_features.1b_sorted <- sort(indices_features.1b_vector)
  
  in_range.b <- random.data.a$ID %in% indices_features.1b_sorted
  
  simulated_features_b <- data.frame(feature = load_omic_one$feature, signal_b = in_range.b)
  mfa_omic_one_load <- merge(x=load_omic_one,y=simulated_features_b, by="feature", all = TRUE)
  
  var_false <- var(mfa_omic_one_load$loading_MFA1[mfa_omic_one_load$signal_b == FALSE], na.rm = TRUE)
  mean_false <- mean(mfa_omic_one_load$loading_MFA1[mfa_omic_one_load$signal_b == FALSE], na.rm = TRUE)
  
  # Suppress elevated factor
  mfa_omic_one_load$loading_MFA1 <- ifelse(
    mfa_omic_one_load$signal_b == TRUE,
    rnorm(length(mfa_omic_one_load$loading_MFA1[mfa_omic_one_load$signal_b == TRUE]), mean_false, var_false), # Action if TRUE
    mfa_omic_one_load$loading_MFA1# Action if FALSE
  )
  
  # FACTOR 2
  indices_features.1a = simulated_data[[paste0("iteration_", i)]]$indices_features.1[1]
  
  indices_features.1a_vector <- unlist(indices_features.1a)  # Convert to atomic vector
  indices_features.1a_sorted <- sort(indices_features.1a_vector)
  
  in_range.a <- random.data.a$ID %in% indices_features.1a_sorted
  
  simulated_features_a <- data.frame(feature = load_omic_one$feature, signal_a = in_range.a)
  mfa_omic_one_load <- merge(x=mfa_omic_one_load,y=simulated_features_a, by="feature", all = TRUE)
  
  var_false <- var(mfa_omic_one_load$loading_MFA2[mfa_omic_one_load$signal_a == FALSE], na.rm = TRUE)
  mean_false <- mean(mfa_omic_one_load$loading_MFA2[mfa_omic_one_load$signal_a == FALSE], na.rm = TRUE)
  
  # Suppress elevated factor
  mfa_omic_one_load$loading_MFA2 <- ifelse(
    mfa_omic_one_load$signal_a == TRUE,
    rnorm(length(mfa_omic_one_load$loading_MFA2[mfa_omic_one_load$signal_a == TRUE]), mean_false, var_false), # Action if TRUE
    mfa_omic_one_load$loading_MFA2# Action if FALSE
  )
  
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
  
  #weights_composite <- list(factor_one_weights = mfa_omic_one_load, factor_two_weights = load_omic_two, all_weights = result_load_all)
  
  # - - - - - - - - - - - - - - - - - - - MFA SCORES - - - - - - - - - - - - - - - - - #
  
  # Factor scores separation
  result_score <- mfa_score
  
  # Factor 1
  
  scores_factor_one <- result_score %>%
    mutate(
      ID = ifelse(grepl("sample_", sample),
                  as.numeric(sub(".*sample_(\\d+)", "\\1", sample)), NA)) %>%
    arrange(ID) %>%
    select(ID, sample, score_MFA1)  # Use `select()` to choose specific columns
  
  # Calculate the standard deviation of gfa_load$loading_GFA1 where signal_b == FALSE
  indices_samples.a = simulated_data[[paste0("iteration_", i)]]$indices_samples[1]
  indices_samples.b = simulated_data[[paste0("iteration_", i)]]$indices_samples[2]
  
  indices_samples.a_vector <- unlist(indices_samples.a)  # Convert to atomic vector
  indices_samples.a_sorted <- sort(indices_samples.a_vector)
  indices_samples.b_vector <- unlist(indices_samples.b)  # Convert to atomic vector
  indices_samples.b_sorted <- sort(indices_samples.b_vector)
  
  in_range_sample_a <- scores_factor_one$ID %in% indices_samples.a_sorted
  in_range_sample_b <- scores_factor_one$ID %in% indices_samples.b_sorted
  
  
  simulated_samples_a <- data.frame(sample = scores_factor_one$sample, signal_a = in_range_sample_a)
  simulated_samples_b <- data.frame(sample = scores_factor_one$sample, signal_b = in_range_sample_b)
  
  mfa_omic_one_score <- merge(merge(x=scores_factor_one,y=simulated_samples_a, by="sample", all = TRUE), simulated_samples_b,  by="sample", all = TRUE)
  
  var_false <- var(
    mfa_omic_one_score$score_MFA1[mfa_omic_one_score$signal_a == FALSE & mfa_omic_one_score$signal_b == FALSE], 
    na.rm = TRUE
  )
  
  mean_false <- mean(
    mfa_omic_one_score$score_MFA1[mfa_omic_one_score$signal_a == FALSE & mfa_omic_one_score$signal_b == FALSE], 
    na.rm = TRUE
  )
  
  # Suppress elevated factor
  mfa_omic_one_score$score_MFA1 <- ifelse(
    mfa_omic_one_score$signal_b == TRUE,
    rnorm(length(mfa_omic_one_score$score_MFA1[mfa_omic_one_score$signal_b == TRUE]), mean_false, var_false), # Action if TRUE
    mfa_omic_one_score$score_MFA1# Action if FALSE
  )
  
  # Factor 2
  
  scores_factor_two <- result_score %>%
    mutate(
      ID = ifelse(grepl("sample_", sample),
                  as.numeric(sub(".*sample_(\\d+)", "\\1", sample)), NA)) %>%
    arrange(ID) %>%
    select(ID, sample, score_MFA2)  # Use `select()` to choose specific columns
  
  mfa_omic_two_score <- merge(merge(x=scores_factor_two,y=simulated_samples_a, by="sample", all = TRUE),simulated_samples_b,  by="sample", all = TRUE)
  
  var_false <- var(
    mfa_omic_two_score$score_MFA2[mfa_omic_two_score$signal_a == FALSE & mfa_omic_two_score$signal_b == FALSE], 
    na.rm = TRUE
  )
  
  mean_false <- mean(
    mfa_omic_two_score$score_MFA2[mfa_omic_two_score$signal_a == FALSE & mfa_omic_two_score$signal_b == FALSE], 
    na.rm = TRUE
  )
  
  # Suppress elevated factor
  mfa_omic_two_score$score_MFA2 <- ifelse(
    mfa_omic_two_score$signal_a == TRUE,
    rnorm(length(mfa_omic_two_score$score_MFA2[mfa_omic_two_score$signal_a == TRUE]), mean_false, var_false), # Action if TRUE
    mfa_omic_two_score$score_MFA2# Action if FALSE
  )
  
  # Select only specific columns from Y
  mfa_omic_two_score_selected <- mfa_omic_two_score %>%
    select(sample, score_MFA2)  # Replace 'column1', 'column2' with desired column names
  
  # Perform the merge, keeping all columns from X
  factor_scores <- merge(
    mfa_omic_one_score, 
    mfa_omic_two_score_selected, 
    by = "sample", 
    all.x = TRUE
  )
  
  
  weights_composite <- list(omic.one_weights = mfa_omic_one_load, omic.two_weights = load_omic_two, original_weights = result_load_all)
  scores_composite <- list(factor_scores = factor_scores, original_scores = result_score)
  
  mfa_result <- list(weights = weights_composite, scores = scores_composite)
  
  return(mfa_result)
}

## ----------------------------------- 3.1 MFA Loading Function ---------------------------------------

mfa_loading <- function(mfa_object, BC_num) {
  # Check if BC_num is greater than 1
  if (BC_num > 1) {
    # Initialize an empty data frame for results
    loadings_mfa_df <- data.frame()
    
    # Extract loadings for each BC_num and combine into a single data frame
    for (i in 1:BC_num) {
      loading_column <- mfa_object$quanti.var$coord[, i]  # Extract the loading for the current factor
      column_name <- paste0("loading_MFA", i)  # Name for the column
      
      # Add the score column to the data frame
      if (ncol(loadings_mfa_df) == 0) {
        loadings_mfa_df <- data.frame(feature = rownames(mfa_object$quanti.var$coord))  # Initialize with feature names
      }
      loadings_mfa_df[[column_name]] <- loading_column
    }
    mfa_df_loading = loadings_mfa_df
    return(mfa_df_loading)
  } else {
    # Handle single BC_num case
    loading_MFA <- mfa_object$quanti.var$coord[, BC_num]
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
      score_column <- mfa_object$ind$coord[, i]  # Extract the score for the current factor
      column_name <- paste0("score_MFA", i)  # Name for the column
      
      # Add the score column to the data frame
      if (ncol(scores_mfa_df) == 0) {
        scores_mfa_df <- data.frame(sample = rownames(mfa_object$ind$coord))  # Initialize with sample names
      }
      scores_mfa_df[[column_name]] <- score_column
    }
    mfa_df_score = scores_mfa_df
    return(mfa_df_score)
  } else {
    # Handle single BC_num case
    score_MFA <- mfa_object$ind$coord[, BC_num]
    mfa_df_score <- as.data.frame(score_MFA)
    mfa_df_score$sample <- rownames(mfa_df_score)
    return(mfa_df_score)
  }
}

# mfa_score = mfa_score(mfa_object, BC_num = 2)
# mfa_loading = mfa_loading(mfa_object, BC_num = 2)
# mfa_res = func_mfa(data, num = 2)

# --------------------------------------- GFA Function -----------------------------------------

## ----------------------------------- 3. GFA Main Function ---------------------------------------
func_gfa <- function(data, num) {
  #set.seed(123)
  
  # Access variables from the shared environment
  n_features_one <- 4000#global_env$n_features_one
  n_features_two <- 3000#global_env$n_features_two
  
  second_omic <- data[, 1:n_features_two]
  first_omic <- data[, (n_features_two + 1):(n_features_one + n_features_two)]
  cdata = as.matrix(cbind(first_omic, second_omic))
  
  gfa_dt = as.data.frame(cdata)
  
  merged_GFA_data = list(t(gfa_dt))
  model_option <- getDefaultOpts()
  model_option$iter.max <- 1000
  model_option$iter.burnin <- 10
  gfa_object <- gfa(t(merged_GFA_data), K= num, opts=model_option)
  #gfa_loading <- list()
  #gfa_score <- list()
  
  gfa_loading <- gfa_loading(gfa_object, num)
  gfa_score <- gfa_score(gfa_object, num)
  
  # Loading separation per factor
  
  result_load <- gfa_loading
  
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
  
  # Calculate the standard deviation of gfa_load$loading_GFA1 where signal_b == FALSE
  indices_features.1b = simulated_data[[paste0("iteration_", i)]]$indices_features.1[2]
  
  indices_features.1b_vector <- unlist(indices_features.1b)  # Convert to atomic vector
  indices_features.1b_sorted <- sort(indices_features.1b_vector)
  
  in_range.b <- random.data.a$ID %in% indices_features.1b_sorted
  
  simulated_features_b <- data.frame(feature = load_omic_one$feature, signal_b = in_range.b)
  gfa_omic_one_load <- merge(x=load_omic_one,y=simulated_features_b, by="feature", all = TRUE)
  
  var_false <- var(gfa_omic_one_load$loading_GFA1[gfa_omic_one_load$signal_b == FALSE], na.rm = TRUE)
  mean_false <- mean(gfa_omic_one_load$loading_GFA1[gfa_omic_one_load$signal_b == FALSE], na.rm = TRUE)
  
  # Suppress elevated factor
  gfa_omic_one_load$loading_GFA1 <- ifelse(
    gfa_omic_one_load$signal_b == TRUE,
    rnorm(length(gfa_omic_one_load$loading_GFA1[mfa_omic_one_load$signal_b == TRUE]), mean_false, var_false), # Action if TRUE
    gfa_omic_one_load$loading_GFA1# Action if FALSE
  )
  
  # FACTOR 2
  indices_features.1a = simulated_data[[paste0("iteration_", i)]]$indices_features.1[1]
  
  indices_features.1a_vector <- unlist(indices_features.1a)  # Convert to atomic vector
  indices_features.1a_sorted <- sort(indices_features.1a_vector)
  
  in_range.a <- random.data.a$ID %in% indices_features.1a_sorted
  
  simulated_features_a <- data.frame(feature = load_omic_one$feature, signal_a = in_range.a)
  gfa_omic_one_load <- merge(x=gfa_omic_one_load,y=simulated_features_a, by="feature", all = TRUE)
  
  var_false <- var(mfa_omic_one_load$loading_MFA2[mfa_omic_one_load$signal_a == FALSE], na.rm = TRUE)
  mean_false <- mean(mfa_omic_one_load$loading_MFA2[mfa_omic_one_load$signal_a == FALSE], na.rm = TRUE)
  
  # Suppress elevated factor
  gfa_omic_one_load$loading_GFA2 <- ifelse(
    gfa_omic_one_load$signal_a == TRUE,
    rnorm(length(gfa_omic_one_load$loading_GFA2[gfa_omic_one_load$signal_a == TRUE]), mean_false, var_false), # Action if TRUE
    gfa_omic_one_load$loading_GFA2# Action if FALSE
  )
  
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
  
  #weights_composite <- list(factor_one_weights = mfa_omic_one_load, factor_two_weights = load_omic_two, all_weights = result_load_all)
  
  # - - - - - - - - - - - - - - - - - - - MFA SCORES - - - - - - - - - - - - - - - - - #
  
  # Factor scores separation
  result_score <- gfa_score
  
  # Factor 1
  
  scores_factor_one <- result_score %>%
    mutate(
      ID = ifelse(grepl("sample_", sample),
                  as.numeric(sub(".*sample_(\\d+)", "\\1", sample)), NA)) %>%
    arrange(ID) %>%
    select(ID, sample, score_GFA1)  # Use `select()` to choose specific columns
  
  # Calculate the standard deviation of gfa_load$loading_GFA1 where signal_b == FALSE
  indices_samples.a = simulated_data[[paste0("iteration_", i)]]$indices_samples[1]
  indices_samples.b = simulated_data[[paste0("iteration_", i)]]$indices_samples[2]
  
  indices_samples.a_vector <- unlist(indices_samples.a)  # Convert to atomic vector
  indices_samples.a_sorted <- sort(indices_samples.a_vector)
  indices_samples.b_vector <- unlist(indices_samples.b)  # Convert to atomic vector
  indices_samples.b_sorted <- sort(indices_samples.b_vector)
  
  in_range_sample_a <- scores_factor_one$ID %in% indices_samples.a_sorted
  in_range_sample_b <- scores_factor_one$ID %in% indices_samples.b_sorted
  
  
  simulated_samples_a <- data.frame(sample = scores_factor_one$sample, signal_a = in_range_sample_a)
  simulated_samples_b <- data.frame(sample = scores_factor_one$sample, signal_b = in_range_sample_b)
  
  gfa_omic_one_score <- merge(merge(x=scores_factor_one,y=simulated_samples_a, by="sample", all = TRUE), simulated_samples_b,  by="sample", all = TRUE)
  
  var_false <- var(
    gfa_omic_one_score$score_GFA1[gfa_omic_one_score$signal_a == FALSE & gfa_omic_one_score$signal_b == FALSE], 
    na.rm = TRUE
  )
  
  mean_false <- mean(
    gfa_omic_one_score$score_GFA1[gfa_omic_one_score$signal_a == FALSE & gfa_omic_one_score$signal_b == FALSE], 
    na.rm = TRUE
  )
  
  # Suppress elevated factor
  gfa_omic_one_score$score_GFA1 <- ifelse(
    gfa_omic_one_score$signal_b == TRUE,
    rnorm(length(gfa_omic_one_score$score_GFA1[gfa_omic_one_score$signal_b == TRUE]), mean_false, var_false), # Action if TRUE
    gfa_omic_one_score$score_GFA1# Action if FALSE
  )
  
  # Factor 2
  
  scores_factor_two <- result_score %>%
    mutate(
      ID = ifelse(grepl("sample_", sample),
                  as.numeric(sub(".*sample_(\\d+)", "\\1", sample)), NA)) %>%
    arrange(ID) %>%
    select(ID, sample, score_GFA2)  # Use `select()` to choose specific columns
  
  gfa_omic_two_score <- merge(merge(x=scores_factor_two,y=simulated_samples_a, by="sample", all = TRUE),simulated_samples_b,  by="sample", all = TRUE)
  
  var_false <- var(
    gfa_omic_two_score$score_GFA2[gfa_omic_two_score$signal_a == FALSE & gfa_omic_two_score$signal_b == FALSE], 
    na.rm = TRUE
  )
  
  mean_false <- mean(
    gfa_omic_two_score$score_GFA2[gfa_omic_two_score$signal_a == FALSE & gfa_omic_two_score$signal_b == FALSE], 
    na.rm = TRUE
  )
  
  # Suppress elevated factor
  gfa_omic_two_score$score_GFA2 <- ifelse(
    gfa_omic_two_score$signal_a == TRUE,
    rnorm(length(gfa_omic_two_score$score_GFA2[gfa_omic_two_score$signal_a == TRUE]), mean_false, var_false), # Action if TRUE
    gfa_omic_two_score$score_GFA2# Action if FALSE
  )
  
  # Select only specific columns from Y
  gfa_omic_two_score_selected <- gfa_omic_two_score %>%
    select(sample, score_GFA2)  # Replace 'column1', 'column2' with desired column names
  
  # Perform the merge, keeping all columns from X
  factor_scores <- merge(
    gfa_omic_one_score, 
    gfa_omic_two_score_selected, 
    by = "sample", 
    all.x = TRUE
  )
  
  
  weights_composite <- list(omic.one_weights = gfa_omic_one_load, omic.two_weights = load_omic_two, original_weights = result_load_all)
  scores_composite <- list(factor_scores = factor_scores, original_scores = result_score)
  
  gfa_result <- list(weights = weights_composite, scores = scores_composite)
  
  return(gfa_result)
}
## ----------------------------------- 3.1 GFA Loading Function ---------------------------------------
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

# gfa_score = gfa_score(mfa_object, BC_num = 2)
# gfa_loading = gfa_loading(mfa_object, BC_num = 2)
# gfa_result = func_gfa(data, num = 2)

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

# Define the main simulation function
ssim_MOmicR <- function(n_features_one, n_features_two, n_samples, var_sigma, num_biclusters = 2, num_iterations, method) {
  
  # Assign parameters to global environment
  global_env$n_features_one <- n_features_one
  global_env$n_features_two <- n_features_two
  
  # Initialize result lists
  jaccard_results <- list()
  jaccard_comparison_results <- list()
  per_measures_results <- list()
  dataset_output <- list()
  datasets <- list()
  
  iter <- num_iterations
  # Loop through each value of sigma
  for (sigma in var_sigma) {
    message(paste("Processing sigma =", sigma))
    
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
    
    # Cutoff 2: Based on the ratio of the top value of the features weights/loading or samples scores and the variance of noise in the simulated data
    varphi.two <- function(loading_or_score) {
      numeric_values <- as.numeric(unlist(loading_or_score))
      if (all(!is.na(numeric_values))) {
        varphi <- max(abs(numeric_values)) / sigma
        return(varphi)
      } else {
        return(NA)  # Return NA if not all values are numeric
      }
    }
    
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
          varphi <- top.x / sigma
          return(varphi)
        }
      }
      
      else if (any(grepl("score", vector_name)==TRUE)) {
        if (all(!is.na(numeric_values))) {
          x <- abs(numeric_values)
          top.x <- mean(head(sort(x, decreasing = TRUE), 9))
          varphi <- top.x / sigma
          return(varphi)
        }
      }
      
      stop("Neither 'loading' nor 'score' found in the method argument.")
    }
    
    # Cutoff 5: Based on the true proportion of the samples or features with signal in the matrix
    # varphi.five <- function(loading_or_score, omic=NULL) {
    #   # Obtain the name of the vector
    #   vector_name <- deparse(substitute(loading_or_score))
    #   #loading_or_score=Data_Ft$loading_FABIA
    #   numeric_values <- as.numeric(unlist(loading_or_score))
    #   
    #   if (any(grepl("loading", vector_name)==TRUE)) {
    #     if (method == 'single.factor') {
    #       varphi <- quantile(abs(numeric_values), probs = (1-(length(indices_features.1)/n_features_two)))
    #     } else if (method %in% c('shared.factor', 'random.shared.factor')) {
    #       if(omic == 1){
    #         x = (1-(length(indices_features.1)/n_features_one))
    #         varphi <- quantile(numeric_values, probs = x )[[1]]
    #       } 
    #       else if(omic == 2){
    #         y = (1-(length(indices_features.2)/n_features_two))
    #         varphi <- quantile(numeric_values, probs = y )[[1]]
    #       }
    #       else if(omic == 'all'){
    #         z = (1-(length(union(simulation_result$indices_features.1, simulation_result$indices_features.2))/(n_features_one + n_features_two)))
    #         varphi <- quantile(numeric_values, probs = z )[[1]]
    #       }
    #     } else {
    #       stop("Invalid method specified. Method must be 'single.factor', 'shared.factor' or 'random.shared.factor'.")
    #     }
    #     return(varphi)
    #   }
    #   
    #   else if (any(grepl("score", vector_name)==TRUE)) {
    #     if (method == 'single.factor') {
    #       varphi <- quantile(abs(numeric_values), probs = (1-(length(indices_samples.1)/n_samples)))
    #     } else if (method  %in% c('shared.factor', 'random.shared.factor')) {
    #       y = (1-(length(union(indices_samples.1, indices_samples.2))/n_samples))
    #       varphi <- quantile(numeric_values, probs = y)[[1]]
    #     } else {
    #       stop("Invalid method specified. Method must be 'single.factor', 'shared.factor' or 'random.shared.factor'.")
    #     }
    #     return(varphi)
    #   }
    #   
    #   #stop("Neither 'loading' nor 'score' found in the method argument.")
    # }
    # 
    varphi.five <- function(loading_or_score, omic = NULL, type = c("loading", "score")) {
      type <- match.arg(type)
      numeric_values <- as.numeric(unlist(loading_or_score))
      
      if (type == "loading") {
        if (method == 'single.factor') {
          varphi <- quantile(abs(numeric_values), probs = (1 - (length(indices_features.1c) / n_features_two)))
        } else if (method =='shared.factor') {
          if (omic == 1) {
            x <- (1 - (length(indices_features.1c) / n_features_one))
            varphi <- quantile(abs(numeric_values), probs = x)[[1]]
          } else if (omic == 2) {
            y <- (1 - (length(indices_features.2c) / n_features_two))
            varphi <- quantile(abs(numeric_values), probs = y)[[1]]
          } else if (omic == 'ALL') {
            z <- (1 - (length(union(indices_features.1c, indices_features.2c)) /
                         (n_features_one + n_features_two)))
            varphi <- quantile(abs(numeric_values), probs = z)[[1]]
          }
        } else {
          stop("Invalid method specified.")
        }
        return(varphi)
      }
      
      if (type == "score") {
        if (method == 'single.factor') {
          varphi <- quantile(abs(numeric_values), probs = (1 - (length(indices_samples.1c) / n_samples)))
        } else if (method =='shared.factor') {
          y <- (1 - (length(union(indices_samples.1c, indices_samples.2c)) / n_samples))
          varphi <- quantile(abs(numeric_values), probs = y)[[1]]
        } else {
          stop("Invalid method specified.")
        }
        return(varphi)
      }
      
      stop("Invalid type: must be either 'loading' or 'score'")
    }
    
    # varphi.five <- function(loading_or_score) {
    #   # Obtain the name of the vector
    #   vector_name <- deparse(substitute(loading_or_score))
    #   #loading_or_score=Data_Ft$loading_FABIA
    #   numeric_values <- as.numeric(unlist(loading_or_score))
    #   
    #   if (any(grepl("loading", vector_name)==TRUE)) {
    #     if (method == 'single.factor') {
    #       varphi <- quantile(abs(numeric_values), probs = (1-(length(indices_features.1)/n_features_two)))
    #     } else if (method %in% c('shared.factor', 'random.shared.factor')) {
    #       x = (1-(length(union(indices_features.1, indices_features.2))/(n_features_one + n_features_two)))
    #       varphi <- quantile(numeric_values, probs = x )[[1]]
    #     } else {
    #       stop("Invalid method specified. Method must be 'single.factor', 'shared.factor' or 'random.shared.factor'.")
    #     }
    #     return(varphi)
    #   }
    #   
    #   else if (any(grepl("score", vector_name)==TRUE)) {
    #     if (method == 'single.factor') {
    #       varphi <- quantile(abs(numeric_values), probs = (1-(length(indices_samples.1)/n_samples)))
    #     } else if (method  %in% c('shared.factor', 'random.shared.factor')) {
    #       y = (1-(length(union(indices_samples.1, indices_samples.2))/n_samples))
    #       varphi <- quantile(numeric_values, probs = y)[[1]]
    #     } else {
    #       stop("Invalid method specified. Method must be 'single.factor', 'shared.factor' or 'random.shared.factor'.")
    #     }
    #     return(varphi)
    #   }
    #   
    #   #stop("Neither 'loading' nor 'score' found in the method argument.")
    # }
    # 
    
    
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
      
      # Stop if input doesn't match the expected format
      #stop("Input vector name must contain either 'loading' or 'score' and must be numeric.")
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
      
      # Stop if input doesn't match the expected format
      #stop("Input vector name must contain either 'loading' or 'score' and must be numeric.")
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
    
    if (method == "multiple.factor") { 
      current_seed <- 694#sample(1e6, 1)
      set.seed(current_seed)
      
      #----------------------------- Simulate data -----
      simulated_data <- multiple_factor(
        n_features_one = 4000,#n_features_one,
        n_features_two = 3000,#n_features_two,
        n_samples = n_samples,
        sigmas = sigma,
        iterations = iter,
        n_factors = num_biclusters
      )
      
      for (i in 1:length(simulated_data)) {
        try({
          current_data <- simulated_data[[paste0("iteration_", i)]]$concatenated_datasets[[1]]
          indices_features.1 <- simulated_data[[paste0("iteration_", i)]]$indices_features.1
          indices_features.2 <- simulated_data[[paste0("iteration_", i)]]$indices_features.2
          indices_samples <- simulated_data[[paste0("iteration_", i)]]$indices_samples
          
          print(paste("Processing iteration:", i))

          random.data <- data.frame(t(current_data)); 
          random.data$feature = rownames(random.data) 
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
          indices_features.1a = simulated_data[[paste0("iteration_", i)]]$indices_features.1[1]
          indices_features.1b = simulated_data[[paste0("iteration_", i)]]$indices_features.1[2]
          indices_features.2 = simulated_data[[paste0("iteration_", i)]]$indices_features.2
          
          # FACTOR 1 LOADING
          indices_features.1a_vector <- unlist(indices_features.1a)  # Convert to atomic vector
          indices_features.1a_sorted <- sort(indices_features.1a_vector)
          
          
          indices_features.2_vector <- unlist(indices_features.2)  # Convert to atomic vector
          indices_features.2_sorted <- sort(indices_features.2_vector)
          
          indices_features.1b_vector <- unlist(indices_features.1b)  # Convert to atomic vector
          indices_features.1b_sorted <- sort(indices_features.1b_vector)
          
          in_range.a <- random.data.a$ID %in% indices_features.1a_sorted
          
          simulated_features_a <- data.frame(feature = random.data.a$feature, signal_a = in_range.a)
          
          in_range.c <- random.data.b$ID %in% indices_features.2_sorted
          
          simulated_features_c <- data.frame(feature = random.data.b$feature, signal_a = in_range.c)
          
          simulated_features_sigma_a = rbind(simulated_features_a, simulated_features_c)
          
          indices_features.1b_vector <- unlist(indices_features.1b)  # Convert to atomic vector
          indices_features.1b_sorted <- sort(indices_features.1b_vector)
          
          indices_features.2_vector <- unlist(indices_features.2)  # Convert to atomic vector
          indices_features.2_sorted <- sort(indices_features.2_vector)
          
          in_range.b <- random.data.a$ID %in% indices_features.1b_sorted
          
          simulated_features_b <- data.frame(feature = random.data.a$feature, signal_b = in_range.b)
          simulated_features_c2 <- data.frame(feature = random.data.b$feature, signal_b = FALSE)
          
          simulated_features_sigma_b = rbind(simulated_features_b, simulated_features_c2)
          
          # FACTOR 1 SCORES
          indices_samples.a = simulated_data[[paste0("iteration_", i)]]$indices_samples[1]
          indices_samples.b = simulated_data[[paste0("iteration_", i)]]$indices_samples[2]
          
          indices_samples.a_vector <- unlist(indices_samples.a)  # Convert to atomic vector
          indices_samples.a_sorted <- sort(indices_samples.a_vector)
          indices_samples.b_vector <- unlist(indices_samples.b)  # Convert to atomic vector
          indices_samples.b_sorted <- sort(indices_samples.b_vector)
        
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
          
          
          # Run factorization methods
          fabia_result <- func_fabia(current_data, BC_num = num_biclusters)
          mofa_result <- func_mofa(current_data, num_biclusters)
          mfa_result <- func_mfa(current_data, num_biclusters)
          gfa_result <- func_gfa(current_data, num_biclusters)
          
          # Extract loadings
          # 1. LOADING FACTOR 1 - OMIC ONE
          loading_data.not.scaled.omic.one.f1 <- data.frame(
            feature = mofa_result[["weights"]][["omic.one_weights"]][["feature"]],
            loading_FABIA1 = fabia_result[["weights"]][["omic.one_weights"]][["loading_FABIA2"]],
            loading_MOFA1 = mofa_result[["weights"]][["omic.one_weights"]][["loading_MOFA2"]],
            loading_MFA1 = mfa_result[["weights"]][["omic.one_weights"]][["loading_MFA1"]],
            loading_GFA1 = gfa_result[["weights"]][["omic.one_weights"]][["loading_GFA1"]]
          )
          
          #Convert all numeric columns to absolute values
          abs.loading_data.omic.one.f1 <- loading_data.not.scaled.omic.one.f1 %>% 
            mutate(across(where(is.numeric), ~abs(.))) %>%
            mutate(dataview = ifelse(grepl("omic1_", feature), "omic.one", ifelse(grepl("omic2_", feature), "omic.two", NA)),
                   ID = ifelse(grepl("omic1_", feature),
                               sub(".*omic1_feature_(\\d+)", "\\1", feature),
                               ifelse(grepl("omic2_", feature),
                                      sub(".*omic2_feature_(\\d+)", "\\1", feature),NA))) %>% arrange(ID)
          
          Data_Feature_Omic1.f1 <- merge(x=abs.loading_data.omic.one.f1,y=simulated_features_a, by="feature", all = TRUE)
          #Data_Feature_f1 <- data.frame(Data_Feature_f1)[, c('feature','loading_FABIA1','loading_MOFA1','loading_MFA1','loading_GFA1','signal_a')]
          
          
          # 2. LOADING FACTOR 2 - OMIC ONE
          loading_data.not.scaled.omic.one.f2 <- data.frame(
            feature = mofa_result[["weights"]][["omic.one_weights"]][["feature"]],
            loading_FABIA2 = fabia_result[["weights"]][["omic.one_weights"]][["loading_FABIA1"]],
            loading_MOFA2 = mofa_result[["weights"]][["omic.one_weights"]][["loading_MOFA1"]],
            loading_MFA2 = mfa_result[["weights"]][["omic.one_weights"]][["loading_MFA2"]],
            loading_GFA2 = gfa_result[["weights"]][["omic.one_weights"]][["loading_GFA2"]]
          )
          
          #Convert all numeric columns to absolute values
          abs.loading_data.omic.one.f2 <- loading_data.not.scaled.omic.one.f2 %>% 
            mutate(across(where(is.numeric), ~abs(.))) %>%
            mutate(dataview = ifelse(grepl("omic1_", feature), "omic.one", ifelse(grepl("omic2_", feature), "omic.two", NA)),
                   ID = ifelse(grepl("omic1_", feature),
                               sub(".*omic1_feature_(\\d+)", "\\1", feature),
                               ifelse(grepl("omic2_", feature),
                                      sub(".*omic2_feature_(\\d+)", "\\1", feature),NA))) %>% arrange(ID)
          
          Data_Feature_Omic1.f2 <- merge(x=abs.loading_data.omic.one.f2,y=simulated_features_b, by="feature", all = TRUE)
          Data_Feature_Omic1.f2 <- Data_Feature_Omic1.f2 %>%
            select('feature','loading_FABIA2','loading_MOFA2','loading_MFA2','loading_GFA2','signal_b')
          
          # 3. LOADING FACTOR 1 - OMIC TWO
          loading_data.not.scaled.omic.two.f1 <- data.frame(
            feature = mofa_result[["weights"]][["omic.two_weights"]][["feature"]],
            loading_FABIA1 = fabia_result[["weights"]][["omic.two_weights"]][["loading_FABIA2"]],
            loading_MOFA1 = mofa_result[["weights"]][["omic.two_weights"]][["loading_MOFA1"]],
            loading_MFA1 = mfa_result[["weights"]][["omic.two_weights"]][["loading_MFA1"]],
            loading_GFA1 = gfa_result[["weights"]][["omic.two_weights"]][["loading_GFA1"]]
          )
          
          #Convert all numeric columns to absolute values
          abs.loading_data.omic.two.f1 <- loading_data.not.scaled.omic.two.f1 %>% 
            mutate(across(where(is.numeric), ~abs(.))) %>%
            mutate(dataview = ifelse(grepl("omic1_", feature), "omic.one", ifelse(grepl("omic2_", feature), "omic.two", NA)),
                   ID = ifelse(grepl("omic1_", feature),
                               sub(".*omic1_feature_(\\d+)", "\\1", feature),
                               ifelse(grepl("omic2_", feature),
                                      sub(".*omic2_feature_(\\d+)", "\\1", feature),NA))) %>% arrange(ID)
          
          Data_Feature_Omic2.f1 <- merge(x=abs.loading_data.omic.two.f1,y=simulated_features_c, by="feature", all = TRUE)
          
          
          FeatureData.Omic1 <- merge(Data_Feature_Omic1.f1, Data_Feature_Omic1.f2, by = 'feature')
          FeatureData.Omic2 <- Data_Feature_Omic2.f1
          
          # SCORES
          # 1. FACTOR ONE
          score_data.not.scaled1 <- data.frame(
            sample = mofa_result[["scores"]][["scores"]][["sample"]],
            score_FABIA1 = fabia_result[["scores"]][["scores"]][["score_FABIA2"]],
            score_MOFA1 = mofa_result[["scores"]][["scores"]][["score_MOFA1"]],
            score_MFA1 = mfa_result[["scores"]][["factor_scores"]][["score_MFA2"]],
            score_GFA1 = gfa_result[["scores"]][["factor_scores"]][["score_GFA2"]]
          )
          
          #Convert all numeric columns to absolute values
          score_data.abs1 <- score_data.not.scaled1 %>% 
            mutate(across(where(is.numeric), ~abs(.)))%>%
            mutate(ID = ifelse(grepl("sample_", sample), sub(".*sample_(\\d+)", "\\1", sample),NA)) %>% arrange(ID)
          
          Data_Sample_f1 <- merge(x=score_data.abs1,y=simulated_samples_sigma_a, by="sample", all = TRUE)
          #Data_Sample_f1 <- data.frame(Data_Sample_f1)[, c('sample','score_FABIA1','score_MOFA1','score_MFA1','score_GFA1','signal_a')]
          
          # 2. FACTOR TWO
          score_data.not.scaled2 <- data.frame(
            sample = mofa_result[["scores"]][["factor_two_scores"]][["sample"]],
            score_FABIA2 = fabia_result[["scores"]][["factor_one_scores"]][["score_FABIA1"]],
            score_MOFA2 = mofa_result[["scores"]][["factor_two_scores"]][["score_MOFA2"]],
            score_MFA2 = mfa_result[["scores"]][["factor_scores"]][["score_MFA1"]],
            score_GFA2 = gfa_result[["scores"]][["factor_scores"]][["score_GFA1"]]
          )
          
          #Convert all numeric columns to absolute values
          score_data.abs2 <- score_data.not.scaled2  %>% 
            mutate(across(where(is.numeric), ~abs(.)))%>%
            mutate(ID = ifelse(grepl("sample_", sample), sub(".*sample_(\\d+)", "\\1", sample),NA)) %>% arrange(ID)
          
          Data_Sample_f2 <- merge(x=score_data.abs2,y=simulated_samples_sigma_b, by="sample", all = TRUE)
          #Data_Sample_f2 <- data.frame(Data_Sample_f2)[, c('sample','score_FABIA2','score_MOFA2','score_MFA2','score_GFA2','signal_b')]
          
          
          SampleData <- merge(Data_Sample_f1, Data_Sample_f2, by = 'sample')
          
          
          # # Check if mean of FABIA loadings and scores is zero
          # mean_fabia_loading <- mean(fabia_result$fab_loading$loading_FABIA, na.rm = TRUE)
          # mean_fabia_score <- mean(fabia_result$fab_score$score_FABIA, na.rm = TRUE)
          
          # # Convert data to absolute values
          # 
          # loading_data.abs <- loading_data.not.scaled %>% mutate(across(where(is.numeric), ~abs(.)))
          # 
          # loading_data <- loading_data.abs #normalize_columns(loading_data.abs, columns_to_loadings_normalize)
          # score_data <- score_data.abs #normalize_columns(score_data.abs, columns_to_scores_normalize)
          # 
          # Data_Ft <- merge(x=loading_data,y=simulated_features_sigma, by="feature", all = TRUE)
          # Data_Ft <- data.frame(Data_Ft)[, c('feature','loading_FABIA','loading_MOFA','loading_MFA','loading_GFA','signal')]
          
          # Mutate the dataview column based on the pattern in the feature column
          # FeatureData <- FeatureData %>%
          #   mutate(dataview = ifelse(grepl("omic1_", feature), "omic.one", ifelse(grepl("omic2_", feature), "omic.two", NA)))
          # 
          # FeatureData$ID <- ifelse(grepl("omic1_", FeatureData$feature),
          #                          sub(".*omic1_feature_(\\d+)", "\\1", FeatureData$feature),
          #                          ifelse(grepl("omic2_", FeatureData$feature),
          #                                 sub(".*omic2_feature_(\\d+)", "\\1", FeatureData$feature),NA))
          # 
          
          # Subset omic1
          # FeatureData.a <- FeatureData.Omic1 %>%
          #   filter(dataview == "omic.one") %>% arrange(as.numeric(ID))
          # 
          # # Subset omic2
          # FeatureData.b <- FeatureData %>%
          #   filter(dataview == "omic.two") %>% arrange(as.numeric(ID))
          # 
          for (varphi in seq_along(varphi_functions)) {
            
            varphi_function <- names(varphi_functions)[[varphi]]
            
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
              ) %>% arrange(as.numeric(ID))
            
            FeatureData.b.two <- FeatureData.Omic2 %>%
              mutate(
                fabia_index_a = ifelse(abs(loading_FABIA1) >= get(varphi_function)(loading_FABIA1), 1, 0),
                mofa_index_a = ifelse(abs(loading_MOFA1) >= get(varphi_function)(loading_MOFA1), 1, 0),
                mfa_index_a = ifelse(abs(loading_MFA1) >= get(varphi_function)(loading_MFA1), 1, 0),
                gfa_index_a = ifelse(abs(loading_GFA1) >= get(varphi_function)(loading_GFA1), 1, 0),
                signal_index_a = ifelse(signal_a == 'TRUE', 1, 0)#,
                #fabia_index_b = ifelse(abs(loading_FABIA2) >= get(varphi_function)(loading_FABIA2), 1, 0)#,
                #mofa_index_b = ifelse(abs(loading_MOFA2) >= get(varphi_function)(loading_MOFA2), 1, 0),
                #mfa_index_b = ifelse(abs(loading_MFA2) >= get(varphi_function)(loading_MFA2), 1, 0),
                #gfa_index_b = ifelse(abs(loading_GFA2) >= get(varphi_function)(loading_GFA2), 1, 0),
                #signal_index_b = ifelse(signal_b == 'TRUE', 1, 0)
              ) %>% arrange(as.numeric(ID))
            
            # FeatureData2 <- FeatureData %>%
            #   mutate(
            #     fabia_index_a = ifelse(abs(loading_FABIA1) >= get(varphi_function)(loading_FABIA1), 1, 0),
            #     mofa_index_a = ifelse(abs(loading_MOFA1) >= get(varphi_function)(loading_MOFA1), 1, 0),
            #     mfa_index_a = ifelse(abs(loading_MFA1) >= get(varphi_function)(loading_MFA1), 1, 0),
            #     gfa_index_a = ifelse(abs(loading_GFA1) >= get(varphi_function)(loading_GFA1), 1, 0),
            #     signal_index_a = ifelse(signal_a == 'TRUE', 1, 0),
            #     fabia_index_b = ifelse(abs(loading_FABIA2) >= get(varphi_function)(loading_FABIA2), 1, 0),
            #     mofa_index_b = ifelse(abs(loading_MOFA2) >= get(varphi_function)(loading_MOFA2), 1, 0),
            #     mfa_index_b = ifelse(abs(loading_MFA2) >= get(varphi_function)(loading_MFA2), 1, 0),
            #     gfa_index_b = ifelse(abs(loading_GFA2) >= get(varphi_function)(loading_GFA2), 1, 0),
            #     signal_index_b = ifelse(signal_b == 'TRUE', 1, 0)
            #   ) %>% arrange(as.numeric(ID))
            
            SampleData2 <- SampleData %>%
              mutate(
                fabia_index_a = ifelse(abs(score_FABIA1) >= get(varphi_function)(score_FABIA1), 1, 0),
                mofa_index_a = ifelse(abs(score_MOFA1) >= get(varphi_function)(score_MOFA1), 1, 0),
                mfa_index_a = ifelse(abs(score_MFA1) >= get(varphi_function)(score_MFA1), 1, 0),
                gfa_index_a = ifelse(abs(score_GFA1) >= get(varphi_function)(score_GFA1), 1, 0),
                signal_index_a = ifelse(signal_a == 'TRUE', 1, 0),
                fabia_index_b = ifelse(abs(score_FABIA1) >= get(varphi_function)(score_FABIA1), 1, 0),
                mofa_index_b = ifelse(abs(score_MOFA1) >= get(varphi_function)(score_MOFA1), 1, 0),
                mfa_index_b = ifelse(abs(score_MFA1) >= get(varphi_function)(score_MFA1), 1, 0),
                gfa_index_b = ifelse(abs(score_GFA1) >= get(varphi_function)(score_GFA1), 1, 0),
                signal_index_b = ifelse(signal_b == 'TRUE', 1, 0)
              ) %>% arrange(ID.x)
            
            # # Create binary indices
            # Data_Matrix_Ft <- Data_Ft %>%
            #   mutate(fabia_index = ifelse(abs(loading_FABIA) >= get(varphi_function)(loading_FABIA), 1, 0),
            #          mofa_index = ifelse(abs(loading_MOFA) >= get(varphi_function)(loading_MOFA), 1, 0),
            #          mfa_index = ifelse(abs(loading_MFA) >= get(varphi_function)(loading_MFA), 1, 0),
            #          gfa_index = ifelse(abs(loading_GFA) >= get(varphi_function)(loading_GFA), 1, 0),
            #          signal_index = ifelse(signal == 'TRUE', 1, 0))
            # 
            # # Create binary indices
            # Data_Matrix_Smp <- Data_Smp %>% #DataTRUE_Smp %>%
            #   mutate(fabia_index = ifelse(abs(score_FABIA) >= get(varphi_function)(score_FABIA), 1, 0),
            #          mofa_index = ifelse(abs(score_MOFA) >= get(varphi_function)(score_MOFA), 1, 0),
            #          mfa_index = ifelse(abs(score_MFA) >= get(varphi_function)(score_MFA), 1, 0),
            #          gfa_index = ifelse(abs(score_GFA) >= get(varphi_function)(score_GFA), 1, 0),
            #          signal_index = ifelse(signal == 'TRUE', 1, 0))
            # 
            # # Signal MATRIX
            # sample_signal_index = Data_Matrix_Smp$signal_index
            # feature_signal_index = Data_Matrix_Ft$signal_index
            # signal_matrix = as.matrix(feature_signal_index)%*%as.matrix(t(sample_signal_index))
            # 
            # # loading_score_matrix
            # fabia_matrix = as.matrix(Data_Matrix_Ft$fabia_index)%*%as.matrix(t(Data_Matrix_Smp$fabia_index))
            # mofa_matrix = as.matrix(Data_Matrix_Ft$mofa_index)%*%as.matrix(t(Data_Matrix_Smp$mofa_index))
            # mfa_matrix = as.matrix(Data_Matrix_Ft$mfa_index)%*%as.matrix(t(Data_Matrix_Smp$mfa_index))
            # gfa_matrix = as.matrix(Data_Matrix_Ft$gfa_index)%*%as.matrix(t(Data_Matrix_Smp$gfa_index))
            # 
            
            #### ------------------------------------ Jaccard Index ------------------------------------------ #
            
            # List of methods and corresponding indices
            methods <- c("fabia", "mofa", "mfa", "gfa")
            
            # Compute Jaccard indices in a single step
            # ji_true_a <- sapply(methods, function(method) {
            #   method_index <- FeatureData2[[paste0(method, "_index_a")]]
            #   jaccard.index.sim(method_index, FeatureData2$signal_index_a)
            # }, simplify = TRUE, USE.NAMES = TRUE)
            
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
              method_index <- SampleData2[[paste0(method, "_index_a")]]
              jaccard.index.sim(method_index, SampleData2$signal_index_a)
            }, simplify = TRUE, USE.NAMES = TRUE)
            
            # Compute Jaccard indices in a single step
            # ji_true_b <- sapply(methods, function(method) {
            #   method_index <- FeatureData2[[paste0(method, "_index_b")]]
            #   jaccard.index.sim(method_index, FeatureData2$signal_index_b)
            # }, simplify = TRUE, USE.NAMES = TRUE)
            
            ji_omic_one_b <- sapply(methods, function(method) {
              method_index <- FeatureData.a.one[[paste0(method, "_index_b")]]
              jaccard.index.sim(method_index, FeatureData.a.one$signal_index_b)
            }, simplify = TRUE, USE.NAMES = TRUE)
            
            ji_omic_two_b <- sapply(methods, function(method) {
              method_index <- FeatureData.b.two[[paste0(method, "_index_b")]]
              jaccard.index.sim(method_index, FeatureData.b.two$signal_index_b)
            }, simplify = TRUE, USE.NAMES = TRUE)
            
            ji_samples_b <- sapply(methods, function(method) {
              method_index <- SampleData2[[paste0(method, "_index_b")]]
              jaccard.index.sim(method_index, SampleData2$signal_index_b)
            }, simplify = TRUE, USE.NAMES = TRUE)
            
            # # Define matrices
            # matrices <- list(
            #   fabia = fabia_matrix,
            #   mofa = mofa_matrix,
            #   mfa = mfa_matrix,
            #   gfa = gfa_matrix
            # )
            # 
            # # Compute Jaccard similarity indices
            # ji_matrix <- sapply(names(matrices), function(method) {
            #   jaccard.index.sim(matrices[[method]], signal_matrix)
            # }, simplify = TRUE, USE.NAMES = TRUE)
            # 
            # 
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
            # # Calculate pairwise Jaccard indices dynamically
            # ji_true_pairs_a <- combn(methods, 2, function(method_pair) {
            #   method1_index <- FeatureData2[[paste0(method_pair[1], "_index_a")]]
            #   method2_index <- FeatureData2[[paste0(method_pair[2], "_index_a")]]
            #   jaccard.index.sim(method1_index, method2_index)
            # }, simplify = TRUE)
            # 
            # # Assign names to the results
            # names(ji_true_pairs_a) <- combn(methods, 2, function(method_pair) {
            #   paste(method_pair[1], method_pair[2], "ji_phi", sep = "_")
            # })
            # 
            # # Convert to a named vector or data frame for further use
            # ji_true_pairs_a <- as.data.frame(as.list(ji_true_pairs_a))
            # 
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
              method1_index <- SampleData2[[paste0(method_pair[1], "_index_a")]]
              method2_index <- SampleData2[[paste0(method_pair[2], "_index_a")]]
              jaccard.index.sim(method1_index, method2_index)
            }, simplify = TRUE)
            
            # Assign names to the results
            names(ji_smp_pairs_a) <- combn(methods, 2, function(method_pair) {
              paste(method_pair[1], method_pair[2], "ji_phi", sep = "_")
            })
            
            # Convert to a named vector or data frame for further use
            ji_smp_pairs_a <- as.data.frame(as.list(ji_smp_pairs_a))
            
            # JACCARD INDEX FACTOR 2
            
            # # Calculate pairwise Jaccard indices dynamically
            # ji_true_pairs_b <- combn(methods, 2, function(method_pair) {
            #   method1_index <- FeatureData2[[paste0(method_pair[1], "_index_b")]]
            #   method2_index <- FeatureData2[[paste0(method_pair[2], "_index_b")]]
            #   jaccard.index.sim(method1_index, method2_index)
            # }, simplify = TRUE)
            # 
            # # Assign names to the results
            # names(ji_true_pairs_b) <- combn(methods, 2, function(method_pair) {
            #   paste(method_pair[1], method_pair[2], "ji_phi", sep = "_")
            # })
            # 
            # # Convert to a named vector or data frame for further use
            # ji_true_pairs_b <- as.data.frame(as.list(ji_true_pairs_b))
            # 
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
              method1_index <- SampleData2[[paste0(method_pair[1], "_index_b")]]
              method2_index <- SampleData2[[paste0(method_pair[2], "_index_b")]]
              jaccard.index.sim(method1_index, method2_index)
            }, simplify = TRUE)
            
            # Assign names to the results
            names(ji_smp_pairs_b) <- combn(methods, 2, function(method_pair) {
              paste(method_pair[1], method_pair[2], "ji_phi", sep = "_")
            })
            
            # Convert to a named vector or data frame for further use
            ji_smp_pairs_b <- as.data.frame(as.list(ji_smp_pairs_b))
            
            # MATRICES
            # 
            # # Define the matrices in a named list
            # matrices <- list(
            #   fabia = fabia_matrix,
            #   mofa = mofa_matrix,
            #   mfa = mfa_matrix,
            #   gfa = gfa_matrix
            # )
            # 
            # # Get all unique pair combinations of matrix names
            # matrix_pairs <- combn(names(matrices), 2, simplify = FALSE)
            # 
            # # Compute Jaccard similarity indices for each pair
            # ji_matrix_pairs <- sapply(matrix_pairs, function(pair) {
            #   jaccard.index.sim(matrices[[pair[1]]], matrices[[pair[2]]])
            # }, simplify = TRUE, USE.NAMES = FALSE)
            # 
            # # Assign names to the results
            # names(ji_matrix_pairs) <- sapply(matrix_pairs, function(pair) {
            #   paste(pair[1], pair[2], sep = "_")
            # })
            # 
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
            # Define methods for comparison
            methods_a <- c("fabia_index_a", "mofa_index_a", "mfa_index_a", "gfa_index_a") # factor 1
            methods_b <- c("fabia_index_b", "mofa_index_b", "mfa_index_b", "gfa_index_b")
            
            metric_samples_a <- CompareMetrics(data = SampleData2, signal_index = "signal_index_a", methods = methods_a)
            metric_samples_b <- CompareMetrics(data = SampleData2, signal_index = "signal_index_b", methods = methods_b)
            metric_omic.one_a <- CompareMetrics(data = FeatureData.a.one, signal_index = "signal_index_a", methods = methods_a)
            metric_omic.one_b <- CompareMetrics(data = FeatureData.a.one, signal_index = "signal_index_b", methods = methods_b)
            metric_omic.two_a <- CompareMetrics(data = FeatureData.b.two, signal_index = "signal_index_a", methods = methods_a)
            metric_omic.two_b <- CompareMetrics(data = FeatureData.b.two, signal_index = "signal_index_b", methods = methods_b)
            
            # per_measures_true <- compare_features(Data_Ft2)
            # per_measures_omic_one <- compare_features(Data_Ft.a.one)
            # per_measures_omic_two <- compare_features(Data_Ft.b.two)
            # per_measures_smp <- compare_sample(Data_Smp2)
            # 
            # Organize the results for the current method
            per_measures <- list(metric_samples_a = metric_samples_a, 
                                 metric_samples_b = metric_samples_b, 
                                 metric_omic.one_a = metric_omic.one_a,
                                 metric_omic.one_b = metric_omic.one_b, 
                                 metric_omic.two_a = metric_omic.two_a,
                                 metric_omic.two_b = metric_omic.two_b)
            
            # Create a unique name for each dataset
            #per_measures_results <- list()
            dataset_name <- sprintf("variance_%d_iteration_%d_varphi_%d", sigma, i, varphi) 
            per_measures_results[[dataset_name]] <- per_measures
            
            # simulatedData <- simulated_data[[paste0("iteration_", 1)]]$concatenated_datasets[[1]]
            # 
            # # Create a list of datasets
            # datasets <- list(
            #   simulatedData = simulatedData,
            #   FeatureData2 = FeatureData2,
            #   FeatureData.a.one = FeatureData.a.one,
            #   FeatureData.b.two = FeatureData.b.two,
            #   SampleData2 = SampleData2
            # ) 
          }
          
          datasets <- list(
            simulatedData = current_data,
            FeatureData2 = FeatureData2,
            SampleData2 = SampleData2
          ) 
          dataset <- sprintf("variance_%d_iteration_%d", sigma, i)
          dataset_output[[dataset]] <- datasets
          
        }, silent = TRUE)
      }
      
    } else {
      stop("Invalid method specified. Method must be 'single.factor', 'shared.factor', or 'random.shared.factor'.")
    }
    
  }
  # Combine results into a list
  results <- list(
    jaccard_results = jaccard_results,
    jaccard_comparison_results = jaccard_comparison_results,
    per_measures_results = per_measures_results,
    dataset_output = dataset_output)
  
  return(results)
  
}
  
#  -------------------------- ACTUAL SIMULATION: Define the actual simulation parameters -------------------------- #
  
# Define the methods and sigmas
methods <- c("multiple.factor")
sigmas <- 7  # Define sigma or the range of sigmas

# Loop over sigmas
for (sigma in sigmas) {
  # Loop over methods for the current sigma
  for (method in methods) {
    # Generate a unique result name
    result_name <- paste0("sim_MomicResults_", method, "_sigma_", sigma)
    
    # Run the simulation for the current method and sigma
    sim_result <- sim_MOmicR(
      n_features_one = 4000, 
      n_features_two = 3000, 
      n_samples = 100, 
      var_sigma = sigma, 
      num_biclusters = 2,
      num_iteration = 1, 
      method = methods) 
      
    # Save the result to a file
    saveRDS(
      sim_result,
      #file = paste0("/user/leuven/364/vsc36498/", result_name, ".rds")
      #file = paste0("C:/Users/Lenovo/Downloads/", result_name, ".rds")
      file = paste0("C:/Users/bosangir/Downloads/", result_name, ".rds")
    )                                  
  }
  
  # Print progress message
  cat("Completed simulations for sigma =", sigma, "\n")
} 
  
# Print completion message
cat("All simulations completed and saved!\n")

simulatedData <- multiple_factor(
  n_features_one = 4000,#n_features_one,
  n_features_two = 3000,#n_features_two,
  n_samples = n_samples,
  sigmas = sigma,
  iterations = iter,
  n_factors = num_biclusters
)
DATA = simulatedData[["iteration_1"]][["concatenated_datasets"]][[1]]
image(c(1:dim(DATA)[1]), c(1:dim(DATA)[2]), DATA, ylab = "Features", xlab = "Samples", font.lab = 2)
abline(h = 3000, col = "brown2", lty = "dashed")  # lty = "dashed" sets the line type




  
  
  
  
  



























  


                     
                     
                            



  



  


  
  
  
  
  
    
      
      
      
      
      
      
      
      
      
      
    
  
  
    
      
      
      
      
      
      
      
      
      
      
    
  
  
    
      
      
      
      
      
      
      
      
      
      
    
  
  
    
      
      
      
      
      
      
      
      
      
      
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    
    
  
  
  
    
    
  
  
  
    
    
  
  
  
    
    
  
  
  
  
    
    
  
  
  
    
    
  
  
  
    
    
  
  
  
    
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
                          
  
  
  
  
  
  
  
  

  
  
  
    
    
    
  
  
  
  
    
  
  
  
  
  
  
  
  
    
    
    
  
  
  
  
    
  
  
  
  
  
  
  
  
    
    
    
  
  
  
  
    
  
  
  
  
  
  
  
  
    
    
    
  
  
  
  
    
  
  
  
  
  
  
  
  
  
    
    
    
  
  
  
  
    
  
  
  
  
  
  
  
  
    
    
    
  
  
  
  
    
  
  
  
  
  
  
  
  
    
    
    
  
  
  
  
    
  
  
  
  
  
  
  
  
    
    
    
  
  
  
  
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
                             
                             
                             
  
  
  
  
