setwd("/user/leuven/364/vsc36498")
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
single.factor <- function(n_features_one, n_features_two, sig_omic_one = c(1500, 1965), n_samples, sigmas, iterations) {
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
      valid_data <- FALSE  # Flag to check if data passes FABIA condition
      
      while (!valid_data) {
        # Generate first omic data
        s_sig_s <- ceiling(n_s / 5.3)
        s_sig_e <- ceiling(n_s / 2.8)
        alpha <- rnorm(n_s, 0, 0.05)
        assigned_indices_samples <- sample(s_sig_s:s_sig_e, length(s_sig_s:s_sig_e))
        alpha[assigned_indices_samples] <- rnorm(length(s_sig_s:s_sig_e), 4.5, 0.05)
        d_sig_s <- sig_omic_one[1]#ceiling(n_features_one / 2.11)
        d_sig_e <- sig_omic_one[2]#ceiling(n_features_one / 1.65)
        beta <- rnorm(n_features_one, 0, 0.05)
        assigned_indices_features <- sample(d_sig_s:d_sig_e, length(d_sig_s:d_sig_e))
        beta[assigned_indices_features] <- rnorm(length(d_sig_s:d_sig_e), 5, 0.05)
        data.1 <- alpha %*% t(beta)
        eps <- rnorm(n_s * n_features_one, 0, sigmas[k])
        omic1_data <- matrix(data.1, n_s, n_features_one) + matrix(eps, n_s, n_features_one)
        colnames(omic1_data) <- paste0('omic1_feature_', seq_len(n_features_one))
        rownames(omic1_data) <- paste0('sample_', seq_len(n_s))
        
        # Generate second omic data
        gamma <- rnorm(n_s, 0, 0.05)
        assigned_indices_samples2 <- sample(s_sig_s:s_sig_e, length(s_sig_s:s_sig_e))
        gamma[assigned_indices_samples2] <- rnorm(length(s_sig_s:s_sig_e), 0, 0.05)
        d_sig2_s <- ceiling(n_features_two / 1.5)
        d_sig2_e <- ceiling(n_features_two / sample(8:9, 1))
        delta <- rnorm(n_features_two, 0, 0.05)
        assigned_indices_features2 <- sample(d_sig2_s:d_sig2_e, length(d_sig2_s:d_sig2_e))
        delta[assigned_indices_features2] <- rnorm(length(d_sig2_s:d_sig2_e), 0, 0.05)
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
        
        # Validate FABIA results
        if (!any(is.na(fab_loading)) && mean(fab_loading, na.rm = TRUE) != 0 &&
            !any(is.na(fab_score)) && mean(fab_score, na.rm = TRUE) != 0) {
          valid_data <- TRUE
        } else {
          message("FABIA validation failed, regenerating concatenated_data...")
        }
      }
      
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
      indices_samples.1 = assigned_indices_samples,
      indices_features.2 = assigned_indices_features2,
      indices_samples.2 = assigned_indices_samples2,
      sample_sig_start = s_sig_s,
      sample_sig_end = s_sig_e,
      feature_sig_start = d_sig_s,
      feature_sig_end = d_sig_e,
      feature_sig2_start = d_sig2_s,
      feature_sig2_end = d_sig2_e
    )
  }
  
  return(all_omic_data)
}

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
        alpha <- rnorm(n_s, 0, 0.05)
        assigned_indices_samples <- sample(s_sig_s:s_sig_e, length(s_sig_s:s_sig_e))
        alpha[assigned_indices_samples] <- rnorm(length(s_sig_s:s_sig_e), 4.5, 0.05)
        d_sig_s <- sig_omic_one[1]#ceiling(n_features_one / 1.11)
        d_sig_e <- ceiling(n_features_one / 1)
        beta <- rnorm(n_features_one, 0, 0.05)
        assigned_indices_features <- sample(d_sig_s:d_sig_e, length(d_sig_s:d_sig_e))
        beta[assigned_indices_features] <- rnorm(length(d_sig_s:d_sig_e), 5, 0.05)
        data.1 <- alpha %*% t(beta)
        eps <- rnorm(n_s * n_features_one, 0, sigmas[k])
        omic1_data <- matrix(data.1, n_s, n_features_one) + matrix(eps, n_s, n_features_one)
        colnames(omic1_data) <- paste0('omic1_feature_', seq_len(n_features_one))
        rownames(omic1_data) <- paste0('sample_', seq_len(n_s))
        
        # Generate second OMIC data
        gamma <- rnorm(n_s, 0, 0.05)
        assigned_indices_samples2 <- sample(s_sig_s:s_sig_e, length(s_sig_s:s_sig_e))
        gamma[assigned_indices_samples2] <- rnorm(length(s_sig_s:s_sig_e), 4.5, 0.05)
        d_sig2_s <- 1#ceiling(n_features_two / 1)
        d_sig2_e <- 335#ceiling(n_features_two + 334)# / sample(8:9, 1))
        delta <- rnorm(n_features_two, 0, 0.05)
        assigned_indices_features2 <- sample(d_sig2_s:d_sig2_e, length(d_sig2_s:d_sig2_e))
        delta[assigned_indices_features2] <- rnorm(length(d_sig2_s:d_sig2_e), 5.5, 0.05)
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
      feature_sig2_end = d_sig2_e
    )
  }
  
  return(all_omic_data)
}

# random.shared.factor function
random.shared.factor <- function(n_features_one, sig_omic_one = c(2250, 2550), n_features_two, sig_omic_two = c(1500, 1850), n_samples, sigmas, iterations) {
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
  
  for (iter in 1:iterations) {  # Outer loop for iterations
    omic.one <- list()
    omic.two <- list()
    
    for (k in 1:length(sigmas)) {  # Inner loop for sigmas
      n_s <- n_samples
      valid_data <- FALSE  # Flag for FABIA validation
      
      while (!valid_data) {
        # First OMIC data
        s_sig_s <- ceiling(n_s / 5.3)
        s_sig_e <- ceiling(n_s / 2.8)
        alpha <- rnorm(n_s, 0, 0.05)
        assigned_indices_samples <- sample(1:n_samples, length(s_sig_s:s_sig_e))
        alpha[assigned_indices_samples] <- rnorm(length(s_sig_s:s_sig_e), 4.5, 0.05)
        d_sig_s <- ceiling(n_features_one / 1.11)
        d_sig_e <- ceiling(n_features_one / 1)
        beta <- rnorm(n_features_one, 0, 0.05)
        assigned_indices_features <- sample(1:n_features_one, length(d_sig_s:d_sig_e))
        beta[assigned_indices_features] <- rnorm(length(d_sig_s:d_sig_e), 5, 0.05)
        data.1 <- alpha %*% t(beta)
        eps <- rnorm(n_s * n_features_one, 0, sigmas[k])
        omic1_data <- matrix(data.1, n_s, n_features_one) + matrix(eps, n_s, n_features_one)
        colnames(omic1_data) <- paste0('omic1_feature_', seq_len(n_features_one))
        rownames(omic1_data) <- paste0('sample_', seq_len(n_s))
        
        # Second OMIC data
        s_sig2_s <- ceiling(n_s / 5.3)
        s_sig2_e <- ceiling(n_s / 2.2)
        gamma <- rnorm(n_s, 0, 0.05)
        assigned_indices_samples2 <- assigned_indices_samples
        replacement_length <- length(assigned_indices_samples2)
        gamma[assigned_indices_samples2] <- rnorm(replacement_length, 4.5, 0.05)
        #gamma[assigned_indices_samples2] <- rnorm(length(s_sig2_s:s_sig2_e), 3, 0.05)
        d_sig2_s <- 1#ceiling(n_features_two / 1)
        d_sig2_e <- 335#ceiling(n_features_two + 334)# / sample(8:9, 1))
        delta <- rnorm(n_features_two, 0, 0.05)
        assigned_indices_features2 <- sample(1:n_features_two, length(d_sig2_s:d_sig2_e))
        delta[assigned_indices_features2] <- rnorm(length(d_sig2_s:d_sig2_e), 5.5, 0.05)
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
      feature_sig2_end = d_sig2_e
    )
  }
  
  return(all_omic_data)
}

# -------------------------- FABIA: Define the function to run FABIA -----------------------------
fab_score <- function(fabia_object, BC_num) { # FABIA factor score function
  score_FABIA <- fabia_object@Z[BC_num, ]
  BC_df_score <- as.data.frame(score_FABIA)
  BC_df_score$sample = rownames(BC_df_score)
  return(BC_df_score)
}

fab_loading <- function(fabia_object, BC_num) { # FABIA features loadings/weights function
  loading_FABIA <- fabia_object@L[,BC_num]
  BC_df_loading <- as.data.frame(loading_FABIA)
  BC_df_loading$feature <- rownames(BC_df_loading)
  return(BC_df_loading)
}

# FABIA function
# Updated func_fabia Function
func_fabia <- function(data, num) {
  valid_data <- FALSE  # Initialize validation flag
  
  while (!valid_data) {
    #set.seed(123)
    
    # Run FABIA
    fabia_object <- fabia(as.matrix(t(data)), p = 1, alpha = 0.0, cyc = 1000, spl = 0.5, spz = 0.5, 
                          random = 1.0, center = 2, norm = 2, lap = 1.0, nL = 1)
    
    # List container for feature weights and factor score results
    fab_loading <- list()
    fab_score <- list()
    
    # Extract feature weights and factor scores
    # Iterate over the number of biclusters
    fab_loading <- fab_loading(fabia_object, num)
    fab_score <- fab_score(fabia_object, num)
    
    # Validate FABIA results
    if (!any(is.na(fab_loading$loading_FABIA)) && mean(fab_loading$loading_FABIA, na.rm = TRUE) != 0 &&
        !any(is.na(fab_score$score_FABIA)) && mean(fab_score$score_FABIA, na.rm = TRUE) != 0) {
      valid_data <- TRUE  # Data passes validation
    } else {
      message("FABIA validation failed, re-run FABIA...")
    }
  }
  
  # FABIA normalized data X
  X <- fabia_object@X
  
  # Return results
  fabia_result <- list(fab_loading = fab_loading, fab_score = fab_score, X = X)
  return(fabia_result)
}
#xx = simulated_data$iteration_1$concatenated_datasets[[1]]
#result =func_fabia(xx, num = 1) 
# func_fabia <- function(data, num) {
#   set.seed(123)
#   # Run FABIA
#   fabia_object <- fabia(as.matrix(data), p = 1, alpha = 0.01, cyc = 1000, spl = 0.5, spz = 0.5, 
#                         random = 1.0, center = 2, norm = 2, lap = 1.0, nL = 1)
#   
#   # List container for feature weights and factor score results
#   fab_loading <- list()
#   fab_score <- list()
#   
#   # Iterate over the number of biclusters
#   fab_loading <- fab_loading(fabia_object, num)
#   fab_score <- fab_score(fabia_object, num)
#   
#   # Create functions to obtain fabia normalized data X 
#   X = fabia_object@X 
#   fabia_result <- list(fab_loading=fab_loading, fab_score=fab_score, X=X)
#   return(fabia_result)
# }

# -------------------------- MOFA: Define the function to run MOFA -----------------------------
# (a) MOFA factor score function # Create sample score data frame 
mofa_score <- function(model, factor_num) {
  score_mofa <- get_factors(model, factors = factor_num, as.data.frame = T)
  score_MOFA <- score_mofa[, c("sample", "value")]
  colnames(score_MOFA) <- c("sample","score_MOFA")
  return(score_MOFA)
} 

# (b) MOFA features loadings/weights function # Create MOFA feature loading data frame
mofa_loading <- function(model, factor_num) {
  loading_mofa <- get_weights(model, factors = factor_num, as.data.frame = T)
  loading_MOFA <- loading_mofa[, c("feature", "value")]
  colnames(loading_MOFA) <- c("feature","loading_MOFA")
  return(loading_MOFA)
}
# MOFA function
func_mofa <- function(data, num) {
  #set.seed(123)
  library(basilisk)
  mofa_data <- data #sim_Omic_X[[1]]
  
  # Access variables from the shared environment
  n_features_one <- global_env$n_features_one
  n_features_two <- global_env$n_features_two
  
  # Subset respective datasets
  first_omic <- mofa_data[, 1:n_features_one]
  second_omic <- mofa_data[, (n_features_one + 1):(n_features_one + n_features_two)]
  
  # Create a data list 'ready' for MOFA
  mofa_data_sim <- list(first_omic = as.matrix(t(first_omic)), second_omic = as.matrix(t(second_omic)))
  
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
  train_opts_sim$convergence_mode <- "slow"
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
      use_basilisk = TRUE)
    model_object <- MOFAobject_sim.trained
    
    # Lists to store loading and scores for all factors
    mofa_loading <- list()
    mofa_score <- list()
    
    # Iterate over the number of factors
    mofa_loading <- mofa_loading(model_object, num)
    mofa_score <- mofa_score(model_object, num)
    
    mofa_result <- list(mofa_loading = mofa_loading, mofa_score = mofa_score)
    return(mofa_result)
    .quality_control(object, verbose = verbose)
    max(tmp, na.rm = TRUE)
  })
}

# -------------------------- MFA: Define the function to run MFA -----------------------------
# (a) MFA factor score function # Create sample score data frame 
mfa_score <- function(mfa_object, num){
  score_MFA = mfa_object$ind$coord[,num]
  score_MFA_df = as.data.frame(score_MFA)
  return(score_MFA_df)
}

# (b) MFA features loadings/weights function # Create MOFA feature loading data frame
mfa_loading <- function(mfa_object, num){
  loading_MFA = mfa_object$quanti.var$coord[,num]
  loading_MFA_df = as.data.frame(loading_MFA)
  return(loading_MFA_df)
}

# MFA function
func_mfa <- function(data, num) {
  #set.seed(123)
  
  # Access variables from the shared environment
  n_features_one <- global_env$n_features_one
  n_features_two <- global_env$n_features_two
  
  mfa_data <- data
  
  # mfa analysis
  mfa_data = data.frame(mfa_data)
  
  # Create a list specifying which columns are quantitative or qualitative
  
  mfa_object = MFA(as.matrix(mfa_data), 
                   group = c(n_features_one, n_features_two), 
                   type = c("s","s"), 
                   name.group = c("first.omic", "second.omic"),
                   graph = FALSE)
  
  mfa_loading <- list()
  mfa_score <- list()
  
  mfa_loading <- mfa_loading(mfa_object, num)
  mfa_score <- mfa_score(mfa_object, num)
  
  mfa_result <- list(mfa_loading = mfa_loading, mfa_score = mfa_score)
  return(mfa_result)
}

# -------------------------- GFA: Define the function to run GFA -----------------------------
# (a) GFA factor score function # Create sample score data frame
gfa_loading <- function(gfa_object, BC_num){
  loading_GFA = as.data.frame(gfa_object$X)
  loading_GFA_df= data.frame(loading_GFA$V)
  rownames(loading_GFA_df)= rownames(loading_GFA)
  return(loading_GFA_df)
}

# (b) MFA features loadings/weights function # Create MOFA feature loading data frame
gfa_score <- function(gfa_object, BC_num){
  score_GFA = as.data.frame(gfa_object$W)
  score_GFA_df= data.frame(score_GFA$V1)
  rownames(score_GFA_df)= rownames(score_GFA)
  return(score_GFA_df)
}

# GFA function
func_gfa <- function(data, num){
  #set.seed(123)
  gfa_dt = as.data.frame(data)
  merged_GFA_data = list(t(gfa_dt))
  model_option <- getDefaultOpts()
  model_option$iter.max <- 1000
  model_option$iter.burnin <- 10
  gfa_object <- gfa(t(merged_GFA_data), K= num, opts=model_option)
  #normalized_data <- normalizeData(merged_GFA_data, type="center")
  #visulization <- visualizeComponents(gfa_object, merged_GFA_data, normalized_data )
  
  gfa_loading <- list()
  gfa_score <- list()
  
  gfa_loading <- gfa_loading(gfa_object, num)
  gfa_score <- gfa_score(gfa_object, num)
  
  gfa_result <- list(gfa_loading = gfa_loading, gfa_score = gfa_score)
  return(gfa_result)
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

# -------------------------- PERFORMANCE: Define the performance metrics: features -----------------------------
# - Performance analysis
# compare_features <- function(data) {
#   calculate_performance_metrics <- function(confusion_matrix) {
#     
#     # Convert the entire confusion matrix to a numeric matrix
#     confusion <- as.matrix(confusion_matrix)
#     
#     # Check if the table has non-zero dimensions
#     if (any(dim(confusion_matrix) == 0)) {
#       # If the table is empty, return NA for all metrics
#       return(list(
#         Pos. = 0,TP = 0,FP = 0,Neg. = 0,TN = 0,FN = 0,Sensitivity = 0,Specificity = 0,Accuracy = 0,AUC = 0))
#     }
#     
#     # Ensure the confusion matrix is 2x2 or larger
#     if (all(colnames(confusion) == '0') && all(rownames(confusion) == '0')) {
#       # If not, initialize TP, FN, FP, TN to zero
#       TP <- 0; FN <- 0; FP <- 0; TN <- confusion[1]
#     } else if (all(colnames(confusion) == '1') && all(rownames(confusion) == '1')) {
#       # If not, initialize TP, FN, FP, TN to zero
#       TP <- confusion[1]; FN <- 0; FP <- 0; TN <- 0
#     } else if (all(colnames(confusion) == '0') && all(rownames(confusion) == '1')) {
#       # If not, initialize TP, FN, FP, TN to zero
#       TP <- 0; FN <- confusion[1]; FP <- 0; TN <- 0
#     } else if (all(colnames(confusion) == '1') && all(rownames(confusion) == '0')) {
#       # If not, initialize TP, FN, FP, TN to zero
#       TP <- 0; FN <- 0; FP <- confusion[1]; TN <- 0
#     } else if (all(colnames(confusion) %in% c('0', '1')) && all(rownames(confusion) == '0')) {
#       # If not, initialize TP, FN, FP, TN to zero
#       TP <- 0; FN <- 0; FP <- confusion[2]; TN <- confusion[1]
#     } else if (all(colnames(confusion) %in% c('0', '1')) && all(rownames(confusion) == '1')) {
#       # If not, initialize TP, FN, FP, TN to zero
#       TP <- confusion[2]; FN <- confusion[1]; FP <- 0; TN <- 0
#     } else if (all(colnames(confusion) == '0') && all(rownames(confusion) %in% c('0', '1'))) {
#       # If not, initialize TP, FN, FP, TN to zero
#       TP <- 0; FN <- confusion[3]; FP <- 0; TN <- confusion[1]
#     } else if (all(colnames(confusion) == '1') && all(rownames(confusion) %in% c('0', '1'))) {
#       # If not, initialize TP, FN, FP, TN to zero
#       TP <- confusion[3]; FN <- 0; FP <- confusion[1]; TN <- 0
#     } else {
#       # If it's 2x2 or larger, extract values
#       TP <- confusion[4]; FN <- confusion[2]; FP <- confusion[3]; TN <- confusion[1]
#     }
#     
#     # Sensitivity (True Positive Rate)
#     sensitivity <- ifelse((TP + FN) == 0, 0, (TP / (TP + FN) * 100))
#     # Specificity (True Negative Rate)
#     specificity <- ifelse((TN + FP) == 0, 0, (TN / (TN + FP) * 100))
#     # Accuracy
#     accuracy <- ifelse(sum(confusion_matrix) == 0, 0, ((TP + TN) / sum(confusion_matrix) * 100))
#     # AUC (Area Under the ROC Curve)
#     tpr <- sensitivity; fpr <- 1 - specificity; auc <- 0.5 * (1 + (tpr - fpr))
#     
#     # Return metrics as a named list
#     metrics <- list(
#       Pos. = (TP + FP),TP = TP,FP = FP,Neg. = (TN + FN),TN = TN,FN = FN,Sensitivity = sensitivity,Specificity = specificity,Accuracy = accuracy,AUC = auc)
#     return(metrics)
#   }
#   
#   # Get the unique values in 'signal_index' column
#   signal_values <- data$signal_index
#   result_true_list <- list()
#   
#   ### Comparisons
#   # Create a 2x2 table for each comparison
#   table_fabia <- as.matrix(table(data$signal_index, data$fabia_index))
#   table_mofa <- as.matrix(table(data$signal_index, data$mofa_index))
#   table_mfa <- as.matrix(table(data$signal_index, data$mfa_index))
#   table_gfa <- as.matrix(table(data$signal_index, data$gfa_index))
#   
#   # Calculate metrics for each comparison
#   metrics_fabia <- calculate_performance_metrics(table_fabia)
#   metrics_mofa <- calculate_performance_metrics(table_mofa)
#   metrics_mfa <- calculate_performance_metrics(table_mfa)
#   metrics_gfa <- calculate_performance_metrics(table_gfa)
#   
#   # Add results to the list
#   result_true_list[['FABIA']] <- data.frame(
#     metrics_fabia,stringsAsFactors = FALSE)
#   
#   result_true_list[['MOFA']] <- data.frame(
#     metrics_mofa,stringsAsFactors = FALSE)
#   
#   result_true_list[['MFA']] <- data.frame(
#     metrics_mfa, stringsAsFactors = FALSE)
#   
#   result_true_list[['GFA']] <- data.frame(
#     metrics_mfa, stringsAsFactors = FALSE)
#   
#   # Combine results into a single-row data frame
#   per_measures_true <- do.call(cbind, result_true_list)
#   
#   ### (b. METHODS Comparison)
#   
#   result_ft_list <- list()
#   
#   # Comparisons
#   # Create a 2x2 table for each comparison
#   table_fabia_mofa <- as.matrix(table(data$fabia_index, data$mofa_index))
#   table_fabia_mfa <- as.matrix(table(data$fabia_index, data$mfa_index))
#   table_fabia_gfa <- as.matrix(table(data$fabia_index, data$gfa_index))
#   table_mofa_mfa <- as.matrix(table(data$mofa_index, data$mfa_index))
#   table_mofa_gfa <- as.matrix(table(data$mofa_index, data$gfa_index))
#   table_mfa_gfa <- as.matrix(table(data$mfa_index, data$gfa_index))
#   
#   # Calculate metrics for each comparison
#   metrics_fabia_mofa <- calculate_performance_metrics(table_fabia_mofa)
#   metrics_fabia_mfa <- calculate_performance_metrics(table_fabia_mfa)
#   metrics_fabia_gfa <- calculate_performance_metrics(table_fabia_gfa)
#   metrics_mofa_mfa <- calculate_performance_metrics(table_mofa_mfa)
#   metrics_mofa_gfa <- calculate_performance_metrics(table_mofa_gfa)
#   metrics_mfa_gfa <- calculate_performance_metrics(table_mfa_gfa)
#   
#   # Add results to the list
#   result_ft_list[['FABIA_vs_MOFA']] <- data.frame(
#     metrics_fabia_mofa,
#     stringsAsFactors = FALSE)
#   
#   result_ft_list[['FABIA_vs_MFA']] <- data.frame(
#     metrics_fabia_mfa,
#     stringsAsFactors = FALSE)
#   
#   result_ft_list[['FABIA_vs_GFA']] <- data.frame(
#     metrics_fabia_gfa,
#     stringsAsFactors = FALSE)
#   
#   result_ft_list[['MOFA_vs_MFA']] <- data.frame(
#     metrics_mofa_mfa,
#     stringsAsFactors = FALSE)
#   
#   result_ft_list[['MOFA_vs_GFA']] <- data.frame(
#     metrics_mofa_gfa,
#     stringsAsFactors = FALSE)
#   
#   result_ft_list[['MFA_vs_GFA']] <- data.frame(
#     metrics_mfa_gfa,
#     stringsAsFactors = FALSE)
#   
#   # Combine results into a single-row data frame
#   per_measures_comparison <- do.call(cbind, result_ft_list)
#   
#   # Store results
#   return(list(per_measures_true=per_measures_true, per_measures_comparison = per_measures_comparison))
#   
# }

compare_features <- function(data) {
  
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
    TP <- cm["1", "1"] %>% ifelse(is.na(.), 0, .)
    TN <- cm["0", "0"] %>% ifelse(is.na(.), 0, .)
    FP <- cm["1", "0"] %>% ifelse(is.na(.), 0, .)
    FN <- cm["0", "1"] %>% ifelse(is.na(.), 0, .)
    
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
  
  # Methods to compare
  methods <- c("fabia_index", "mofa_index", "mfa_index", "gfa_index")
  
  ### TRUE Comparisons: Compare methods to the ground truth (signal_index)
  for (method in methods) {
    result_true_list[[method]] <- calculate_metrics(data$signal_index, data[[method]])
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

# Define the performance metrics: samples
# compare_sample <- function(data) {
#   calculate_performance_metrics <- function(confusion_matrix) {
#     
#     # Convert the entire confusion matrix to a numeric matrix
#     confusion <- as.matrix(confusion_matrix)
#     
#     # Check if the table has non-zero dimensions
#     if (any(dim(confusion_matrix) == 0)) {
#       # If the table is empty, return NA for all metrics
#       return(list(
#         Pos. = 0,TP = 0,FP = 0,Neg. = 0,TN = 0,FN = 0,Sensitivity = 0,Specificity = 0,Accuracy = 0,AUC = 0))
#     }
#     
#     # Ensure the confusion matrix is 2x2 or larger
#     if (all(colnames(confusion) == '0') && all(rownames(confusion) == '0')) {
#       # If not, initialize TP, FN, FP, TN to zero
#       TP <- 0; FN <- 0; FP <- 0; TN <- confusion[1]
#     } else if (all(colnames(confusion) == '1') && all(rownames(confusion) == '1')) {
#       # If not, initialize TP, FN, FP, TN to zero
#       TP <- confusion[1]; FN <- 0; FP <- 0; TN <- 0
#     } else if (all(colnames(confusion) == '0') && all(rownames(confusion) == '1')) {
#       # If not, initialize TP, FN, FP, TN to zero
#       TP <- 0; FN <- confusion[1]; FP <- 0; TN <- 0
#     } else if (all(colnames(confusion) == '1') && all(rownames(confusion) == '0')) {
#       # If not, initialize TP, FN, FP, TN to zero
#       TP <- 0; FN <- 0; FP <- confusion[1]; TN <- 0
#     } else if (all(colnames(confusion) %in% c('0', '1')) && all(rownames(confusion) == '0')) {
#       # If not, initialize TP, FN, FP, TN to zero
#       TP <- 0; FN <- 0; FP <- confusion[2]; TN <- confusion[1]
#     } else if (all(colnames(confusion) %in% c('0', '1')) && all(rownames(confusion) == '1')) {
#       # If not, initialize TP, FN, FP, TN to zero
#       TP <- confusion[2]; FN <- confusion[1]; FP <- 0; TN <- 0
#     } else if (all(colnames(confusion) == '0') && all(rownames(confusion) %in% c('0', '1'))) {
#       # If not, initialize TP, FN, FP, TN to zero
#       TP <- 0; FN <- confusion[3]; FP <- 0; TN <- confusion[1]
#     } else if (all(colnames(confusion) == '1') && all(rownames(confusion) %in% c('0', '1'))) {
#       # If not, initialize TP, FN, FP, TN to zero
#       TP <- confusion[3]; FN <- 0; FP <- confusion[1]; TN <- 0
#     } else {
#       # If it's 2x2 or larger, extract values
#       TP <- confusion[4]; FN <- confusion[2]; FP <- confusion[3]; TN <- confusion[1]
#     }
#     
#     # Sensitivity (True Positive Rate)
#     sensitivity <- ifelse((TP + FN) == 0, 0, (TP / (TP + FN) * 100))
#     # Specificity (True Negative Rate)
#     specificity <- ifelse((TN + FP) == 0, 0, (TN / (TN + FP) * 100))
#     # Accuracy
#     accuracy <- ifelse(sum(confusion_matrix) == 0, 0, ((TP + TN) / sum(confusion_matrix) * 100))
#     # AUC (Area Under the ROC Curve)
#     tpr <- sensitivity; fpr <- 1 - specificity; auc <- 0.5 * (1 + (tpr - fpr))
#     
#     # Return metrics as a named list
#     metrics <- list(
#       Pos. = (TP + FP),TP = TP,FP = FP,Neg. = (TN + FN),TN = TN,FN = FN,Sensitivity = sensitivity,Specificity = specificity,Accuracy = accuracy,AUC = auc)
#     return(metrics)
#   }
#   
#   # Get the unique values in 'signal_index' column
#   signal_values <- data$signal_index
#   result_true_list <- list()
#   
#   ### (a. TRUE Comparisons)
#   
#   # Create a 2x2 table for each comparison
#   table_fabia <- as.matrix(table(data$signal_index, data$fabia_index))
#   table_mofa <- as.matrix(table(data$signal_index, data$mofa_index))
#   table_mfa <- as.matrix(table(data$signal_index, data$mfa_index))
#   table_gfa <- as.matrix(table(data$signal_index, data$gfa_index))
#   
#   # Calculate metrics for each comparison
#   metrics_fabia <- calculate_performance_metrics(table_fabia)
#   metrics_mofa <- calculate_performance_metrics(table_mofa)
#   metrics_mfa <- calculate_performance_metrics(table_mfa)
#   metrics_gfa <- calculate_performance_metrics(table_gfa)
#   
#   # Add results to the list
#   result_true_list[['FABIA']] <- data.frame(
#     metrics_fabia,stringsAsFactors = FALSE)
#   result_true_list[['MOFA']] <- data.frame(
#     metrics_mofa,stringsAsFactors = FALSE)
#   result_true_list[['MFA']] <- data.frame(
#     metrics_mfa, stringsAsFactors = FALSE)
#   result_true_list[['GFA']] <- data.frame(
#     metrics_mfa, stringsAsFactors = FALSE)
#   
#   # Combine results into a single-row data frame
#   per_measures_true_smp <- do.call(cbind, result_true_list)
#   
#   ### (b. METHODS Comparisons)
#   
#   result_smp_list <- list()
#   
#   # Create a 2x2 table for each comparison
#   table_fabia_mofa <- as.matrix(table(data$fabia_index, data$mofa_index))
#   table_fabia_mfa <- as.matrix(table(data$fabia_index, data$mfa_index))
#   table_fabia_gfa <- as.matrix(table(data$fabia_index, data$gfa_index))
#   table_mofa_mfa <- as.matrix(table(data$mofa_index, data$mfa_index))
#   table_mofa_gfa <- as.matrix(table(data$mofa_index, data$gfa_index))
#   table_mfa_gfa <- as.matrix(table(data$mfa_index, data$gfa_index))
#   
#   # Calculate metrics for each comparison
#   metrics_fabia_mofa <- calculate_performance_metrics(table_fabia_mofa)
#   metrics_fabia_mfa <- calculate_performance_metrics(table_fabia_mfa)
#   metrics_fabia_gfa <- calculate_performance_metrics(table_fabia_gfa)
#   metrics_mofa_mfa <- calculate_performance_metrics(table_mofa_mfa)
#   metrics_mofa_gfa <- calculate_performance_metrics(table_mofa_gfa)
#   metrics_mfa_gfa <- calculate_performance_metrics(table_mfa_gfa)
#   
#   # Add results to the list
#   result_smp_list[['FABIA_vs_MOFA']] <- data.frame(
#     metrics_fabia_mofa,
#     stringsAsFactors = FALSE)
#   
#   result_smp_list[['FABIA_vs_MFA']] <- data.frame(
#     metrics_fabia_mfa,
#     stringsAsFactors = FALSE)
#   
#   result_smp_list[['FABIA_vs_GFA']] <- data.frame(
#     metrics_fabia_gfa,
#     stringsAsFactors = FALSE)
#   
#   result_smp_list[['MOFA_vs_MFA']] <- data.frame(
#     metrics_mofa_mfa,
#     stringsAsFactors = FALSE)
#   
#   result_smp_list[['MOFA_vs_GFA']] <- data.frame(
#     metrics_mofa_gfa,
#     stringsAsFactors = FALSE)
#   
#   result_smp_list[['MFA_vs_GFA']] <- data.frame(
#     metrics_mfa_gfa,
#     stringsAsFactors = FALSE)
#   
#   # Combine results into a single-row data frame
#   per_measures_comparison_smp <- do.call(cbind, result_smp_list)
#   
#   # Store results
#   return(list(per_measures_true_smp = per_measures_true_smp, per_measures_comparison_smp = per_measures_comparison_smp))
# }

compare_sample <- function(data) {
  
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
    TP <- cm["1", "1"] %>% ifelse(is.na(.), 0, .)
    TN <- cm["0", "0"] %>% ifelse(is.na(.), 0, .)
    FP <- cm["1", "0"] %>% ifelse(is.na(.), 0, .)
    FN <- cm["0", "1"] %>% ifelse(is.na(.), 0, .)
    
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
  
  # Methods to compare
  methods <- c("fabia_index", "mofa_index", "mfa_index", "gfa_index")
  
  ### TRUE Comparisons: Compare methods to the ground truth (signal_index)
  for (method in methods) {
    result_true_list[[method]] <- calculate_metrics(data$signal_index, data[[method]])
  }
  
  ### Pairwise Comparisons Between Methods
  for (i in 1:(length(methods) - 1)) {
    for (j in (i + 1):length(methods)) {
      name <- paste(methods[i], "vs", methods[j], sep = "_")
      result_comparison_list[[name]] <- calculate_metrics(data[[methods[i]]], data[[methods[j]]])
    }
  }
  
  # Combine results into separate data frames
  per_measures_true_smp <- do.call(rbind, lapply(names(result_true_list), function(x) {
    cbind(Method = x, result_true_list[[x]])
  }))
  
  per_measures_comparison_smp <- do.call(rbind, lapply(names(result_comparison_list), function(x) {
    cbind(Method_Comparison = x, result_comparison_list[[x]])
  }))
  
  # Return results as a list of two datasets
  return(list(
    per_measures_true_smp = per_measures_true_smp,
    per_measures_comparison_smp = per_measures_comparison_smp
  ))
}


## Function to re-scale a vector to (-1,1) scale
# rescale_minusone_to_one <- function(vector){
#   scaled_vector <- vector / max(abs(vector), 1e-6) # Using 1e-6 to avoid division by zero
#   return(scaled_vector)
# }
# ## Function to re-scale a vector to (0,1) scale
# scale_to_01 <- function(vector) {
#   abs_vector <- abs(vector) # added code for absolute on '04_04_2024'; abs_vector is the first step
#   min_val <- min(abs_vector)
#   max_val <- max(abs_vector)
#   scaled_vector <- (abs_vector - min_val) / (max_val - min_val)
#   return(scaled_vector)
# }

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
varphi.five <- function(loading_or_score, omic = NULL, type = c("loading", "score")) {
  type <- match.arg(type)
  numeric_values <- as.numeric(unlist(loading_or_score))
  
  if (type == "loading") {
    if (method == 'single.factor') {
      varphi <- quantile(abs(numeric_values), probs = (1 - (length(indices_features.1) / n_features_two)))
    } else if (method %in% c('shared.factor', 'random.shared.factor')) {
      if (omic == 1) {
        x <- (1 - (length(indices_features.1) / n_features_one))
        varphi <- quantile(numeric_values, probs = x)[[1]]
      } else if (omic == 2) {
        y <- (1 - (length(indices_features.2) / n_features_two))
        varphi <- quantile(numeric_values, probs = y)[[1]]
      } else if (omic == 'all') {
        z <- (1 - (length(union(simulation_result$indices_features.1, simulation_result$indices_features.2)) /
                     (n_features_one + n_features_two)))
        varphi <- quantile(numeric_values, probs = z)[[1]]
      }
    } else {
      stop("Invalid method specified.")
    }
    return(varphi)
  }
  
  if (type == "score") {
    if (method == 'single.factor') {
      varphi <- quantile(abs(numeric_values), probs = (1 - (length(indices_samples.1) / n_samples)))
    } else if (method %in% c('shared.factor', 'random.shared.factor')) {
      y <- (1 - (length(union(indices_samples.1, indices_samples.2)) / n_samples))
      varphi <- quantile(numeric_values, probs = y)[[1]]
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

# Varphi functions
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
# -------------------------- MAIN: Define the main simulation function -----------------------------


global_env <- new.env() # shared_env -   # Create a specific environment to hold the variables

# Define the main simulation function
sim_OmicR <- function(n_features_one, n_features_two, n_samples, var_sigma, num_biclusters = 1, num_iterations, method) {
  
  global_env$n_features_one <- n_features_one   # Assign parameters to global environment
  global_env$n_features_two <- n_features_two  # Assign parameters to global environment
  
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
          varphi <- quantile(abs(numeric_values), probs = (1 - (length(indices_features.1) / n_features_two)))
        } else if (method %in% c('shared.factor', 'random.shared.factor')) {
          if (omic == 1) {
            x <- (1 - (length(indices_features.1) / n_features_one))
            varphi <- quantile(numeric_values, probs = x)[[1]]
          } else if (omic == 2) {
            y <- (1 - (length(indices_features.2) / n_features_two))
            varphi <- quantile(numeric_values, probs = y)[[1]]
          } else if (omic == 'all') {
            z <- (1 - (length(union(simulation_result$indices_features.1, simulation_result$indices_features.2)) /
                         (n_features_one + n_features_two)))
            varphi <- quantile(numeric_values, probs = z)[[1]]
          }
        } else {
          stop("Invalid method specified.")
        }
        return(varphi)
      }
      
      if (type == "score") {
        if (method == 'single.factor') {
          varphi <- quantile(abs(numeric_values), probs = (1 - (length(indices_samples.1) / n_samples)))
        } else if (method %in% c('shared.factor', 'random.shared.factor')) {
          y <- (1 - (length(union(indices_samples.1, indices_samples.2)) / n_samples))
          varphi <- quantile(numeric_values, probs = y)[[1]]
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
    
    # Check the method and simulate data
    ## --------------------- Single Factor -------------------------------
    if (method == "single.factor") {
      #current_seed <- sample(1e6, 1)
      #set.seed(current_seed)
      
      # Simulate data
      simulated_data <- single.factor(
        n_features_one = n_features_one,
        n_features_two = n_features_two,
        n_samples = n_samples,
        sigmas = sigma,
        iterations = iter
      )
      
      for (i in 1:length(simulated_data)) {
        try({
          current_data <- simulated_data[[paste0("iteration_", i)]]$concatenated_datasets[[1]]
          indices_features.1 <- simulated_data[[paste0("iteration_", i)]]$indices_features.1
          indices_features.2 <- simulated_data[[paste0("iteration_", i)]]$indices_features.2
          indices_samples.1 <- simulated_data[[paste0("iteration_", i)]]$indices_samples.1
          indices_samples.2 <- simulated_data[[paste0("iteration_", i)]]$indices_samples.2
          
          print(paste("Processing iteration:", i))
          # Run factorization methods
          fabia_result <- func_fabia(data = current_data, num = num_biclusters)
          mofa_result <- func_mofa(data = current_data, num = num_biclusters)
          mfa_result <- func_mfa(data = current_data, num = num_biclusters)
          gfa_result <- func_gfa(data = current_data, num = num_biclusters)
          
          random.data <- data.frame(t(current_data)); random.data$feature = rownames(random.data) 
          random.data$ID <- ifelse(grepl("omic1_", random.data$feature),
                                   sub(".*omic1_feature_(\\d+)", "\\1", random.data$feature),
                                   ifelse(grepl("omic2_", random.data$feature),
                                          sub(".*omic2_feature_(\\d+)", "\\1", random.data$feature),NA))
          
          # Create dataview i.e that will help separate the two merged dataset
          random.data <- random.data %>%
            mutate(dataview = ifelse(grepl("omic1_", feature), "omic.one", ifelse(grepl("omic2_", feature), "omic.two", NA)))
          
          # Filtering the dataset
          random.data.a <- random.data %>% filter(dataview == "omic.one")
          random.data.b <- random.data %>% filter(dataview == "omic.two")
          
          # Assigning 'indices_features.1' to 'random.data.a'
          indices_features.1 = simulated_data[[paste0("iteration_", i)]]$indices_features.1
          in_range.a <- random.data.a$ID %in% indices_features.1
          simulated_features_a <- data.frame(feature = random.data.a$feature, signal = in_range.a)
          
          # Assigning 'indices_features.1' to 'random.data.a'
          indices_features.2 = simulated_data[[paste0("iteration_", i)]]$indices_features.2
          in_range.b <- rep(FALSE, n_features_two)
          simulated_features_b <- data.frame(feature = random.data.b$feature, signal = in_range.b)
          
          # Merge 'simulated_features_a' and 'simulated_features_b'
          simulated_features_sigma = rbind(simulated_features_a, simulated_features_b)
          
          # samples
          indices_samples1 <- paste("sample_", simulated_data[[paste0("iteration_", i)]]$indices_samples.1, sep = "")
          indices_samples2 <- paste("sample_", simulated_data[[paste0("iteration_", i)]]$indices_samples.2, sep = "")
          indices_samples <- unique(c(indices_samples1, indices_samples2))
          in_range_sample <- rownames(current_data) %in% rownames(current_data[indices_samples,])
          simulated_samples_sigma <- data.frame(sample = rownames(current_data),
                                                signal = in_range_sample)
          
          # Extract loadings and scores
          loading_data.not.scaled <- data.frame(
            feature = mofa_result$mofa_loading$feature,
            loading_FABIA = fabia_result$fab_loading$loading_FABIA,
            loading_MOFA = mofa_result$mofa_loading$loading_MOFA,
            loading_MFA = mfa_result$mfa_loading$loading_MFA,
            loading_GFA = gfa_result$gfa_loading$loading_GFA.V
          )
          
          score_data.not.scaled <- data.frame(
            sample = mofa_result$mofa_score$sample,
            score_FABIA = fabia_result$fab_score$score_FABIA,
            score_MOFA = mofa_result$mofa_score$score_MOFA,
            score_MFA = mfa_result$mfa_score$score_MFA,
            score_GFA = gfa_result$gfa_score$score_GFA.V
          )
          
          # Check if mean of FABIA loadings and scores is zero
          mean_fabia_loading <- mean(fabia_result$fab_loading$loading_FABIA, na.rm = TRUE)
          mean_fabia_score <- mean(fabia_result$fab_score$score_FABIA, na.rm = TRUE)
          
          #Convert all numeric columns to absolute values
          score_data.abs <- score_data.not.scaled %>% mutate(across(where(is.numeric), ~abs(.)))
          loading_data.abs <- loading_data.not.scaled %>% mutate(across(where(is.numeric), ~abs(.)))
          
          loading_data <- loading_data.abs #normalize_columns(loading_data.abs, columns_to_loadings_normalize)
          score_data <- score_data.abs #normalize_columns(score_data.abs, columns_to_scores_normalize)
          
          Data_Ft <- merge(x=loading_data,y=simulated_features_sigma, by="feature", all = TRUE)
          #Data_Ft <- data.frame(Data_Ft)[, c('feature','loading_FABIA','loading_MOFA','loading_MFA','signal','loading_GFA')]
          
          # Mutate the dataview column based on the pattern in the feature column
          Data_Ft <- Data_Ft %>%
            mutate(dataview = ifelse(grepl("omic1_", feature), "omic.one", ifelse(grepl("omic2_", feature), "omic.two", NA)))
          
          Data_Ft$ID <- ifelse(grepl("omic1_", Data_Ft$feature),
                               sub(".*omic1_feature_(\\d+)", "\\1", Data_Ft$feature),
                               ifelse(grepl("omic2_", Data_Ft$feature),
                                      sub(".*omic2_feature_(\\d+)", "\\1", Data_Ft$feature),NA))
          
          Data_Smp <- merge(x=score_data,y=simulated_samples_sigma, by="sample", all = TRUE)
          #Data_Smp <- data.frame(Data_Smp)[, c('sample', 'score_FABIA','score_MOFA','score_MFA','signal','score_GFA')]
          
          # Subset omic1
          Data_Ft.a <- Data_Ft %>%
            filter(dataview == "omic.one")
          
          # Subset omic2
          Data_Ft.b <- Data_Ft %>%
            filter(dataview == "omic.two")
          
          for (varphi in seq_along(varphi_functions)) {
            
            varphi_function <- names(varphi_functions)[[varphi]]
            
            # Apply the function to Data_Ft.a
            if (varphi_function == "varphi.five"){
              # Apply the function to Data_Ft.a
              Data_Ft.a.one <- Data_Ft.a %>%
                mutate(
                  fabia_index = ifelse(
                    abs(loading_FABIA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_FABIA,
                      omic = 1
                    )), 1, 0),
                  
                  mofa_index = ifelse(
                    abs(loading_MOFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_MOFA,
                      omic = 1
                    )), 1, 0),
                  
                  mfa_index = ifelse(
                    abs(loading_MFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_MFA,
                      omic = 1
                    )), 1, 0),
                  
                  gfa_index = ifelse(
                    abs(loading_GFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_GFA,
                      omic = 1
                    )), 1, 0),
                  
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
              
            }else{
              # Apply the function to Data_Ft.a
              Data_Ft.a.one <- Data_Ft.a %>%
                mutate(
                  fabia_index = ifelse(abs(loading_FABIA) >= get(varphi_function)(loading_FABIA), 1, 0),
                  mofa_index = ifelse(abs(loading_MOFA) >= get(varphi_function)(loading_MOFA), 1, 0),
                  mfa_index = ifelse(abs(loading_MFA) >= get(varphi_function)(loading_MFA), 1, 0),
                  gfa_index = ifelse(abs(loading_GFA) >= get(varphi_function)(loading_GFA), 1, 0),
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
            }
            
            if (varphi_function == "varphi.five"){
              Data_Ft.b.two <- Data_Ft.b %>%
                mutate(
                  fabia_index = ifelse(
                    abs(loading_FABIA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_FABIA,
                      omic = 2
                    )), 1, 0),
                  
                  mofa_index = ifelse(
                    abs(loading_MOFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_MOFA,
                      omic = 2
                    )), 1, 0),
                  
                  mfa_index = ifelse(
                    abs(loading_MFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_MFA,
                      omic = 2
                    )), 1, 0),
                  
                  gfa_index = ifelse(
                    abs(loading_GFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_GFA,
                      omic = 2
                    )), 1, 0),
                  
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
            } else{
              Data_Ft.b.two <- Data_Ft.b %>%
                mutate(
                  fabia_index = ifelse(abs(loading_FABIA) >= get(varphi_function)(loading_FABIA), 1, 0),
                  mofa_index = ifelse(abs(loading_MOFA) >= get(varphi_function)(loading_MOFA), 1, 0),
                  mfa_index = ifelse(abs(loading_MFA) >= get(varphi_function)(loading_MFA), 1, 0),
                  gfa_index = ifelse(abs(loading_GFA) >= get(varphi_function)(loading_GFA), 1, 0),
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
            }
            
            if (varphi_function == "varphi.five"){
              Data_Ft2 <- Data_Ft %>%
                mutate(
                  fabia_index = ifelse(
                    abs(loading_FABIA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_FABIA,
                      omic = "ALL"
                    )), 1, 0),
                  
                  mofa_index = ifelse(
                    abs(loading_MOFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_MOFA,
                      omic = "ALL"
                    )), 1, 0),
                  
                  mfa_index = ifelse(
                    abs(loading_MFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_MFA,
                      omic = "ALL"
                    )), 1, 0),
                  
                  gfa_index = ifelse(
                    abs(loading_GFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_GFA,
                      omic = "ALL"
                    )), 1, 0),
                  
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
            } else{
              Data_Ft2 <- Data_Ft %>%
                mutate(
                  fabia_index = ifelse(abs(loading_FABIA) >= get(varphi_function)(loading_FABIA), 1, 0),
                  mofa_index = ifelse(abs(loading_MOFA) >= get(varphi_function)(loading_MOFA), 1, 0),
                  mfa_index = ifelse(abs(loading_MFA) >= get(varphi_function)(loading_MFA), 1, 0),
                  gfa_index = ifelse(abs(loading_GFA) >= get(varphi_function)(loading_GFA), 1, 0),
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
            }
            
            if (varphi_function == "varphi.five"){
              Data_Smp2 <- Data_Smp %>%
                mutate(
                  fabia_index = ifelse(
                    abs(score_FABIA) >= do.call(get(varphi_function), list(
                      loading_or_score = score_FABIA
                    )), 1, 0),
                  
                  mofa_index = ifelse(
                    abs(score_MOFA) >= do.call(get(varphi_function), list(
                      loading_or_score = score_MOFA
                    )), 1, 0),
                  
                  mfa_index = ifelse(
                    abs(score_MFA) >= do.call(get(varphi_function), list(
                      loading_or_score = score_MFA
                    )), 1, 0),
                  
                  gfa_index = ifelse(
                    abs(score_GFA) >= do.call(get(varphi_function), list(
                      loading_or_score = score_GFA
                    )), 1, 0),
                  
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
              
            } else{
              Data_Smp2 <- Data_Smp %>%
                mutate(
                  fabia_index = ifelse(abs(score_FABIA) >= get(varphi_function)(score_FABIA), 1, 0),
                  mofa_index = ifelse(abs(score_MOFA) >= get(varphi_function)(score_MOFA), 1, 0),
                  mfa_index = ifelse(abs(score_MFA) >= get(varphi_function)(score_MFA), 1, 0),
                  gfa_index = ifelse(abs(score_GFA) >= get(varphi_function)(score_GFA), 1, 0),
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
            }
            
            # Create binary indices
            Data_Matrix_Ft <- Data_Ft2
            
            # Create binary indices
            Data_Matrix_Smp <- Data_Smp2
            
            # Signal MATRIX
            Data_Matrix_Smp <- Data_Matrix_Smp[order(as.numeric(gsub("sample_", "", Data_Matrix_Smp$sample))), ]
            
            # Sort Data_Matrix_Ft
            Data_Matrix_Ft <- Data_Matrix_Ft[order(
              # Step 1: Sort by prefix ('omic1' first, 'omic2' second)
              grepl("^omic2", Data_Matrix_Ft$feature),  
              # Step 2: Sort by the numeric part of 'feature'
              as.numeric(gsub(".*_(\\d+)$", "\\1", Data_Matrix_Ft$feature))
            ), ]
            
            sample_signal_index = Data_Matrix_Smp$signal_index
            feature_signal_index = Data_Matrix_Ft$signal_index
            signal_matrix = as.matrix(feature_signal_index)%*%as.matrix(t(sample_signal_index))
            
            # loading_score_matrix
            fabia_matrix = as.matrix(Data_Matrix_Ft$fabia_index)%*%as.matrix(t(Data_Matrix_Smp$fabia_index))
            mofa_matrix = as.matrix(Data_Matrix_Ft$mofa_index)%*%as.matrix(t(Data_Matrix_Smp$mofa_index))
            mfa_matrix = as.matrix(Data_Matrix_Ft$mfa_index)%*%as.matrix(t(Data_Matrix_Smp$mfa_index))
            gfa_matrix = as.matrix(Data_Matrix_Ft$gfa_index)%*%as.matrix(t(Data_Matrix_Smp$gfa_index))
            
            fabia_raw_matrix = as.matrix(Data_Matrix_Ft$loading_FABIA)%*%as.matrix(t(Data_Matrix_Smp$score_FABIA))
            mofa_raw_matrix = as.matrix(Data_Matrix_Ft$loading_MOFA)%*%as.matrix(t(Data_Matrix_Smp$score_MOFA))
            mfa_raw_matrix = as.matrix(Data_Matrix_Ft$loading_MFA)%*%as.matrix(t(Data_Matrix_Smp$score_MFA))
            gfa_raw_matrix = as.matrix(Data_Matrix_Ft$loading_GFA)%*%as.matrix(t(Data_Matrix_Smp$score_GFA))
            
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ JACCARD INDEX ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            # List of methods and corresponding indices
            methods <- c("fabia", "mofa", "mfa", "gfa")
            
            # Compute Jaccard indices in a single step
            ji_true <- sapply(methods, function(method) {
              method_index <- Data_Ft2[[paste0(method, "_index")]]
              jaccard.index.sim(method_index, Data_Ft2$signal_index)
            }, simplify = TRUE, USE.NAMES = TRUE)
            
            ji_omic_one <- sapply(methods, function(method) {
              method_index <- Data_Ft.a.one[[paste0(method, "_index")]]
              jaccard.index.sim(method_index, Data_Ft.a.one$signal_index)
            }, simplify = TRUE, USE.NAMES = TRUE)
            
            ji_omic_two <- sapply(methods, function(method) {
              method_index <- Data_Ft.b.two[[paste0(method, "_index")]]
              jaccard.index.sim(method_index, Data_Ft.b.two$signal_index)
            }, simplify = TRUE, USE.NAMES = TRUE)
            
            ji_samples <- sapply(methods, function(method) {
              method_index <- Data_Smp2[[paste0(method, "_index")]]
              jaccard.index.sim(method_index, Data_Smp2$signal_index)
            }, simplify = TRUE, USE.NAMES = TRUE)
            
            # Define matrices
            matrices <- list(
              fabia = fabia_matrix,
              mofa = mofa_matrix,
              mfa = mfa_matrix,
              gfa = gfa_matrix
            )
            
            # Compute Jaccard similarity indices
            ji_matrix <- sapply(names(matrices), function(method) {
              jaccard.index.sim(matrices[[method]], signal_matrix)
            }, simplify = TRUE, USE.NAMES = TRUE)
            
            
            # Organize the results for the current method
            jaccard_loading <- list(ji_true=ji_true, ji_omic_one=ji_omic_one, ji_omic_two=ji_omic_two, ji_samples=ji_samples, ji_matrix=ji_matrix)
            
            # Create a unique name for each dataset
            dataset_name <- sprintf("variance_%d_iteration_%d_varphi_%d", sigma, i, varphi) 
            jaccard_results[[dataset_name]] <- jaccard_loading
            
            # COMPARISONS BY METHODS
            methods <- c("fabia", "mofa", "mfa", "gfa")  # Define the list of methods
            
            # ALL DATA
            # Calculate pairwise Jaccard indices dynamically
            ji_true_pairs <- combn(methods, 2, function(method_pair) {
              method1_index <- Data_Ft2[[paste0(method_pair[1], "_index")]]
              method2_index <- Data_Ft2[[paste0(method_pair[2], "_index")]]
              jaccard.index.sim(method1_index, method2_index)
            }, simplify = TRUE)
            
            # Assign names to the results
            names(ji_true_pairs) <- combn(methods, 2, function(method_pair) {
              paste(method_pair[1], method_pair[2], "ji_phi", sep = "_")
            })
            
            # Convert to a named vector or data frame for further use
            ji_true_pairs <- as.data.frame(as.list(ji_true_pairs))
            
            # OMIC ONE
            # Calculate pairwise Jaccard indices dynamically
            ji_omic_one_pairs <- combn(methods, 2, function(method_pair) {
              method1_index <- Data_Ft.a.one[[paste0(method_pair[1], "_index")]]
              method2_index <- Data_Ft.a.one[[paste0(method_pair[2], "_index")]]
              jaccard.index.sim(method1_index, method2_index)
            }, simplify = TRUE)
            
            # Assign names to the results
            names(ji_omic_one_pairs) <- combn(methods, 2, function(method_pair) {
              paste(method_pair[1], method_pair[2], "ji_phi", sep = "_")
            })
            
            # Convert to a named vector or data frame for further use
            ji_omic_one_pairs <- as.data.frame(as.list(ji_omic_one_pairs))
            
            # OMIC TWO
            # Calculate pairwise Jaccard indices dynamically
            ji_omic_two_pairs <- combn(methods, 2, function(method_pair) {
              method1_index <- Data_Ft.b.two[[paste0(method_pair[1], "_index")]]
              method2_index <- Data_Ft.b.two[[paste0(method_pair[2], "_index")]]
              jaccard.index.sim(method1_index, method2_index)
            }, simplify = TRUE)
            
            # Assign names to the results
            names(ji_omic_two_pairs) <- combn(methods, 2, function(method_pair) {
              paste(method_pair[1], method_pair[2], "ji_phi", sep = "_")
            })
            
            # Convert to a named vector or data frame for further use
            ji_omic_two_pairs <- as.data.frame(as.list(ji_omic_two_pairs))
            
            # SAMPLES
            # Calculate pairwise Jaccard indices dynamically
            ji_smp_pairs <- combn(methods, 2, function(method_pair) {
              method1_index <- Data_Smp2[[paste0(method_pair[1], "_index")]]
              method2_index <- Data_Smp2[[paste0(method_pair[2], "_index")]]
              jaccard.index.sim(method1_index, method2_index)
            }, simplify = TRUE)
            
            # Assign names to the results
            names(ji_smp_pairs) <- combn(methods, 2, function(method_pair) {
              paste(method_pair[1], method_pair[2], "ji_phi", sep = "_")
            })
            
            # Convert to a named vector or data frame for further use
            ji_smp_pairs <- as.data.frame(as.list(ji_smp_pairs))
            
            # MATRICES
            
            # Define the matrices in a named list
            matrices <- list(
              fabia = fabia_matrix,
              mofa = mofa_matrix,
              mfa = mfa_matrix,
              gfa = gfa_matrix
            )
            
            # Get all unique pair combinations of matrix names
            matrix_pairs <- combn(names(matrices), 2, simplify = FALSE)
            
            # Compute Jaccard similarity indices for each pair
            ji_matrix_pairs <- sapply(matrix_pairs, function(pair) {
              jaccard.index.sim(matrices[[pair[1]]], matrices[[pair[2]]])
            }, simplify = TRUE, USE.NAMES = FALSE)
            
            # Assign names to the results
            names(ji_matrix_pairs) <- sapply(matrix_pairs, function(pair) {
              paste(pair[1], pair[2], sep = "_")
            })
            
            # Organize the results for the current method
            jaccard_comparison <- list(ji_true_pairs=ji_true_pairs, ji_omic_one_pairs=ji_omic_one_pairs, 
                                       ji_omic_two_pairs=ji_omic_two_pairs, ji_smp_pairs=ji_smp_pairs,
                                       ji_matrix_pairs=ji_matrix_pairs)
            
            # Create a unique name for each dataset
            dataset_name <- sprintf("variance_%d_iteration_%d_varphi_%d", sigma, i, varphi) 
            jaccard_comparison_results[[dataset_name]] <- jaccard_comparison
            
            # COMPARISON METRICS: SENSITIVITY, SPECIFICITY; ACCURACY
            per_measures_true <- compare_features(Data_Ft2)
            per_measures_omic_one <- compare_features(Data_Ft.a.one)
            per_measures_omic_two <- compare_features(Data_Ft.b.two)
            per_measures_smp <- compare_sample(Data_Smp2)
            
            # Organize the results for the current method
            per_measures <- list(per_measures_true=per_measures_true, per_measures_omic_one=per_measures_omic_one, per_measures_omic_two=per_measures_omic_two, per_measures_smp=per_measures_smp)
            
            # Create a unique name for each dataset
            dataset_name <- sprintf("variance_%d_iteration_%d_varphi_%d", sigma, i, varphi) 
            per_measures_results[[dataset_name]] <- per_measures
            
          }
          #simulatedData <- simulated_data[[paste0("iteration_", i)]]$concatenated_datasets[[1]]
          # Create a list of datasets
          datasets <- list(
            simulatedData = current_data,
            Data_Ft = Data_Ft,
            Data_Smp = Data_Smp
          ) 
          dataset <- sprintf("variance_%d_iteration_%d", sigma, i)
          dataset_output[[dataset]] <- datasets
          
        }, silent = TRUE)
        
      }
      
      ## --------------------------- Shared Factor -------------------
    } else if (method == "shared.factor") {
      #current_seed <- sample(1e6, 1)
      #set.seed(current_seed)
      
      # Simulate data
      simulated_data <- shared.factor(
        n_features_one = n_features_one,
        n_features_two = n_features_two,
        n_samples = n_samples,
        sigmas = sigma,
        iterations = iter 
      )
      
      for (i in 1:length(simulated_data)) {
        try({
          current_data <- simulated_data[[paste0("iteration_", i)]]$concatenated_datasets[[1]]
          indices_features.1 <- simulated_data[[paste0("iteration_", i)]]$indices_features.1
          indices_features.2 <- simulated_data[[paste0("iteration_", i)]]$indices_features.2
          indices_samples.1 <- simulated_data[[paste0("iteration_", i)]]$indices_samples.1
          indices_samples.2 <- simulated_data[[paste0("iteration_", i)]]$indices_samples.2
          
          print(paste("Processing iteration:", i))
          # Run factorization methods
          fabia_result <- func_fabia(data = current_data, num = num_biclusters)
          mofa_result <- func_mofa(data = current_data, num = num_biclusters)
          mfa_result <- func_mfa(data = current_data, num = num_biclusters)
          gfa_result <- func_gfa(data = current_data, num = num_biclusters)
          
          random.data <- data.frame(t(current_data)); random.data$feature = rownames(random.data) 
          random.data$ID <- ifelse(grepl("omic1_", random.data$feature),
                                   sub(".*omic1_feature_(\\d+)", "\\1", random.data$feature),
                                   ifelse(grepl("omic2_", random.data$feature),
                                          sub(".*omic2_feature_(\\d+)", "\\1", random.data$feature),NA))
          
          # Create dataview i.e that will help separate the two merged dataset
          random.data <- random.data %>%
            mutate(dataview = ifelse(grepl("omic1_", feature), "omic.one", ifelse(grepl("omic2_", feature), "omic.two", NA)))
          
          # Filtering the dataset
          random.data.a <- random.data %>% filter(dataview == "omic.one")
          random.data.b <- random.data %>% filter(dataview == "omic.two")
          
          # Assigning 'indices_features.1' to 'random.data.a'
          indices_features.1 = simulated_data[[paste0("iteration_", i)]]$indices_features.1
          in_range.a <- random.data.a$ID %in% indices_features.1
          simulated_features_a <- data.frame(feature = random.data.a$feature, signal = in_range.a)
          
          # Assigning 'indices_features.1' to 'random.data.a'
          indices_features.2 = simulated_data[[paste0("iteration_", i)]]$indices_features.2
          in_range.b <- random.data.b$ID %in% indices_features.2
          simulated_features_b <- data.frame(feature = random.data.b$feature, signal = in_range.b)
          
          # Merge 'simulated_features_a' and 'simulated_features_b'
          simulated_features_sigma = rbind(simulated_features_a, simulated_features_b)
          
          # samples
          indices_samples1 <- paste("sample_", simulated_data[[paste0("iteration_", i)]]$indices_samples.1, sep = "")
          indices_samples2 <- paste("sample_", simulated_data[[paste0("iteration_", i)]]$indices_samples.2, sep = "")
          indices_samples <- unique(c(indices_samples1, indices_samples2))
          in_range_sample <- rownames(current_data) %in% rownames(current_data[indices_samples,])
          simulated_samples_sigma <- data.frame(sample = rownames(current_data),
                                                signal = in_range_sample)
          
          # Extract loadings and scores
          loading_data.not.scaled <- data.frame(
            feature = mofa_result$mofa_loading$feature,
            loading_FABIA = fabia_result$fab_loading$loading_FABIA,
            loading_MOFA = mofa_result$mofa_loading$loading_MOFA,
            loading_MFA = mfa_result$mfa_loading$loading_MFA,
            loading_GFA = gfa_result$gfa_loading$loading_GFA.V
          )
          
          score_data.not.scaled <- data.frame(
            sample = mofa_result$mofa_score$sample,
            score_FABIA = fabia_result$fab_score$score_FABIA,
            score_MOFA = mofa_result$mofa_score$score_MOFA,
            score_MFA = mfa_result$mfa_score$score_MFA,
            score_GFA = gfa_result$gfa_score$score_GFA.V
          )
          
          # Check if mean of FABIA loadings and scores is zero
          mean_fabia_loading <- mean(fabia_result$fab_loading$loading_FABIA, na.rm = TRUE)
          mean_fabia_score <- mean(fabia_result$fab_score$score_FABIA, na.rm = TRUE)
          
          # Convert data to absolute values
          #Convert all numeric columns to absolute values
          score_data.abs <- score_data.not.scaled %>% mutate(across(where(is.numeric), ~abs(.)))
          loading_data.abs <- loading_data.not.scaled %>% mutate(across(where(is.numeric), ~abs(.)))
          
          loading_data <- loading_data.abs #normalize_columns(loading_data.abs, columns_to_loadings_normalize)
          score_data <- score_data.abs #normalize_columns(score_data.abs, columns_to_scores_normalize)
          
          Data_Ft <- merge(x=loading_data,y=simulated_features_sigma, by="feature", all = TRUE)
          #Data_Ft <- data.frame(Data_Ft)[, c('feature','loading_FABIA','loading_MOFA','loading_MFA','loading_GFA','signal')]
          
          # Mutate the dataview column based on the pattern in the feature column
          Data_Ft <- Data_Ft %>%
            mutate(dataview = ifelse(grepl("omic1_", feature), "omic.one", ifelse(grepl("omic2_", feature), "omic.two", NA)))
          
          Data_Ft$ID <- ifelse(grepl("omic1_", Data_Ft$feature),
                               sub(".*omic1_feature_(\\d+)", "\\1", Data_Ft$feature),
                               ifelse(grepl("omic2_", Data_Ft$feature),
                                      sub(".*omic2_feature_(\\d+)", "\\1", Data_Ft$feature),NA))
          
          Data_Smp <- merge(x=score_data,y=simulated_samples_sigma, by="sample", all = TRUE)
          #Data_Smp <- data.frame(Data_Smp)[, c('sample', 'score_FABIA','score_MOFA','score_MFA','signal','score_GFA')]
          
          # Subset omic1
          Data_Ft.a <- Data_Ft %>%
            filter(dataview == "omic.one")
          
          # Subset omic2
          Data_Ft.b <- Data_Ft %>%
            filter(dataview == "omic.two")
          
          for (varphi in seq_along(varphi_functions)) {
            
            varphi_function <- names(varphi_functions)[[varphi]]
            
            # Apply the function to Data_Ft.a
            if (varphi_function == "varphi.five"){
              # Apply the function to Data_Ft.a
              Data_Ft.a.one <- Data_Ft.a %>%
                mutate(
                  fabia_index = ifelse(
                    abs(loading_FABIA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_FABIA,
                      omic = 1
                    )), 1, 0),
                  
                  mofa_index = ifelse(
                    abs(loading_MOFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_MOFA,
                      omic = 1
                    )), 1, 0),
                  
                  mfa_index = ifelse(
                    abs(loading_MFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_MFA,
                      omic = 1
                    )), 1, 0),
                  
                  gfa_index = ifelse(
                    abs(loading_GFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_GFA,
                      omic = 1
                    )), 1, 0),
                  
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
              
            }else{
              # Apply the function to Data_Ft.a
              Data_Ft.a.one <- Data_Ft.a %>%
                mutate(
                  fabia_index = ifelse(abs(loading_FABIA) >= get(varphi_function)(loading_FABIA), 1, 0),
                  mofa_index = ifelse(abs(loading_MOFA) >= get(varphi_function)(loading_MOFA), 1, 0),
                  mfa_index = ifelse(abs(loading_MFA) >= get(varphi_function)(loading_MFA), 1, 0),
                  gfa_index = ifelse(abs(loading_GFA) >= get(varphi_function)(loading_GFA), 1, 0),
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
            }
            
            if (varphi_function == "varphi.five"){
              Data_Ft.b.two <- Data_Ft.b %>%
                mutate(
                  fabia_index = ifelse(
                    abs(loading_FABIA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_FABIA,
                      omic = 2
                    )), 1, 0),
                  
                  mofa_index = ifelse(
                    abs(loading_MOFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_MOFA,
                      omic = 2
                    )), 1, 0),
                  
                  mfa_index = ifelse(
                    abs(loading_MFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_MFA,
                      omic = 2
                    )), 1, 0),
                  
                  gfa_index = ifelse(
                    abs(loading_GFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_GFA,
                      omic = 2
                    )), 1, 0),
                  
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
            } else{
              Data_Ft.b.two <- Data_Ft.b %>%
                mutate(
                  fabia_index = ifelse(abs(loading_FABIA) >= get(varphi_function)(loading_FABIA), 1, 0),
                  mofa_index = ifelse(abs(loading_MOFA) >= get(varphi_function)(loading_MOFA), 1, 0),
                  mfa_index = ifelse(abs(loading_MFA) >= get(varphi_function)(loading_MFA), 1, 0),
                  gfa_index = ifelse(abs(loading_GFA) >= get(varphi_function)(loading_GFA), 1, 0),
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
            }
            
            if (varphi_function == "varphi.five"){
              Data_Ft2 <- Data_Ft %>%
                mutate(
                  fabia_index = ifelse(
                    abs(loading_FABIA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_FABIA,
                      omic = "ALL"
                    )), 1, 0),
                  
                  mofa_index = ifelse(
                    abs(loading_MOFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_MOFA,
                      omic = "ALL"
                    )), 1, 0),
                  
                  mfa_index = ifelse(
                    abs(loading_MFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_MFA,
                      omic = "ALL"
                    )), 1, 0),
                  
                  gfa_index = ifelse(
                    abs(loading_GFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_GFA,
                      omic = "ALL"
                    )), 1, 0),
                  
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
            } else{
              Data_Ft2 <- Data_Ft %>%
                mutate(
                  fabia_index = ifelse(abs(loading_FABIA) >= get(varphi_function)(loading_FABIA), 1, 0),
                  mofa_index = ifelse(abs(loading_MOFA) >= get(varphi_function)(loading_MOFA), 1, 0),
                  mfa_index = ifelse(abs(loading_MFA) >= get(varphi_function)(loading_MFA), 1, 0),
                  gfa_index = ifelse(abs(loading_GFA) >= get(varphi_function)(loading_GFA), 1, 0),
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
            }
            
            if (varphi_function == "varphi.five"){
              Data_Smp2 <- Data_Smp %>%
                mutate(
                  fabia_index = ifelse(
                    abs(score_FABIA) >= do.call(get(varphi_function), list(
                      loading_or_score = score_FABIA
                    )), 1, 0),
                  
                  mofa_index = ifelse(
                    abs(score_MOFA) >= do.call(get(varphi_function), list(
                      loading_or_score = score_MOFA
                    )), 1, 0),
                  
                  mfa_index = ifelse(
                    abs(score_MFA) >= do.call(get(varphi_function), list(
                      loading_or_score = score_MFA
                    )), 1, 0),
                  
                  gfa_index = ifelse(
                    abs(score_GFA) >= do.call(get(varphi_function), list(
                      loading_or_score = score_GFA
                    )), 1, 0),
                  
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
              
            } else{
              Data_Smp2 <- Data_Smp %>%
                mutate(
                  fabia_index = ifelse(abs(score_FABIA) >= get(varphi_function)(score_FABIA), 1, 0),
                  mofa_index = ifelse(abs(score_MOFA) >= get(varphi_function)(score_MOFA), 1, 0),
                  mfa_index = ifelse(abs(score_MFA) >= get(varphi_function)(score_MFA), 1, 0),
                  gfa_index = ifelse(abs(score_GFA) >= get(varphi_function)(score_GFA), 1, 0),
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
            }
            
            
            # Create binary indices
            Data_Matrix_Ft <- Data_Ft2
            
            # Create binary indices
            Data_Matrix_Smp <- Data_Smp2
            
            # Signal MATRIX
            Data_Matrix_Smp <- Data_Matrix_Smp[order(as.numeric(gsub("sample_", "", Data_Matrix_Smp$sample))), ]
            
            # Sort Data_Matrix_Ft
            Data_Matrix_Ft <- Data_Matrix_Ft[order(
              # Step 1: Sort by prefix ('omic1' first, 'omic2' second)
              grepl("^omic2", Data_Matrix_Ft$feature),  
              # Step 2: Sort by the numeric part of 'feature'
              as.numeric(gsub(".*_(\\d+)$", "\\1", Data_Matrix_Ft$feature))
            ), ]
            
            # Signal MATRIX
            sample_signal_index = Data_Matrix_Smp$signal_index
            feature_signal_index = Data_Matrix_Ft$signal_index
            signal_matrix = as.matrix(feature_signal_index)%*%as.matrix(t(sample_signal_index))
            
            # loading_score_matrix
            fabia_matrix = as.matrix(Data_Matrix_Ft$fabia_index)%*%as.matrix(t(Data_Matrix_Smp$fabia_index))
            mofa_matrix = as.matrix(Data_Matrix_Ft$mofa_index)%*%as.matrix(t(Data_Matrix_Smp$mofa_index))
            mfa_matrix = as.matrix(Data_Matrix_Ft$mfa_index)%*%as.matrix(t(Data_Matrix_Smp$mfa_index))
            gfa_matrix = as.matrix(Data_Matrix_Ft$gfa_index)%*%as.matrix(t(Data_Matrix_Smp$gfa_index))
            
            fabia_raw_matrix = as.matrix(Data_Matrix_Ft$loading_FABIA)%*%as.matrix(t(Data_Matrix_Smp$score_FABIA))
            mofa_raw_matrix = as.matrix(Data_Matrix_Ft$loading_MOFA)%*%as.matrix(t(Data_Matrix_Smp$score_MOFA))
            mfa_raw_matrix = as.matrix(Data_Matrix_Ft$loading_MFA)%*%as.matrix(t(Data_Matrix_Smp$score_MFA))
            gfa_raw_matrix = as.matrix(Data_Matrix_Ft$loading_GFA)%*%as.matrix(t(Data_Matrix_Smp$score_GFA))
            
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ JACCARD INDEX ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            # List of methods and corresponding indices
            methods <- c("fabia", "mofa", "mfa", "gfa")
            
            # Compute Jaccard indices in a single step
            ji_true <- sapply(methods, function(method) {
              method_index <- Data_Ft2[[paste0(method, "_index")]]
              jaccard.index.sim(method_index, Data_Ft2$signal_index)
            }, simplify = TRUE, USE.NAMES = TRUE)
            
            ji_omic_one <- sapply(methods, function(method) {
              method_index <- Data_Ft.a.one[[paste0(method, "_index")]]
              jaccard.index.sim(method_index, Data_Ft.a.one$signal_index)
            }, simplify = TRUE, USE.NAMES = TRUE)
            
            ji_omic_two <- sapply(methods, function(method) {
              method_index <- Data_Ft.b.two[[paste0(method, "_index")]]
              jaccard.index.sim(method_index, Data_Ft.b.two$signal_index)
            }, simplify = TRUE, USE.NAMES = TRUE)
            
            ji_samples <- sapply(methods, function(method) {
              method_index <- Data_Smp2[[paste0(method, "_index")]]
              jaccard.index.sim(method_index, Data_Smp2$signal_index)
            }, simplify = TRUE, USE.NAMES = TRUE)
            
            # Define matrices
            matrices <- list(
              fabia = fabia_matrix,
              mofa = mofa_matrix,
              mfa = mfa_matrix,
              gfa = gfa_matrix
            )
            
            # Compute Jaccard similarity indices
            ji_matrix <- sapply(names(matrices), function(method) {
              jaccard.index.sim(matrices[[method]], signal_matrix)
            }, simplify = TRUE, USE.NAMES = TRUE)
            
            
            # Organize the results for the current method
            jaccard_loading <- list(ji_true=ji_true, ji_omic_one=ji_omic_one, ji_omic_two=ji_omic_two, 
                                    ji_samples=ji_samples, ji_matrix=ji_matrix)
            
            # Create a unique name for each dataset
            dataset_name <- sprintf("variance_%d_iteration_%d_varphi_%d", sigma, i, varphi) 
            jaccard_results[[dataset_name]] <- jaccard_loading
            
            # COMPARISONS BY METHODS
            methods <- c("fabia", "mofa", "mfa", "gfa")  # Define the list of methods
            
            # ALL DATA
            # Calculate pairwise Jaccard indices dynamically
            ji_true_pairs <- combn(methods, 2, function(method_pair) {
              method1_index <- Data_Ft2[[paste0(method_pair[1], "_index")]]
              method2_index <- Data_Ft2[[paste0(method_pair[2], "_index")]]
              jaccard.index.sim(method1_index, method2_index)
            }, simplify = TRUE)
            
            # Assign names to the results
            names(ji_true_pairs) <- combn(methods, 2, function(method_pair) {
              paste(method_pair[1], method_pair[2], "ji_phi", sep = "_")
            })
            
            # Convert to a named vector or data frame for further use
            ji_true_pairs <- as.data.frame(as.list(ji_true_pairs))
            
            # OMIC ONE
            # Calculate pairwise Jaccard indices dynamically
            ji_omic_one_pairs <- combn(methods, 2, function(method_pair) {
              method1_index <- Data_Ft.a.one[[paste0(method_pair[1], "_index")]]
              method2_index <- Data_Ft.a.one[[paste0(method_pair[2], "_index")]]
              jaccard.index.sim(method1_index, method2_index)
            }, simplify = TRUE)
            
            # Assign names to the results
            names(ji_omic_one_pairs) <- combn(methods, 2, function(method_pair) {
              paste(method_pair[1], method_pair[2], "ji_phi", sep = "_")
            })
            
            # Convert to a named vector or data frame for further use
            ji_omic_one_pairs <- as.data.frame(as.list(ji_omic_one_pairs))
            
            # OMIC TWO
            # Calculate pairwise Jaccard indices dynamically
            ji_omic_two_pairs <- combn(methods, 2, function(method_pair) {
              method1_index <- Data_Ft.b.two[[paste0(method_pair[1], "_index")]]
              method2_index <- Data_Ft.b.two[[paste0(method_pair[2], "_index")]]
              jaccard.index.sim(method1_index, method2_index)
            }, simplify = TRUE)
            
            # Assign names to the results
            names(ji_omic_two_pairs) <- combn(methods, 2, function(method_pair) {
              paste(method_pair[1], method_pair[2], "ji_phi", sep = "_")
            })
            
            # Convert to a named vector or data frame for further use
            ji_omic_two_pairs <- as.data.frame(as.list(ji_omic_two_pairs))
            
            # SAMPLES
            # Calculate pairwise Jaccard indices dynamically
            ji_smp_pairs <- combn(methods, 2, function(method_pair) {
              method1_index <- Data_Smp2[[paste0(method_pair[1], "_index")]]
              method2_index <- Data_Smp2[[paste0(method_pair[2], "_index")]]
              jaccard.index.sim(method1_index, method2_index)
            }, simplify = TRUE)
            
            # Assign names to the results
            names(ji_smp_pairs) <- combn(methods, 2, function(method_pair) {
              paste(method_pair[1], method_pair[2], "ji_phi", sep = "_")
            })
            
            # Convert to a named vector or data frame for further use
            ji_smp_pairs <- as.data.frame(as.list(ji_smp_pairs))
            
            # MATRICES
            
            # Define the matrices in a named list
            matrices <- list(
              fabia = fabia_matrix,
              mofa = mofa_matrix,
              mfa = mfa_matrix,
              gfa = gfa_matrix
            )
            
            # Get all unique pair combinations of matrix names
            matrix_pairs <- combn(names(matrices), 2, simplify = FALSE)
            
            # Compute Jaccard similarity indices for each pair
            ji_matrix_pairs <- sapply(matrix_pairs, function(pair) {
              jaccard.index.sim(matrices[[pair[1]]], matrices[[pair[2]]])
            }, simplify = TRUE, USE.NAMES = FALSE)
            
            # Assign names to the results
            names(ji_matrix_pairs) <- sapply(matrix_pairs, function(pair) {
              paste(pair[1], pair[2], sep = "_")
            })
            
            # Organize the results for the current method
            jaccard_comparison <- list(ji_true_pairs=ji_true_pairs, ji_omic_one_pairs=ji_omic_one_pairs, 
                                       ji_omic_two_pairs=ji_omic_two_pairs, ji_smp_pairs=ji_smp_pairs,
                                       ji_matrix_pairs=ji_matrix_pairs)
            
            # Create a unique name for each dataset
            dataset_name <- sprintf("variance_%d_iteration_%d_varphi_%d", sigma, i, varphi) 
            jaccard_comparison_results[[dataset_name]] <- jaccard_comparison
            
            # COMPARISON METRICS: SENSITIVITY, SPECIFICITY; ACCURACY
            per_measures_true <- compare_features(Data_Ft2)
            per_measures_omic_one <- compare_features(Data_Ft.a.one)
            per_measures_omic_two <- compare_features(Data_Ft.b.two)
            per_measures_smp <- compare_sample(Data_Smp2)
            
            # Organize the results for the current method
            per_measures <- list(per_measures_true=per_measures_true, per_measures_omic_one=per_measures_omic_one, 
                                 per_measures_omic_two=per_measures_omic_two, per_measures_smp=per_measures_smp)
            
            # Create a unique name for each dataset
            dataset_name <- sprintf("variance_%d_iteration_%d_varphi_%d", sigma, i, varphi) 
            per_measures_results[[dataset_name]] <- per_measures
            
          }
          #simulatedData <- simulated_data[[paste0("iteration_", i)]]$concatenated_datasets[[1]]
          # Create a list of datasets
          datasets <- list(
            simulatedData = current_data,
            Data_Ft = Data_Ft,
            Data_Smp = Data_Smp
          ) 
          dataset <- sprintf("variance_%d_iteration_%d", sigma, i)
          dataset_output[[dataset]] <- datasets
          
        }, silent = TRUE)
      }
      
      ## ------------------------- Random Shared Factor -------------------------------
    } else if (method == "random.shared.factor") {
      #current_seed <- sample(1e6, 1)
      #set.seed(current_seed)
      
      # Simulate data
      simulated_data <- random.shared.factor(
        n_features_one = n_features_one,
        n_features_two = n_features_two,
        n_samples = n_samples,
        sigmas = sigma,
        iterations = iter
      )
      
      for (i in 1:length(simulated_data)) {
        try({
          current_data <- simulated_data[[paste0("iteration_", i)]]$concatenated_datasets[[1]]
          indices_features.1 <- simulated_data[[paste0("iteration_", i)]]$indices_features.1
          indices_features.2 <- simulated_data[[paste0("iteration_", i)]]$indices_features.2
          indices_samples.1 <- simulated_data[[paste0("iteration_", i)]]$indices_samples.1
          indices_samples.2 <- simulated_data[[paste0("iteration_", i)]]$indices_samples.2
          
          print(paste("Processing iteration:", i))
          # Run factorization methods
          fabia_result <- func_fabia(current_data, num_biclusters)
          mofa_result <- func_mofa(current_data, num_biclusters)
          mfa_result <- func_mfa(current_data, num_biclusters)
          gfa_result <- func_gfa(current_data, num_biclusters)
          
          if (!is.null(current_data) && nrow(current_data) > 0) {
            mean_cvs_features <- mean(apply(t(current_data), 2, function(x) sd(x) / abs(mean(x))))
            mean_cvs_samples <- mean(apply(current_data, 2, function(x) sd(x) / abs(mean(x))))
            
          } else {
            stop("Simulated datasets are empty or invalid.")
          }
          
          random.data <- data.frame(t(current_data)); random.data$feature = rownames(random.data) 
          random.data$ID <- ifelse(grepl("omic1_", random.data$feature),
                                   sub(".*omic1_feature_(\\d+)", "\\1", random.data$feature),
                                   ifelse(grepl("omic2_", random.data$feature),
                                          sub(".*omic2_feature_(\\d+)", "\\1", random.data$feature),NA))
          
          # Create dataview i.e that will help separate the two merged dataset
          random.data <- random.data %>%
            mutate(dataview = ifelse(grepl("omic1_", feature), "omic.one", ifelse(grepl("omic2_", feature), "omic.two", NA)))
          
          # Filtering the dataset
          random.data.a <- random.data %>% filter(dataview == "omic.one")
          random.data.b <- random.data %>% filter(dataview == "omic.two")
          
          # Assigning 'indices_features.1' to 'random.data.a'
          indices_features.1 = simulated_data[[paste0("iteration_", i)]]$indices_features.1
          in_range.a <- random.data.a$ID %in% indices_features.1
          simulated_features_a <- data.frame(feature = random.data.a$feature, signal = in_range.a)
          
          # Assigning 'indices_features.1' to 'random.data.a'
          indices_features.2 = simulated_data[[paste0("iteration_", i)]]$indices_features.2
          in_range.b <- random.data.b$ID %in% indices_features.2
          simulated_features_b <- data.frame(feature = random.data.b$feature, signal = in_range.b)
          
          # Merge 'simulated_features_a' and 'simulated_features_b'
          simulated_features_sigma = rbind(simulated_features_a, simulated_features_b)
          
          # samples
          indices_samples1 <- paste("sample_", simulated_data[[paste0("iteration_", i)]]$indices_samples.1, sep = "")
          indices_samples2 <- paste("sample_", simulated_data[[paste0("iteration_", i)]]$indices_samples.2, sep = "")
          indices_samples <- unique(c(indices_samples1, indices_samples2))
          in_range_sample <- rownames(current_data) %in% rownames(current_data[indices_samples,])
          simulated_samples_sigma <- data.frame(sample = rownames(current_data),
                                                signal = in_range_sample)
          
          # Extract loadings and scores
          loading_data.not.scaled <- data.frame(
            feature = mofa_result$mofa_loading$feature,
            loading_FABIA = fabia_result$fab_loading$loading_FABIA,
            loading_MOFA = mofa_result$mofa_loading$loading_MOFA,
            loading_MFA = mfa_result$mfa_loading$loading_MFA,
            loading_GFA = gfa_result$gfa_loading$loading_GFA.V
          )
          
          score_data.not.scaled <- data.frame(
            sample = mofa_result$mofa_score$sample,
            score_FABIA = fabia_result$fab_score$score_FABIA,
            score_MOFA = mofa_result$mofa_score$score_MOFA,
            score_MFA = mfa_result$mfa_score$score_MFA,
            score_GFA = gfa_result$gfa_score$score_GFA.V
          )
          
          # Check if mean of FABIA loadings and scores is zero
          mean_fabia_loading <- mean(fabia_result$fab_loading$loading_FABIA, na.rm = TRUE)
          mean_fabia_score <- mean(fabia_result$fab_score$score_FABIA, na.rm = TRUE)
          
          # Convert data to absolute values
          #Convert all numeric columns to absolute values
          score_data.abs <- score_data.not.scaled %>% mutate(across(where(is.numeric), ~abs(.)))
          loading_data.abs <- loading_data.not.scaled %>% mutate(across(where(is.numeric), ~abs(.)))
          
          loading_data <- loading_data.abs #normalize_columns(loading_data.abs, columns_to_loadings_normalize)
          score_data <- score_data.abs #normalize_columns(score_data.abs, columns_to_scores_normalize)
          
          Data_Ft <- merge(x=loading_data,y=simulated_features_sigma, by="feature", all = TRUE)
          Data_Ft <- data.frame(Data_Ft)[, c('feature','loading_FABIA','loading_MOFA','loading_MFA','signal','loading_GFA')]
          
          # Mutate the dataview column based on the pattern in the feature column
          Data_Ft <- Data_Ft %>%
            mutate(dataview = ifelse(grepl("omic1_", feature), "omic.one", ifelse(grepl("omic2_", feature), "omic.two", NA)))
          
          Data_Ft$ID <- ifelse(grepl("omic1_", Data_Ft$feature),
                               sub(".*omic1_feature_(\\d+)", "\\1", Data_Ft$feature),
                               ifelse(grepl("omic2_", Data_Ft$feature),
                                      sub(".*omic2_feature_(\\d+)", "\\1", Data_Ft$feature),NA))
          
          Data_Smp <- merge(x=score_data,y=simulated_samples_sigma, by="sample", all = TRUE)
          Data_Smp <- data.frame(Data_Smp)[, c('sample', 'score_FABIA','score_MOFA','score_MFA','signal','score_GFA')]
          
          # Subset omic1
          Data_Ft.a <- Data_Ft %>%
            filter(dataview == "omic.one")
          
          # Subset omic2
          Data_Ft.b <- Data_Ft %>%
            filter(dataview == "omic.two")
          
          for (varphi in seq_along(varphi_functions)) {
            
            varphi_function <- names(varphi_functions)[[varphi]]
            
            # Apply the function to Data_Ft.a
            # Apply the function to Data_Ft.a
            if (varphi_function == "varphi.five"){
              # Apply the function to Data_Ft.a
              Data_Ft.a.one <- Data_Ft.a %>%
                mutate(
                  fabia_index = ifelse(
                    abs(loading_FABIA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_FABIA,
                      omic = 1
                    )), 1, 0),
                  
                  mofa_index = ifelse(
                    abs(loading_MOFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_MOFA,
                      omic = 1
                    )), 1, 0),
                  
                  mfa_index = ifelse(
                    abs(loading_MFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_MFA,
                      omic = 1
                    )), 1, 0),
                  
                  gfa_index = ifelse(
                    abs(loading_GFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_GFA,
                      omic = 1
                    )), 1, 0),
                  
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
              
            }else{
              # Apply the function to Data_Ft.a
              Data_Ft.a.one <- Data_Ft.a %>%
                mutate(
                  fabia_index = ifelse(abs(loading_FABIA) >= get(varphi_function)(loading_FABIA), 1, 0),
                  mofa_index = ifelse(abs(loading_MOFA) >= get(varphi_function)(loading_MOFA), 1, 0),
                  mfa_index = ifelse(abs(loading_MFA) >= get(varphi_function)(loading_MFA), 1, 0),
                  gfa_index = ifelse(abs(loading_GFA) >= get(varphi_function)(loading_GFA), 1, 0),
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
            }
            
            if (varphi_function == "varphi.five"){
              Data_Ft.b.two <- Data_Ft.b %>%
                mutate(
                  fabia_index = ifelse(
                    abs(loading_FABIA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_FABIA,
                      omic = 2
                    )), 1, 0),
                  
                  mofa_index = ifelse(
                    abs(loading_MOFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_MOFA,
                      omic = 2
                    )), 1, 0),
                  
                  mfa_index = ifelse(
                    abs(loading_MFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_MFA,
                      omic = 2
                    )), 1, 0),
                  
                  gfa_index = ifelse(
                    abs(loading_GFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_GFA,
                      omic = 2
                    )), 1, 0),
                  
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
            } else{
              Data_Ft.b.two <- Data_Ft.b %>%
                mutate(
                  fabia_index = ifelse(abs(loading_FABIA) >= get(varphi_function)(loading_FABIA), 1, 0),
                  mofa_index = ifelse(abs(loading_MOFA) >= get(varphi_function)(loading_MOFA), 1, 0),
                  mfa_index = ifelse(abs(loading_MFA) >= get(varphi_function)(loading_MFA), 1, 0),
                  gfa_index = ifelse(abs(loading_GFA) >= get(varphi_function)(loading_GFA), 1, 0),
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
            }
            
            if (varphi_function == "varphi.five"){
              Data_Ft2 <- Data_Ft %>%
                mutate(
                  fabia_index = ifelse(
                    abs(loading_FABIA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_FABIA,
                      omic = "ALL"
                    )), 1, 0),
                  
                  mofa_index = ifelse(
                    abs(loading_MOFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_MOFA,
                      omic = "ALL"
                    )), 1, 0),
                  
                  mfa_index = ifelse(
                    abs(loading_MFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_MFA,
                      omic = "ALL"
                    )), 1, 0),
                  
                  gfa_index = ifelse(
                    abs(loading_GFA) >= do.call(get(varphi_function), list(
                      loading_or_score = loading_GFA,
                      omic = "ALL"
                    )), 1, 0),
                  
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
            } else{
              Data_Ft2 <- Data_Ft %>%
                mutate(
                  fabia_index = ifelse(abs(loading_FABIA) >= get(varphi_function)(loading_FABIA), 1, 0),
                  mofa_index = ifelse(abs(loading_MOFA) >= get(varphi_function)(loading_MOFA), 1, 0),
                  mfa_index = ifelse(abs(loading_MFA) >= get(varphi_function)(loading_MFA), 1, 0),
                  gfa_index = ifelse(abs(loading_GFA) >= get(varphi_function)(loading_GFA), 1, 0),
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
            }
            
            if (varphi_function == "varphi.five"){
              Data_Smp2 <- Data_Smp %>%
                mutate(
                  fabia_index = ifelse(
                    abs(score_FABIA) >= do.call(get(varphi_function), list(
                      loading_or_score = score_FABIA
                    )), 1, 0),
                  
                  mofa_index = ifelse(
                    abs(score_MOFA) >= do.call(get(varphi_function), list(
                      loading_or_score = score_MOFA
                    )), 1, 0),
                  
                  mfa_index = ifelse(
                    abs(score_MFA) >= do.call(get(varphi_function), list(
                      loading_or_score = score_MFA
                    )), 1, 0),
                  
                  gfa_index = ifelse(
                    abs(score_GFA) >= do.call(get(varphi_function), list(
                      loading_or_score = score_GFA
                    )), 1, 0),
                  
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
              
            } else{
              Data_Smp2 <- Data_Smp %>%
                mutate(
                  fabia_index = ifelse(abs(score_FABIA) >= get(varphi_function)(score_FABIA), 1, 0),
                  mofa_index = ifelse(abs(score_MOFA) >= get(varphi_function)(score_MOFA), 1, 0),
                  mfa_index = ifelse(abs(score_MFA) >= get(varphi_function)(score_MFA), 1, 0),
                  gfa_index = ifelse(abs(score_GFA) >= get(varphi_function)(score_GFA), 1, 0),
                  signal_index = ifelse(signal == 'TRUE', 1, 0)
                )
            }
            
            
            # Create binary indices
            Data_Matrix_Ft <- Data_Ft2
            
            # Create binary indices
            Data_Matrix_Smp <- Data_Smp2
            
            # Signal MATRIX
            sample_signal_index = Data_Matrix_Smp$signal_index
            feature_signal_index = Data_Matrix_Ft$signal_index
            signal_matrix = as.matrix(feature_signal_index)%*%as.matrix(t(sample_signal_index))
            
            # loading_score_matrix
            fabia_matrix = as.matrix(Data_Matrix_Ft$fabia_index)%*%as.matrix(t(Data_Matrix_Smp$fabia_index))
            mofa_matrix = as.matrix(Data_Matrix_Ft$mofa_index)%*%as.matrix(t(Data_Matrix_Smp$mofa_index))
            mfa_matrix = as.matrix(Data_Matrix_Ft$mfa_index)%*%as.matrix(t(Data_Matrix_Smp$mfa_index))
            gfa_matrix = as.matrix(Data_Matrix_Ft$gfa_index)%*%as.matrix(t(Data_Matrix_Smp$gfa_index))
            
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ JACCARD INDEX ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            # List of methods and corresponding indices
            methods <- c("fabia", "mofa", "mfa", "gfa")
            
            # Compute Jaccard indices in a single step
            ji_true <- sapply(methods, function(method) {
              method_index <- Data_Ft2[[paste0(method, "_index")]]
              jaccard.index.sim(method_index, Data_Ft2$signal_index)
            }, simplify = TRUE, USE.NAMES = TRUE)
            
            ji_omic_one <- sapply(methods, function(method) {
              method_index <- Data_Ft.a.one[[paste0(method, "_index")]]
              jaccard.index.sim(method_index, Data_Ft.a.one$signal_index)
            }, simplify = TRUE, USE.NAMES = TRUE)
            
            ji_omic_two <- sapply(methods, function(method) {
              method_index <- Data_Ft.b.two[[paste0(method, "_index")]]
              jaccard.index.sim(method_index, Data_Ft.b.two$signal_index)
            }, simplify = TRUE, USE.NAMES = TRUE)
            
            ji_samples <- sapply(methods, function(method) {
              method_index <- Data_Smp2[[paste0(method, "_index")]]
              jaccard.index.sim(method_index, Data_Smp2$signal_index)
            }, simplify = TRUE, USE.NAMES = TRUE)
            
            # Define matrices
            matrices <- list(
              fabia = fabia_matrix,
              mofa = mofa_matrix,
              mfa = mfa_matrix,
              gfa = gfa_matrix
            )
            
            # Compute Jaccard similarity indices
            ji_matrix <- sapply(names(matrices), function(method) {
              jaccard.index.sim(matrices[[method]], signal_matrix)
            }, simplify = TRUE, USE.NAMES = TRUE)
            
            
            # Organize the results for the current method
            jaccard_loading <- list(ji_true=ji_true, ji_omic_one=ji_omic_one, ji_omic_two=ji_omic_two, ji_samples=ji_samples, ji_matrix=ji_matrix)
            
            # Create a unique name for each dataset
            dataset_name <- sprintf("variance_%d_iteration_%d_varphi_%d", sigma, i, varphi) 
            jaccard_results[[dataset_name]] <- jaccard_loading
            
            # COMPARISONS BY METHODS
            methods <- c("fabia", "mofa", "mfa", "gfa")  # Define the list of methods
            
            # ALL DATA
            # Calculate pairwise Jaccard indices dynamically
            ji_true_pairs <- combn(methods, 2, function(method_pair) {
              method1_index <- Data_Ft2[[paste0(method_pair[1], "_index")]]
              method2_index <- Data_Ft2[[paste0(method_pair[2], "_index")]]
              jaccard.index.sim(method1_index, method2_index)
            }, simplify = TRUE)
            
            # Assign names to the results
            names(ji_true_pairs) <- combn(methods, 2, function(method_pair) {
              paste(method_pair[1], method_pair[2], "ji_phi", sep = "_")
            })
            
            # Convert to a named vector or data frame for further use
            ji_true_pairs <- as.data.frame(as.list(ji_true_pairs))
            
            # OMIC ONE
            # Calculate pairwise Jaccard indices dynamically
            ji_omic_one_pairs <- combn(methods, 2, function(method_pair) {
              method1_index <- Data_Ft.a.one[[paste0(method_pair[1], "_index")]]
              method2_index <- Data_Ft.a.one[[paste0(method_pair[2], "_index")]]
              jaccard.index.sim(method1_index, method2_index)
            }, simplify = TRUE)
            
            # Assign names to the results
            names(ji_omic_one_pairs) <- combn(methods, 2, function(method_pair) {
              paste(method_pair[1], method_pair[2], "ji_phi", sep = "_")
            })
            
            # Convert to a named vector or data frame for further use
            ji_omic_one_pairs <- as.data.frame(as.list(ji_omic_one_pairs))
            
            # OMIC TWO
            # Calculate pairwise Jaccard indices dynamically
            ji_omic_two_pairs <- combn(methods, 2, function(method_pair) {
              method1_index <- Data_Ft.b.two[[paste0(method_pair[1], "_index")]]
              method2_index <- Data_Ft.b.two[[paste0(method_pair[2], "_index")]]
              jaccard.index.sim(method1_index, method2_index)
            }, simplify = TRUE)
            
            # Assign names to the results
            names(ji_omic_two_pairs) <- combn(methods, 2, function(method_pair) {
              paste(method_pair[1], method_pair[2], "ji_phi", sep = "_")
            })
            
            # Convert to a named vector or data frame for further use
            ji_omic_two_pairs <- as.data.frame(as.list(ji_omic_two_pairs))
            
            # SAMPLES
            # Calculate pairwise Jaccard indices dynamically
            ji_smp_pairs <- combn(methods, 2, function(method_pair) {
              method1_index <- Data_Smp2[[paste0(method_pair[1], "_index")]]
              method2_index <- Data_Smp2[[paste0(method_pair[2], "_index")]]
              jaccard.index.sim(method1_index, method2_index)
            }, simplify = TRUE)
            
            # Assign names to the results
            names(ji_smp_pairs) <- combn(methods, 2, function(method_pair) {
              paste(method_pair[1], method_pair[2], "ji_phi", sep = "_")
            })
            
            # Convert to a named vector or data frame for further use
            ji_smp_pairs <- as.data.frame(as.list(ji_smp_pairs))
            
            # MATRICES
            
            # Define the matrices in a named list
            matrices <- list(
              fabia = fabia_matrix,
              mofa = mofa_matrix,
              mfa = mfa_matrix,
              gfa = gfa_matrix
            )
            
            # Get all unique pair combinations of matrix names
            matrix_pairs <- combn(names(matrices), 2, simplify = FALSE)
            
            # Compute Jaccard similarity indices for each pair
            ji_matrix_pairs <- sapply(matrix_pairs, function(pair) {
              jaccard.index.sim(matrices[[pair[1]]], matrices[[pair[2]]])
            }, simplify = TRUE, USE.NAMES = FALSE)
            
            # Assign names to the results
            names(ji_matrix_pairs) <- sapply(matrix_pairs, function(pair) {
              paste(pair[1], pair[2], sep = "_")
            })
            
            # Organize the results for the current method
            jaccard_comparison <- list(ji_true_pairs=ji_true_pairs, ji_omic_one_pairs=ji_omic_one_pairs, 
                                       ji_omic_two_pairs=ji_omic_two_pairs, ji_smp_pairs=ji_smp_pairs,
                                       ji_matrix_pairs=ji_matrix_pairs)
            
            # Create a unique name for each dataset
            dataset_name <- sprintf("variance_%d_iteration_%d_varphi_%d", sigma, i, varphi) 
            jaccard_comparison_results[[dataset_name]] <- jaccard_comparison
            
            # COMPARISON METRICS: SENSITIVITY, SPECIFICITY; ACCURACY
            per_measures_true <- compare_features(Data_Ft2)
            per_measures_omic_one <- compare_features(Data_Ft.a.one)
            per_measures_omic_two <- compare_features(Data_Ft.b.two)
            per_measures_smp <- compare_sample(Data_Smp2)
            
            # Organize the results for the current method
            per_measures <- list(per_measures_true=per_measures_true, per_measures_omic_one=per_measures_omic_one, per_measures_omic_two=per_measures_omic_two, per_measures_smp=per_measures_smp)
            
            # Create a unique name for each dataset
            dataset_name <- sprintf("variance_%d_iteration_%d_varphi_%d", sigma, i, varphi) 
            per_measures_results[[dataset_name]] <- per_measures
            simulatedData <- simulated_data[[paste0("iteration_", 1)]]$concatenated_datasets[[1]]
            # Create a list of datasets
            datasets <- list(
              simulatedData = simulatedData,
              Data_Ft = Data_Ft,
              Data_Smp = Data_Smp
            )          }
          
          datasets <- list(
            simulatedData = current_data,
            Data_Ft = Data_Ft,
            Data_Smp = Data_Smp
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
    per_measures_results = per_measures_results#,
    #dataset_output = dataset_output
    )
  
  return(results)
  
}

# -------------------------- ACTUAL SIMULATION: Define the actual simulation parameters -------------------------- #

# Define the methods and sigmas
methods <- c("single.factor")#'single.factor', 'shared.factor', 'random.shared.factor'
sigmas <- 15  # Define sigma or the range of sigmas

# Loop over sigmas
for (sigma in sigmas) {
  # Loop over methods for the current sigma
  for (method in methods) {
    # Generate a unique result name
    result_name <- paste0("sim_omicResults_", method, "_sigma_", sigma)

    # Run the simulation for the current method and sigma
    sim_result <- sim_OmicR(
      n_features_one = 4000,
      n_features_two = 3000,
      n_samples = 100,
      var_sigma = sigma,
      num_biclusters = 1,
      num_iterations = 100,
      method = method
    )

    # Save the result to a file
    saveRDS(
      sim_result,
      file = paste0("/user/leuven/364/vsc36498/", result_name, ".rds")
      #file = paste0("C:/Users/Lenovo/Downloads/", result_name, ".rds")
      #file = paste0("C:/Users/bosangir/Downloads/", result_name, ".rds")
    )
  }

  # Print progress message
  cat("Completed simulations for sigma =", sigma, "\n")
}

# Print completion message
cat("All simulations completed and saved!\n")
# library(doParallel)
# library(foreach)
# 
# # Define your simulation parameters
# methods <- c("single.factor")  # Add more methods if needed
# sigmas <- c(2:3)                 # Use a vector if testing multiple values
# 
# # Register parallel backend
# num_cores <- parallel::detectCores() - 1  # Leave one core free
# cl <- makeCluster(num_cores)
# registerDoParallel(cl)
# 
# # Export sim_OmicR to each worker
# clusterExport(cl, varlist = c("sim_OmicR"))
# 
# # Run simulations in parallel
# results <- foreach(sigma = sigmas, .combine = 'list') %:%
#   foreach(method = methods, .combine = 'list') %dopar% {
#     
#     result_name <- paste0("sim_omicResults_", method, "_sigma_", sigma)
#     
#     sim_result <- sim_OmicR(
#       n_features_one = 4000,
#       n_features_two = 3000,
#       n_samples = 100,
#       var_sigma = sigma,
#       num_biclusters = 1,
#       num_iterations = 50,
#       method = method
#     )
#     
#     saveRDS(
#       sim_result,
#       file = paste0("/user/leuven/364/vsc36498/", result_name, ".rds")
#     )
#     
#     return(result_name)  # Return name for reference/logging
#   }
# 
# # Stop the cluster
# stopCluster(cl)
# 
# # Print completion summary
# cat("All simulations completed and saved!\n")
# #print(unlist(results))
