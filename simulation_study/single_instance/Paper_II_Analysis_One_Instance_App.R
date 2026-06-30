#setwd("/user/leuven/364/vsc36498")
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
        
        alpha <- rnorm(n_s, 0, 0.5)
        
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
          list_alphas[[paste0("alpha", i)]] <- rnorm(n_samples, 0, 0.5)
        }
        
        # Assign corresponding values to alpha variables based on assigned_indices_samples
        for (i in seq_along(assigned_indices_samples)) {
          indices <- assigned_indices_samples[[i]]
          list_alphas[[i]][indices] <- rnorm(length(indices), (3 + 0.5*i), 1.0)  # Adjust values dynamically
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
          list_betas[[i]] <- rnorm(n_features_one, 0, 0.5)
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
            list_betas[[i]][indices_ns] <- rnorm(length(indices_ns), mean = (4.0 + 0.5 * i), sd = 1.0)
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
          list_deltas[[i]] <- rnorm(n_features_two, 0, 0.5)
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
            list_deltas[[i]][indices_ns] <- rnorm(length(indices_ns), mean = (5.0 + 0.0 * i), sd = 1.0)
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
      factor_scores = list_alphas,
      loading_omic1 = list_betas,
      loading_omic2 = list_deltas
    )
    }

return(all_omic_data)
}

simulatedData <- multiple_factor(n_features_one = 4000, n_features_two = 3000, n_samples = 100, sigmas = 7, n_factors = 2, iterations = 1)
   
dataset_overall <- simulatedData$iteration_1$concatenated_datasets[[1]]
dataset_omic1 <- dataset_overall[,c(1:4000)]
dataset_omic2 <- dataset_overall[,c(4001:7000)]
dataset = dataset_omic1
image(c(1:dim(dataset)[1]), c(1:dim(dataset)[2]), dataset, ylab = "Features", xlab = "Samples", font.lab = 2)
dataset = dataset_omic2
image(c(1:dim(dataset)[1]), c(1:dim(dataset)[2]), dataset, ylab = "Features", xlab = "Samples", font.lab = 2)

# ----------------------------------------- Visualize the raw/true signals ------------------------------
# Factor scores
plot_factor_with_indices <- function(
    scores, signal_idx, factor_id = 1,
    shuffle = c("signal","all","none"), seed = 123,
    palette = c(Noise = "grey35", Signal = "#FFD400"),
    x_breaks = NULL, x_limits = NULL
){
  shuffle <- match.arg(shuffle)
  s <- if (is.matrix(scores) || is.data.frame(scores)) as.numeric(scores) else as.numeric(scores)
  n <- length(s)
  
  if (is.numeric(signal_idx)) {
    sig_num <- as.integer(signal_idx)
  } else {
    sig_num <- suppressWarnings(as.integer(gsub("[^0-9]+","", signal_idx)))
    sig_num <- sig_num[!is.na(sig_num)]
  }
  sig_num <- intersect(sig_num, seq_len(n))
  
  df <- data.frame(
    sample_id = seq_len(n),
    sample_name = paste0("sample_", seq_len(n)),
    score = s,
    signal = ifelse(seq_len(n) %in% sig_num, "Signal", "Noise")
  )
  
  set.seed(seed)
  df$plot_x <- df$sample_id
  if (shuffle == "all") {
    df$plot_x <- sample(df$plot_x, size = n)
  } else if (shuffle == "signal") {
    i <- df$signal == "Signal"
    df$plot_x[i] <- sample(df$plot_x[i], size = sum(i))
  }
  df$xnum <- as.numeric(df$plot_x)  # numeric axis
  
  library(ggplot2)
  
  p_scatter <- ggplot(df, aes(x = xnum, y = score, colour = signal)) +
    geom_point(size = 2, alpha = 0.9) +
    scale_colour_manual(values = palette, drop = FALSE) +
    labs(title = paste0("Factor ", factor_id), x = "Samples", y = "Factor score") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  
  if (!is.null(x_breaks) || !is.null(x_limits)) {
    p_scatter <- p_scatter + scale_x_continuous(breaks = x_breaks, limits = x_limits)
  }
  
  p_violin <- ggplot(df, aes(x = signal, y = score, fill = signal)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = 19, outlier.alpha = 0.6) +
    scale_fill_manual(values = palette, drop = FALSE) +
    labs(title = "", x = NULL, y = "Factor score") + # Distribution by signal
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  
  if (requireNamespace("patchwork", quietly = TRUE)) {
    print(p_scatter + p_violin + patchwork::plot_layout(widths = c(3, 2)))
  } else {
    print(p_scatter); print(p_violin)
  }
  invisible(list(data = df, scatter = p_scatter, violin = p_violin))
}

# Factor 1
plot_factor_with_indices(
  scores     = simulatedData$iteration_1$factor_scores[[1]],
  signal_idx = simulatedData$iteration_1$indices_samples[[1]],
  factor_id  = 1,
  shuffle    = "signal",
  x_breaks   = seq(0, 100, 20),
  x_limits   = c(0, 100)
)

# Factor 2
plot_factor_with_indices(
  scores     = simulatedData$iteration_1$factor_scores[[2]],
  signal_idx = simulatedData$iteration_1$indices_samples[[2]],
  factor_id  = 2,
  shuffle    = "signal",
  x_breaks   = seq(0, 100, 20),
  x_limits   = c(0, 100)
)

# Loadings
# install.packages("ggplot2") # if needed

plot_overlay_violin_scatter <- function(
    values,                 # loadings or scores (numeric vector or 1-col matrix)
    signal_idx,             # indices (numeric) or names like "feature_18"
    factor_id = 1,
    kind = c("Features","Samples"),   # x-axis label
    shuffle = c("none","signal","all"),
    seed = 123,
    palette = c(Noise = "grey55", Signal = "#FFD400"),
    x_limits = NULL,                 # e.g., c(0,100)
    x_breaks = NULL,                 # e.g., seq(0,100,20)
    violin_alpha = 0.25,             # transparency of overlaid violins
    violin_width_factor = 0.9,       # fraction of each group's x-range to span
    box_width_factor = 0.18          # fraction of violin width for the boxplot
){
  shuffle <- match.arg(shuffle)
  kind <- match.arg(kind)
  
  v <- if (is.matrix(values) || is.data.frame(values)) as.numeric(values) else as.numeric(values)
  n <- length(v)
  
  # accept numeric indices or names with numbers
  if (is.numeric(signal_idx)) {
    sig_num <- as.integer(signal_idx)
  } else {
    sig_num <- suppressWarnings(as.integer(gsub("[^0-9]+","", signal_idx)))
    sig_num <- sig_num[!is.na(sig_num)]
  }
  sig_num <- intersect(sig_num, seq_len(n))
  
  df <- data.frame(
    id     = seq_len(n),
    value  = v,
    signal = ifelse(seq_len(n) %in% sig_num, "Signal", "Noise"),
    stringsAsFactors = FALSE
  )
  
  set.seed(seed)
  df$plot_x <- df$id
  if (shuffle == "all") {
    df$plot_x <- sample(df$plot_x, n)
  } else if (shuffle == "signal") {
    i <- df$signal == "Signal"
    df$plot_x[i] <- sample(df$plot_x[i], sum(i))
  }
  
  # split groups
  d_noise  <- df[df$signal == "Noise", ]
  d_signal <- df[df$signal == "Signal", ]
  
  # centers and widths in *x-axis units* so violins sit over the clouds
  cx_noise  <- stats::median(d_noise$plot_x)
  cx_signal <- stats::median(d_signal$plot_x)
  
  rng_noise  <- diff(range(d_noise$plot_x))
  rng_signal <- diff(range(d_signal$plot_x))
  
  w_noise  <- max(1, rng_noise  * violin_width_factor)
  w_signal <- max(1, rng_signal * violin_width_factor)
  
  d_noise$x_violin  <- cx_noise
  d_signal$x_violin <- cx_signal
  
  library(ggplot2)
  
  p <- ggplot(df, aes(x = plot_x, y = value, colour = signal)) +
    geom_point(size = 1.9, alpha = 0.9) +
    scale_colour_manual(values = palette, drop = FALSE) +
    scale_fill_manual(values   = palette, drop = FALSE) +
    labs(title = paste0("Factor ", factor_id),
         x = kind, y = if (kind == "Features") "Loading" else "Factor score") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom")
  
  if (!is.null(x_limits) || !is.null(x_breaks)) {
    p <- p + scale_x_continuous(limits = x_limits, breaks = x_breaks)
  }
  
  # overlay violins + boxplots exactly where each group lives on x
  p <- p +
    geom_violin(
      data = d_noise,
      aes(x = x_violin, y = value, fill = signal, group = 1),
      width = w_noise, alpha = violin_alpha, colour = NA, inherit.aes = FALSE
    ) +
    geom_boxplot(
      data = d_noise,
      aes(x = x_violin, y = value, fill = signal, group = 1),
      width = w_noise * box_width_factor, outlier.shape = 19, outlier.alpha = 0.5,
      inherit.aes = FALSE
    ) +
    geom_violin(
      data = d_signal,
      aes(x = x_violin, y = value, fill = signal, group = 1),
      width = w_signal, alpha = violin_alpha, colour = NA, inherit.aes = FALSE
    ) +
    geom_boxplot(
      data = d_signal,
      aes(x = x_violin, y = value, fill = signal, group = 1),
      width = w_signal * box_width_factor, outlier.shape = 19, outlier.alpha = 0.5,
      inherit.aes = FALSE
    )
  
  print(p)
  invisible(p)
}

# Dataset 1
plot_overlay_violin_scatter(
  values     = simulatedData$iteration_1$loading_omic1[[1]],
  signal_idx = simulatedData$iteration_1$indices_features.1[[1]],
  factor_id  = 1,
  kind       = "Features",
  shuffle    = "none",
  x_limits   = c(0, 4000),
  x_breaks   = seq(0, 4000, 1000),
  violin_alpha = 0.30
)

plot_overlay_violin_scatter(
  values     = simulatedData$iteration_1$loading_omic1[[2]],
  signal_idx = simulatedData$iteration_1$indices_features.1[[2]],
  factor_id  = 1,
  kind       = "Features",
  shuffle    = "none",
  x_limits   = c(0, 4000),
  x_breaks   = seq(0, 4000, 1000),
  violin_alpha = 0.30
)

# Dataset 2
plot_overlay_violin_scatter(
  values     = simulatedData$iteration_1$loading_omic2[[1]],
  signal_idx = simulatedData$iteration_1$indices_features.2[[1]],
  factor_id  = 1,
  kind       = "Features",
  shuffle    = "none",
  x_limits   = c(0, 3000),
  x_breaks   = seq(0, 3000, 1000),
  violin_alpha = 0.30
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

first_omic <- dataset_omic1
second_omic <- dataset_omic2

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
factor_mofa2$sample_id <- as.numeric(gsub("sample_", "", factor_mofa2$sample))

# sort in ascending order by sample_id
factor_mofa2 <- factor_mofa2[order(factor_mofa2$sample_id), ]

# Factor 1
plot_factor_with_indices(
  scores     = factor_mofa2$score_MOFA1,
  signal_idx = simulatedData$iteration_1$indices_samples[[1]],
  factor_id  = 1,
  shuffle    = "signal",
  x_breaks   = seq(0, 100, 20),
  x_limits   = c(0, 100)
)

# Factor 2
plot_factor_with_indices(
  scores     = factor_mofa2$score_MOFA2,
  signal_idx = simulatedData$iteration_1$indices_samples[[2]],
  factor_id  = 2,
  shuffle    = "signal",
  x_breaks   = seq(0, 100, 20),
  x_limits   = c(0, 100)
)

#  Weights
# Weights/Loading 
loading_mofa<- get_weights(MOFAobject.sim.trained, factors = 'all', as.data.frame = T)
loading_mofa2 <- loading_mofa %>%
  spread(factor, value)
colnames(loading_mofa2) =  c('feature', 'view', 'loading_MOFA1', 'loading_MOFA2')

omic1 <- loading_mofa2 %>%
  filter(view == 'Omic 1') %>%
  select(feature, view, loading_MOFA1, loading_MOFA2)
omic1$feature_id1 <- as.numeric(gsub("omic1_feature_", "", omic1$feature))

omic1 <- omic1[order(omic1$feature_id1), ]

# Loadings
# Dataset 1
plot_overlay_violin_scatter(
  values     = omic1$loading_MOFA1,
  signal_idx = simulatedData$iteration_1$indices_features.1[[1]],
  factor_id  = 1,
  kind       = "Features",
  shuffle    = "none",
  x_limits   = c(0, 4000),
  x_breaks   = seq(0, 4000, 1000),
  violin_alpha = 0.30
)

plot_overlay_violin_scatter(
  values     = omic1$loading_MOFA2,
  signal_idx = simulatedData$iteration_1$indices_features.1[[2]],
  factor_id  = 1,
  kind       = "Features",
  shuffle    = "none",
  x_limits   = c(0, 4000),
  x_breaks   = seq(0, 4000, 1000),
  violin_alpha = 0.30
)

omic2 <- loading_mofa2 %>%
  filter(view == 'Omic 2') %>%
  select(feature, view, loading_MOFA1, loading_MOFA2)
omic2$feature_id2 <- as.numeric(gsub("omic2_feature_", "", omic2$feature))
omic2 <- omic2[order(omic2$feature_id2), ]


# Dataset 2
plot_overlay_violin_scatter(
  values     = omic2$loading_MOFA1,
  signal_idx = simulatedData$iteration_1$indices_features.2[[1]],
  factor_id  = 1,
  kind       = "Features",
  shuffle    = "none",
  x_limits   = c(0, 3000),
  x_breaks   = seq(0, 3000, 1000),
  violin_alpha = 0.30
)

plot_overlay_violin_scatter <- function(
    values,                 # loadings or scores (numeric vector or 1-col matrix)
    signal_idx = NULL,      # NULL/empty => no signal
    factor_id = 1,
    kind = c("Features","Samples"),
    shuffle = c("none","signal","all"),
    seed = 123,
    palette = c(Noise = "grey55", Signal = "#FFD400"),
    x_limits = NULL, x_breaks = NULL,
    violin_alpha = 0.25,
    violin_width_factor = 0.9,
    box_width_factor = 0.18
){
  shuffle <- match.arg(shuffle)
  kind <- match.arg(kind)
  
  # vectorize values
  v <- if (is.matrix(values) || is.data.frame(values)) as.numeric(values) else as.numeric(values)
  n <- length(v)
  
  # --- handle signal indices robustly ---
  if (is.null(signal_idx) || length(signal_idx) == 0) {
    sig_num <- integer(0)                # <-- no signal
  } else if (is.numeric(signal_idx)) {
    sig_num <- as.integer(signal_idx)
  } else {
    sig_num <- suppressWarnings(as.integer(gsub("[^0-9]+","", signal_idx)))
    sig_num <- sig_num[!is.na(sig_num)]
  }
  sig_num <- intersect(sig_num, seq_len(n))
  has_signal <- length(sig_num) > 0
  
  df <- data.frame(
    id     = seq_len(n),
    value  = v,
    signal = ifelse(seq_len(n) %in% sig_num, "Signal", "Noise"),
    stringsAsFactors = FALSE
  )
  
  set.seed(seed)
  df$plot_x <- df$id
  if (shuffle == "all") {
    df$plot_x <- sample(df$plot_x, n)
  } else if (shuffle == "signal" && has_signal) {
    i <- df$signal == "Signal"
    df$plot_x[i] <- sample(df$plot_x[i], sum(i))
  }
  
  d_noise  <- df[df$signal == "Noise", ]
  d_signal <- if (has_signal) df[df$signal == "Signal", ] else df[0, ]
  
  # centers/widths in x-units
  cx_noise <- stats::median(d_noise$plot_x)
  rng_noise <- diff(range(d_noise$plot_x))
  w_noise <- max(1, rng_noise * violin_width_factor)
  d_noise$x_violin <- cx_noise
  
  if (has_signal) {
    cx_signal <- stats::median(d_signal$plot_x)
    rng_signal <- diff(range(d_signal$plot_x))
    w_signal <- max(1, rng_signal * violin_width_factor)
    d_signal$x_violin <- cx_signal
  }
  
  library(ggplot2)
  
  p <- ggplot(df, aes(x = plot_x, y = value, colour = signal)) +
    geom_point(size = 1.9, alpha = 0.9) +
    # drop=TRUE hides the 'Signal' legend entry when absent
    scale_colour_manual(values = palette, drop = TRUE) +
    scale_fill_manual(values   = palette, drop = TRUE) +
    labs(title = paste0("Factor ", factor_id),
         x = kind, y = if (kind == "Features") "Loading" else "Factor score") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom")
  
  if (!is.null(x_limits) || !is.null(x_breaks)) {
    p <- p + scale_x_continuous(limits = x_limits, breaks = x_breaks)
  }
  
  # Noise violin + box
  p <- p +
    geom_violin(
      data = d_noise,
      aes(x = x_violin, y = value, fill = signal, group = 1),
      width = w_noise, alpha = violin_alpha, colour = NA, inherit.aes = FALSE
    ) +
    geom_boxplot(
      data = d_noise,
      aes(x = x_violin, y = value, fill = signal, group = 1),
      width = w_noise * box_width_factor, outlier.shape = 19, outlier.alpha = 0.5,
      inherit.aes = FALSE
    )
  
  # Signal layers only if present
  if (has_signal) {
    p <- p +
      geom_violin(
        data = d_signal,
        aes(x = x_violin, y = value, fill = signal, group = 1),
        width = w_signal, alpha = violin_alpha, colour = NA, inherit.aes = FALSE
      ) +
      geom_boxplot(
        data = d_signal,
        aes(x = x_violin, y = value, fill = signal, group = 1),
        width = w_signal * box_width_factor, outlier.shape = 19, outlier.alpha = 0.5,
        inherit.aes = FALSE
      )
  }
  
  print(p)
  invisible(p)
}

plot_overlay_violin_scatter <- function(
    values, signal_idx = NULL,
    factor_id = 1,
    kind = c("Features","Samples"),
    shuffle = c("none","signal","all"),
    seed = 123,
    palette = c(Noise = "grey55", Signal = "#FFD400"),
    x_limits = NULL, x_breaks = NULL,
    y_limits = NULL, y_breaks = NULL,   # <-- NEW
    violin_alpha = 0.25,
    violin_width_factor = 0.9,
    box_width_factor = 0.18
){
  shuffle <- match.arg(shuffle); kind <- match.arg(kind)
  v <- if (is.matrix(values) || is.data.frame(values)) as.numeric(values) else as.numeric(values)
  n <- length(v)
  
  # signal handling (NULL => all Noise)
  if (is.null(signal_idx) || length(signal_idx) == 0) {
    sig_num <- integer(0)
  } else if (is.numeric(signal_idx)) {
    sig_num <- as.integer(signal_idx)
  } else {
    sig_num <- suppressWarnings(as.integer(gsub("[^0-9]+","", signal_idx)))
    sig_num <- sig_num[!is.na(sig_num)]
  }
  sig_num <- intersect(sig_num, seq_len(n))
  has_signal <- length(sig_num) > 0
  
  df <- data.frame(
    id = seq_len(n), value = v,
    signal = ifelse(seq_len(n) %in% sig_num, "Signal", "Noise")
  )
  
  set.seed(seed)
  df$plot_x <- df$id
  if (shuffle == "all") df$plot_x <- sample(df$plot_x, n)
  if (shuffle == "signal" && has_signal) {
    i <- df$signal == "Signal"; df$plot_x[i] <- sample(df$plot_x[i], sum(i))
  }
  
  d_noise  <- df[df$signal == "Noise", ]
  d_signal <- if (has_signal) df[df$signal == "Signal", ] else df[0, ]
  
  cx_noise <- stats::median(d_noise$plot_x)
  rng_noise <- diff(range(d_noise$plot_x))
  w_noise <- max(1, rng_noise * violin_width_factor)
  d_noise$x_violin <- cx_noise
  
  if (has_signal) {
    cx_signal <- stats::median(d_signal$plot_x)
    rng_signal <- diff(range(d_signal$plot_x))
    w_signal <- max(1, rng_signal * violin_width_factor)
    d_signal$x_violin <- cx_signal
  }
  
  library(ggplot2)
  p <- ggplot(df, aes(x = plot_x, y = value, colour = signal)) +
    geom_point(size = 1.9, alpha = 0.9) +
    scale_colour_manual(values = palette, drop = TRUE) +
    scale_fill_manual(values   = palette, drop = TRUE) +
    labs(title = paste0("Factor ", factor_id),
         x = kind, y = if (kind == "Features") "Loading" else "Factor score") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom")
  
  if (!is.null(x_limits) || !is.null(x_breaks)) {
    p <- p + scale_x_continuous(limits = x_limits, breaks = x_breaks)
  }
  if (!is.null(y_breaks)) p <- p + scale_y_continuous(breaks = y_breaks)
  if (!is.null(y_limits)) p <- p + coord_cartesian(ylim = y_limits)
  
  # noise violin+box
  p <- p +
    geom_violin(
      data = d_noise, aes(x = x_violin, y = value, fill = signal, group = 1),
      width = w_noise, alpha = violin_alpha, colour = NA, inherit.aes = FALSE
    ) +
    geom_boxplot(
      data = d_noise, aes(x = x_violin, y = value, fill = signal, group = 1),
      width = w_noise * box_width_factor, outlier.shape = 19, outlier.alpha = 0.5,
      inherit.aes = FALSE
    )
  
  # signal layers only if present
  if (has_signal) {
    p <- p +
      geom_violin(
        data = d_signal, aes(x = x_violin, y = value, fill = signal, group = 1),
        width = w_signal, alpha = violin_alpha, colour = NA, inherit.aes = FALSE
      ) +
      geom_boxplot(
        data = d_signal, aes(x = x_violin, y = value, fill = signal, group = 1),
        width = w_signal * box_width_factor, outlier.shape = 19, outlier.alpha = 0.5,
        inherit.aes = FALSE
      )
  }
  
  print(p); invisible(p)
}

plot_overlay_violin_scatter(
  values     = omic2$loading_MOFA2,
  signal_idx = NULL,
  factor_id  = 1,
  kind       = "Features",
  shuffle    = "none",
  y_limits   = c(-0.1, 0.4),
  x_limits   = c(0, 3000),
  x_breaks   = seq(0, 3000, 1000),
  violin_alpha = 0.30
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

# sort in ascending order by sample_id
df_fabia_factor <- df_fabia_factor[order(df_fabia_factor$sample_id), ]

# Factor 1
plot_factor_with_indices(
  scores     = (-1)*(df_fabia_factor$score_FABIA2),
  signal_idx = simulatedData$iteration_1$indices_samples[[1]],
  factor_id  = 1,
  shuffle    = "signal",
  x_breaks   = seq(0, 100, 20),
  x_limits   = c(0, 100)
)

# Factor 2
plot_factor_with_indices(
  scores     = (df_fabia_factor$score_FABIA1),
  signal_idx = simulatedData$iteration_1$indices_samples[[2]],
  factor_id  = 2,
  shuffle    = "signal",
  x_breaks   = seq(0, 100, 20),
  x_limits   = c(0, 100)
)

# Weights/Loading
loading_simfabia = data.frame(FABIAobject.sim.trained@L)
loading_simfabia$feature = rownames(loading_simfabia)
colnames(loading_simfabia) <- c("loading_FABIA1", "loading_FABIA2","features")
features_info.sim = data.frame(rownames(fabia.sim.scaledData)) 
colnames(features_info.sim) <- c("features") 
features_info.sim$ID <- as.integer(row.names(features_info.sim))
features_info.sim$view <- NA  # Create a new column with NAs
features_info.sim$view[1:4000] <- "Omic 1"  # Assign 'x' to rows 1 to 50
features_info.sim$view[4001:7000] <- 'Omic 2'  # Assign 'y' to rows 51 to 100

loading_simfabia = merge(loading_simfabia, features_info.sim, by = "features")

omic1_simfabia <- loading_simfabia %>%
  filter(view == 'Omic 1') %>%
  select(features, view, loading_FABIA1,loading_FABIA2, ID)
omic1_simfabia <- omic1_simfabia[order(omic1_simfabia$ID), ]

omic2_simfabia <- loading_simfabia %>%
  filter(view == 'Omic 2') %>%
  select(features, view, loading_FABIA1,loading_FABIA2, ID)
omic2_simfabia <- omic2_simfabia[order(omic2_simfabia$ID), ]

# Dataset 1
plot_overlay_violin_scatter(
  values     = -1*(omic1_simfabia$loading_FABIA2),
  signal_idx = simulatedData$iteration_1$indices_features.1[[1]],
  factor_id  = 1,
  kind       = "Features",
  shuffle    = "none",
  x_limits   = c(0, 4000),
  x_breaks   = seq(0, 4000, 1000),
  violin_alpha = 0.30
)

plot_overlay_violin_scatter(
  values     = omic1_simfabia$loading_FABIA1,
  signal_idx = simulatedData$iteration_1$indices_features.1[[2]],
  factor_id  = 1,
  kind       = "Features",
  shuffle    = "none",
  x_limits   = c(0, 4000),
  x_breaks   = seq(0, 4000, 1000),
  violin_alpha = 0.30
)

# Dataset 2
plot_overlay_violin_scatter(
  values     = (-1)*(omic2_simfabia$loading_FABIA2),
  signal_idx = simulatedData$iteration_1$indices_features.2[[1]],
  factor_id  = 1,
  kind       = "Features",
  shuffle    = "none",
  x_limits   = c(0, 3000),
  x_breaks   = seq(0, 3000, 1000),
  violin_alpha = 0.30
)

plot_overlay_violin_scatter(
  values     = omic2_simfabia$loading_FABIA2,
  signal_idx = NULL,
  factor_id  = 1,
  kind       = "Features",
  shuffle    = "none",
  y_limits   = c(-0.03, 0.1),
  x_limits   = c(0, 3000),
  x_breaks   = seq(0, 3000, 1000),
  violin_alpha = 0.30
)


# merge the two datasets
sim_data_fabia <- rbind(t(first_omic), t(second_omic))

# FABIA MODEL
set.seed(123)
FABIAobject.sim.trained_sparse <- fabia(sim_data_fabia,
                                 p = 2,           # number of hidden factors = number of biclusters, default = 5
                                 alpha = 0,     # sparseness loadings (0 - 1.0); default = 0.1
                                 cyc = 1000,      # number of iterations; default = 500
                                 spl = 0.5,       # sparseness prior loadings (0 - 2.0); default = 0 (Laplace)
                                 spz = 0.5,       # sparseness factors (0.5 - 2.0); default = 0.5 (Laplace)
                                 random = 1.0,    # random initialization of loadings in [-random,random]; default = 1.0.
                                 center = 2,      # data centering: 1 (mean), 2 (median), > 2 (mode), 0 (no); default = 2
                                 norm = 2,        # data normalization: 1 (0.75-0.25 quantile), >1 (var=1), 0 (no); default = 1
                                 lap = 1.0,       # minimal value of the variational parameter, default = 1
                                 nL = 1           # maximal number of biclusters at which a row element can participate; default = 0 (no limit)
)

fabia.sim.scaledData_sparse <- FABIAobject.sim.trained_sparse@X; 
fabia.sim.scaledData_sparse <- as.matrix(fabia.sim.scaledData_sparse)

# samples
samples_info.sim_sparse = data.frame(rownames(t(fabia.sim.scaledData_sparse))); colnames(samples_info.sim_sparse) <- c("samples")

# fabia scores
data_factors_fabia.sim_sparse = data.frame(t(FABIAobject.sim.trained_sparse@Z)) 
data_factors_fabia.sim_sparse$samples = rownames(data_factors_fabia.sim_sparse)
colnames(data_factors_fabia.sim_sparse) <- c("score_FABIA1", "score_FABIA2", "sample")
df_fabia_factor_sparse <- data_factors_fabia.sim_sparse
df_fabia_factor_sparse$sample_id <- as.numeric(gsub("sample_", "", df_fabia_factor_sparse$sample))
df_fabia_factor_sparse = merge(df_fabia_factor_sparse, group_data, by = "sample_id")

# sort in ascending order by sample_id
df_fabia_factor_sparse <- df_fabia_factor_sparse[order(df_fabia_factor_sparse$sample_id), ]

# Factor 1
plot_factor_with_indices(
  scores     = (-1)*(df_fabia_factor_sparse$score_FABIA2),
  signal_idx = simulatedData$iteration_1$indices_samples[[2]],
  factor_id  = 1,
  shuffle    = "signal",
  x_breaks   = seq(0, 100, 20),
  x_limits   = c(0, 100)
)

# Factor 2
plot_factor_with_indices(
  scores     = (df_fabia_factor$score_FABIA1),
  signal_idx = simulatedData$iteration_1$indices_samples[[2]],
  factor_id  = 2,
  shuffle    = "signal",
  x_breaks   = seq(0, 100, 20),
  x_limits   = c(0, 100)
)

# Weights/Loading
loading_simfabia = data.frame(FABIAobject.sim.trained@L)
loading_simfabia$feature = rownames(loading_simfabia)
colnames(loading_simfabia) <- c("loading_FABIA1", "loading_FABIA2","features")
features_info.sim = data.frame(rownames(fabia.sim.scaledData)) 
colnames(features_info.sim) <- c("features") 
features_info.sim$ID <- as.integer(row.names(features_info.sim))
features_info.sim$view <- NA  # Create a new column with NAs
features_info.sim$view[1:4000] <- "Omic 1"  # Assign 'x' to rows 1 to 50
features_info.sim$view[4001:7000] <- 'Omic 2'  # Assign 'y' to rows 51 to 100

loading_simfabia = merge(loading_simfabia, features_info.sim, by = "features")

omic1_simfabia <- loading_simfabia %>%
  filter(view == 'Omic 1') %>%
  select(features, view, loading_FABIA1,loading_FABIA2, ID)
omic1_simfabia <- omic1_simfabia[order(omic1_simfabia$ID), ]

omic2_simfabia <- loading_simfabia %>%
  filter(view == 'Omic 2') %>%
  select(features, view, loading_FABIA1,loading_FABIA2, ID)
omic2_simfabia <- omic2_simfabia[order(omic2_simfabia$ID), ]

# Dataset 1
plot_overlay_violin_scatter(
  values     = -1*(omic1_simfabia$loading_FABIA2),
  signal_idx = simulatedData$iteration_1$indices_features.1[[1]],
  factor_id  = 1,
  kind       = "Features",
  shuffle    = "none",
  x_limits   = c(0, 4000),
  x_breaks   = seq(0, 4000, 1000),
  violin_alpha = 0.30
)

plot_overlay_violin_scatter(
  values     = omic1_simfabia$loading_FABIA1,
  signal_idx = simulatedData$iteration_1$indices_features.1[[2]],
  factor_id  = 1,
  kind       = "Features",
  shuffle    = "none",
  x_limits   = c(0, 4000),
  x_breaks   = seq(0, 4000, 1000),
  violin_alpha = 0.30
)

# Dataset 2
plot_overlay_violin_scatter(
  values     = (-1)*(omic2_simfabia$loading_FABIA2),
  signal_idx = simulatedData$iteration_1$indices_features.2[[1]],
  factor_id  = 1,
  kind       = "Features",
  shuffle    = "none",
  x_limits   = c(0, 3000),
  x_breaks   = seq(0, 3000, 1000),
  violin_alpha = 0.30
)

plot_overlay_violin_scatter(
  values     = omic2_simfabia$loading_FABIA2,
  signal_idx = NULL,
  factor_id  = 1,
  kind       = "Features",
  shuffle    = "none",
  y_limits   = c(-0.03, 0.1),
  x_limits   = c(0, 3000),
  x_breaks   = seq(0, 3000, 1000),
  violin_alpha = 0.30
)

## -------------- MFA ----------------------
set.seed(123)

# first_omic: 100 x 4000 (samples x features)
# second_omic: 100 x 3000

library(FactoMineR)

# Name the views and keep them samples x features
V <- list(omic1 = first_omic,
          omic2 = second_omic)

# add tiny jitter to zero-variance columns so nothing gets dropped
jitter_cols <- function(M, eps = 1e-8) {
  sds <- apply(M, 2, sd)
  bad <- which(!is.finite(sds) | sds == 0)
  if (length(bad)) {
    M[, bad] <- M[, bad, drop = FALSE] +
      matrix(rnorm(nrow(M) * length(bad), 0, eps), nrow(M))
  }
  M
}
V <- lapply(V, jitter_cols)

# bind blocks and define group sizes
X_mfa <- do.call(cbind, V)                          # 100 x 7000
group_sizes <- vapply(V, ncol, 1L)                  # c(4000, 3000)
stopifnot(sum(group_sizes) == ncol(X_mfa))

# FactoMineR wants a data.frame
X_df <- as.data.frame(X_mfa)
if (!is.null(rownames(X_mfa))) rownames(X_df) <- rownames(X_mfa)

k <- 2  # number of factors you want

fit <- MFA(
  X_df,
  group = group_sizes,
  type  = rep("s", length(group_sizes)),   # <-- quantitative blocks
  ncp   = k,
  graph = FALSE
)

# Scores (samples x factors)
Z <- as.matrix(fit$ind$coord)
if (ncol(Z) > 0) colnames(Z) <- paste0("F", seq_len(ncol(Z)))

# Loadings (features x factors) – use quanti.var for numeric blocks
Lq <- fit$quanti.var$coord
dim_cols <- grep("^Dim\\.", colnames(Lq))
L <- as.matrix(Lq[, dim_cols, drop = FALSE])
colnames(L) <- paste0("F", seq_len(ncol(L)))
rownames(L) <- unlist(Map(function(nm, M) paste0(nm, "::", colnames(M)),
                          names(V), V), use.names = FALSE)

results$MFA <- list(scores = Z, loadings = L, meta = list(n_factors = ncol(Z)))

F <- results$MFA$scores   # features x factors
rn <- rownames(F)

F_omic1 <- F[grepl("^sample_", rn), , drop = FALSE]

# Factor 1
plot_factor_with_indices(
  scores     = F_omic1[, "F1"],
  signal_idx = simulatedData$iteration_1$indices_samples[[1]],
  factor_id  = 1,
  shuffle    = "signal",
  x_breaks   = seq(0, 100, 20),
  x_limits   = c(0, 100)
)

# Factor 2
plot_factor_with_indices(
  scores     = F_omic1[, "F2"],
  signal_idx = simulatedData$iteration_1$indices_samples[[2]],
  factor_id  = 2,
  shuffle    = "signal",
  x_breaks   = seq(0, 100, 20),
  x_limits   = c(0, 100)
)

L <- results$MFA$loadings   # features x factors
rn <- rownames(L)

L_omic1 <- L[grepl("^omic1::", rn), , drop = FALSE]
L_omic2 <- L[grepl("^omic2::", rn), , drop = FALSE]

# omic1 features (ncol(first_omic) = 4000)
plot_overlay_violin_scatter(
  values     = L_omic1[, "F1"],
  signal_idx = simulatedData$iteration_1$indices_features.1[[1]],  # or NULL
  factor_id  = 1, kind = "Features",
  x_limits   = c(0, ncol(first_omic)),
  x_breaks   = seq(0, ncol(first_omic), 1000),
  violin_alpha = 0.30
)

# omic1 features (ncol(first_omic) = 4000)
plot_overlay_violin_scatter(
  values     = L_omic1[, "F2"],
  signal_idx = simulatedData$iteration_1$indices_features.1[[2]],  # or NULL
  factor_id  = 1, kind = "Features",
  x_limits   = c(0, ncol(first_omic)),
  x_breaks   = seq(0, ncol(first_omic), 1000),
  violin_alpha = 0.30
)

# omic2 features (ncol(second_omic) = 3000)
plot_overlay_violin_scatter(
  values     = L_omic2[, "F1"],
  signal_idx = simulatedData$iteration_1$indices_features.2[[1]],  # or NULL
  factor_id  = 1, kind = "Features",
  x_limits   = c(0, ncol(second_omic)),
  x_breaks   = seq(0, ncol(second_omic), 1000),
  violin_alpha = 0.30
)

plot_overlay_violin_scatter(
  values     = L_omic2[, "F2"],
  signal_idx = NULL,
  factor_id  = 1, kind = "Features",
  x_limits   = c(0, ncol(second_omic)),
  x_breaks   = seq(0, ncol(second_omic), 1000),
  violin_alpha = 0.30
)

## -------------------- GFA-----------------
# -------------
# GFA (multiview) — single run, fixed labeling & single-matrix W support
# -------------

library(GFA)

# --- 1) Assemble views: samples x features ---
V <- list(omic1 = first_omic,
          omic2 = second_omic)

# --- 2) Coerce to numeric matrices and add names ---
V <- lapply(V, function(M) {
  M <- as.matrix(M)
  storage.mode(M) <- "double"
  M
})

# identical, non-NULL sample IDs for all views
N <- nrow(V[[1]])
stopifnot(all(vapply(V, nrow, 1L) == N))
sample_ids <- rownames(V[[1]])
if (is.null(sample_ids)) sample_ids <- sprintf("sample_%03d", seq_len(N))
for (i in seq_along(V)) rownames(V[[i]]) <- sample_ids

# feature names (needed for loadings later)
for (i in seq_along(V)) {
  if (is.null(colnames(V[[i]]))) {
    colnames(V[[i]]) <- sprintf("%s_f%04d", names(V)[i], seq_len(ncol(V[[i]])))
  }
}

# --- 3) Stabilize zero-variance cols + scale (keeps dimnames) ---
fix_and_scale <- function(M, eps = 1e-8) {
  sds <- apply(M, 2, sd)
  bad <- which(!is.finite(sds) | sds == 0)
  if (length(bad)) {
    M[, bad] <- M[, bad, drop = FALSE] +
      matrix(rnorm(nrow(M) * length(bad), 0, eps), nrow(M))
  }
  # scale() preserves row/col names
  S <- scale(M, center = TRUE, scale = TRUE)
  # guard against possible all-NA columns after scaling (rare)
  S[, apply(S, 2, function(x) all(is.finite(x))), drop = FALSE]
}
Vz <- lapply(V, fix_and_scale)

# --- 4) Final sanity before GFA (these must be FALSE/TRUE as indicated) ---
stopifnot(!any(vapply(Vz, function(M) is.null(rownames(M)), logical(1))))     # FALSE
stopifnot(length(unique(lapply(Vz, rownames))) == 1)                           # TRUE
stopifnot(!any(vapply(Vz, function(M) is.null(colnames(M)), logical(1))))     # FALSE

# --- 5) Fit GFA ---
opts <- GFA::getDefaultOpts()
opts$iter.burnin <- 500
opts$iter.max    <- 1000
k <- 2

fit <- GFA::gfa(Vz, K = k, opts = opts)

# --- 6) Extract scores (samples x k) ---
Z <- if (is.list(fit$X)) as.matrix(fit$X[[1]]) else as.matrix(fit$X)
if (is.null(Z)) Z <- matrix(0, nrow = N, ncol = 0, dimnames = list(sample_ids, NULL))
if (NCOL(Z) > 0 && is.null(colnames(Z))) colnames(Z) <- paste0("F", seq_len(NCOL(Z)))

# --- 7) Extract loadings (features x k), stacked per view ---
build_block <- function(view_name, Wi, Vview) {
  Wi <- as.matrix(Wi)
  # orient to features x k
  if (nrow(Wi) == ncol(Vview)) {
    Wi <- t(Wi)                  # was k x features -> features x k
  } else if (ncol(Wi) != ncol(Vview)) {
    # unexpected; return empty block with correct rownames
    Wi <- matrix(0, nrow = ncol(Vview), ncol = 0)
  }
  rownames(Wi) <- paste0(view_name, "::", colnames(Vview))
  if (NCOL(Wi) > 0 && is.null(colnames(Wi))) colnames(Wi) <- paste0("F", seq_len(NCOL(Wi)))
  Wi
}

if (is.list(fit$W) && length(fit$W) == length(Vz)) {
  W_blocks <- Map(build_block, names(Vz), fit$W, Vz)
  L <- do.call(rbind, W_blocks)
} else if (is.matrix(fit$W)) {
  feat_per_view <- vapply(Vz, ncol, 1L)
  cuts <- c(0, cumsum(feat_per_view))
  Wi_all <- as.matrix(fit$W)
  if (nrow(Wi_all) == sum(feat_per_view)) {
    blocks <- lapply(seq_along(Vz), function(i) {
      rows <- (cuts[i] + 1):cuts[i + 1]
      Wi <- Wi_all[rows, , drop = FALSE]
      rownames(Wi) <- paste0(names(Vz)[i], "::", colnames(Vz[[i]]))
      if (NCOL(Wi) > 0 && is.null(colnames(Wi))) colnames(Wi) <- paste0("F", seq_len(NCOL(Wi)))
      Wi
    })
    L <- do.call(rbind, blocks)
  } else if (ncol(Wi_all) == sum(feat_per_view)) {
    Wi_all <- t(Wi_all)
    blocks <- lapply(seq_along(Vz), function(i) {
      rows <- (cuts[i] + 1):cuts[i + 1]
      Wi <- Wi_all[rows, , drop = FALSE]
      rownames(Wi) <- paste0(names(Vz)[i], "::", colnames(Vz[[i]]))
      if (NCOL(Wi) > 0 && is.null(colnames(Wi))) colnames(Wi) <- paste0("F", seq_len(NCOL(Wi)))
      Wi
    })
    L <- do.call(rbind, blocks)
  } else {
    # fallback: labeled empty
    L <- do.call(rbind, lapply(seq_along(Vz), function(i) {
      M <- matrix(0, nrow = ncol(Vz[[i]]), ncol = 0)
      rownames(M) <- paste0(names(Vz)[i], "::", colnames(Vz[[i]]))
      M
    }))
  }
} else {
  L <- do.call(rbind, lapply(seq_along(Vz), function(i) {
    M <- matrix(0, nrow = ncol(Vz[[i]]), ncol = 0)
    rownames(M) <- paste0(names(Vz)[i], "::", colnames(Vz[[i]]))
    M
  }))
}

results$GFA <- list(
  scores = Z,
  loadings = L,
  meta = list(n_factors = NCOL(Z), iter.max = opts$iter.max,
              iter.burnin = opts$iter.burnin, K = k)
)


L <- results$GFA$loadings  # features x factors (stacked)
idx <- list(
  omic1 = grepl("^omic1::", rownames(L)),
  omic2 = grepl("^omic2::", rownames(L))
)

# L2 norm per view per factor
activity <- sapply(idx, function(ii) colSums(L[ii, , drop = FALSE]^2)^0.5)
activity <- t(activity)  # rows: views, cols: factors
round(activity, 3)

Z_mfa <- results$MFA$scores
Z_gfa <- results$GFA$scores
k <- min(ncol(Z_mfa), ncol(Z_gfa))
C <- cor(Z_mfa[, 1:k, drop=FALSE], Z_gfa[, 1:k, drop=FALSE])

# best one-to-one match
if (!requireNamespace("clue", quietly = TRUE)) install.packages("clue")
perm <- clue::solve_LSAP(abs(C), maximum = TRUE)

# reorder & sign-flip GFA to match MFA
flip <- sign(diag(C[, perm]))
Z_gfa_aligned <- sweep(Z_gfa[, perm, drop=FALSE], 2, flip, `*`)
L_gfa_aligned <- sweep(results$GFA$loadings[, perm, drop=FALSE], 2, flip, `*`)

# Factor 1
plot_factor_with_indices(
  scores     = Z_gfa_aligned[, "F1"],
  signal_idx = simulatedData$iteration_1$indices_samples[[1]],
  factor_id  = 1,
  shuffle    = "signal",
  x_breaks   = seq(0, 100, 20),
  x_limits   = c(0, 100)
)

# Factor 2
plot_factor_with_indices(
  scores     = Z_gfa_aligned[, "F2"],
  signal_idx = simulatedData$iteration_1$indices_samples[[2]],
  factor_id  = 2,
  shuffle    = "signal",
  x_breaks   = seq(0, 100, 20),
  x_limits   = c(0, 100)
)

# omic1 factor 1
plot_overlay_violin_scatter(
  values     = L_gfa_aligned[grepl("^omic1::", rownames(L_gfa_aligned)), "F1"],
  signal_idx = simulatedData$iteration_1$indices_features.1[[1]],
  factor_id  = 1, kind = "Features",
  x_limits   = c(0, ncol(first_omic)),
  x_breaks   = seq(0, ncol(first_omic), 1000),
  y_limits   = NULL,
  violin_alpha = 0.30
)

# omic1 factor 1
plot_overlay_violin_scatter(
  values     = L_gfa_aligned[grepl("^omic1::", rownames(L_gfa_aligned)), "F2"],
  signal_idx = simulatedData$iteration_1$indices_features.1[[2]],
  factor_id  = 1, kind = "Features",
  x_limits   = c(0, ncol(first_omic)),
  x_breaks   = seq(0, ncol(first_omic), 1000),
  y_limits   = NULL,
  violin_alpha = 0.30
)

# omic2 factor 1
plot_overlay_violin_scatter(
  values     = L_gfa_aligned[grepl("^omic2::", rownames(L_gfa_aligned)), "F1"],
  signal_idx = simulatedData$iteration_1$indices_features.2[[1]],
  factor_id  = 1, kind = "Features",
  x_limits   = c(0, ncol(second_omic)),
  x_breaks   = seq(0, ncol(second_omic), 1000),
  violin_alpha = 0.30
)

# omic2 factor 1
plot_overlay_violin_scatter(
  values     = L_gfa_aligned[grepl("^omic2::", rownames(L_gfa_aligned)), "F2"],
  signal_idx = NULL,
  factor_id  = 1, kind = "Features",
  x_limits   = c(0, ncol(second_omic)),
  x_breaks   = seq(0, ncol(second_omic), 1000),
  violin_alpha = 0.30
)


# -------------- Comparison Plots (Factor & Loading) ----------------------
# Scores
simMOFAscores <- factor_mofa2
simFABIAscores <- df_fabia_factor
simMFAscores <- results$MFA$scores
simGFAscores <- results$GFA$scores

# Extract
simMFAscores <- data.frame(results$MFA$scores)
simGFAscores <- data.frame(results$GFA$scores)

# Rename MFA scores: F1, F2, ... -> score_MFA1, score_MFA2, ...
if (!is.null(simMFAscores) && NCOL(simMFAscores) > 0) {
  colnames(simMFAscores) <- paste0("score_MFA", seq_len(NCOL(simMFAscores)))
  simMFAscores$sample <- rownames(simMFAscores)
}


# Rename GFA scores: F1, F2, ... -> score_GFA1, score_GFA2, ...
if (!is.null(simGFAscores) && NCOL(simGFAscores) > 0) {
  colnames(simGFAscores) <- paste0("score_GFA", seq_len(NCOL(simGFAscores)))
  simGFAscores$sample <- rownames(simGFAscores)
}


# Rename GT scores: alpha1, alpha2, ... -> score_GT1, score_GT2, ...
simGTscores <- data.frame(simulatedData$iteration_1$factor_scores)
if (!is.null(simGTscores) && NCOL(simGTscores) > 0) {
  colnames(simGTscores) <- paste0("score_GT", seq_len(NCOL(simGTscores)))
  simGTscores$sample <- paste0("sample_", seq_len(NROW(simGTscores)))
}

# ----------------- merge data using common column ------------
merge_many <- function(dfs, by, all = TRUE) {
  # dfs: a list of data.frames (or tibbles)
  # by:  name of the column to merge on (character)
  # all: whether to keep all rows (like full join). If FALSE, keeps only common rows.
  
  if (!is.list(dfs)) {
    stop("`dfs` must be a list of data frames.")
  }
  if (length(dfs) == 0) {
    stop("`dfs` is empty.")
  }
  
  # check column exists in all
  missing_col <- vapply(dfs, function(x) !(by %in% names(x)), logical(1))
  if (any(missing_col)) {
    stop(
      "These data frames don't have the merge column: ",
      paste0(which(missing_col), collapse = ", ")
    )
  }
  
  # start merging
  out <- dfs[[1]]
  if (length(dfs) > 1) {
    for (i in 2:length(dfs)) {
      out <- merge(out, dfs[[i]], by = by, all = all)
    }
  }
  out
}

simFactorscores <- merge_many(list(simFABIAscores, 
                               simMOFAscores,
                               simMFAscores,
                               simGFAscores,
                               simGTscores), by = "sample", all = FALSE)

# Loading

# ---------------------------
# Factorization mega-pipeline
# ---------------------------
# Input data:
# - X can be either:
#     * a single numeric matrix (features x samples)
#     * a named list of numeric matrices (each features x samples). All must share samples.
# - 'methods' any subset of c("FABIA","MOFA","GFA","MFA")
#
# Returns a named list with one entry per method:
#   $<METHOD>$scores   : matrix (samples x k)
#   $<METHOD>$loadings : matrix (features(±view) x k)  (row names are features; if multiple views, rows are stacked with "view::feature")
#   $<METHOD>$meta     : list with method-specific bits
#
# You can set k (target #factors), scaling, and MFA group sizes.
#
# Required CRAN packages:
#   fabia, MOFA2, GFA, FactoMineR, Matrix
# --------------------------------------

# ---------------------------
# Helper: build comparison object (loadings per method, same feature order)
# ---------------------------
build_loading_comparison <- function(results) {
  stopifnot(length(results) >= 1)
  # union of feature names across methods
  feat_all <- unique(unlist(lapply(results, function(r) rownames(r$loadings))))
  out <- lapply(results, function(r) {
    L <- r$loadings
    # align to union (pad missing with 0)
    M <- matrix(0, nrow = length(feat_all), ncol = ncol(L),
                dimnames = list(feat_all, colnames(L)))
    M[rownames(L), ] <- L
    M
  })
  return(out)  # list: method -> (features x k)
}

run_factor_pipeline <- function(
    X,
    methods = c("FABIA","MOFA","GFA","MFA"),
    k = 10,
    scale_features = TRUE,
    mfa_group_sizes = NULL,
    verbose = TRUE
) {
  .stop_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE)) stop("Package '", p, "' is required.")
  .stack_views <- function(lst) {
    do.call(rbind, lapply(names(lst), function(v) {
      M <- lst[[v]]
      rn <- rownames(M); if (is.null(rn)) rn <- paste0("f", seq_len(nrow(M)))
      rownames(M) <- paste0(v, "::", rn)
      M
    }))
  }
  .as_list_views <- function(X) if (is.list(X)) X else list(omic1 = X)
  .nzcol_mask <- function(M) {
    # keep columns with non-zero variance (ignoring NAs)
    if (!is.matrix(M)) M <- as.matrix(M)
    if (ncol(M) == 0) return(rep(FALSE, 0))
    apply(M, 2, function(v) {
      vv <- v[is.finite(v)]
      if (!length(vv)) return(FALSE)
      sd(vv) > 0
    })
  }
  
  methods <- toupper(methods)
  
  # ---- input checks & harmonize ----
  if (!(is.matrix(X) || is.list(X))) stop("X must be a matrix or a list of matrices.")
  V <- .as_list_views(X)
  sample_names <- colnames(V[[1]])
  if (is.null(sample_names)) stop("All matrices must have column (sample) names.")
  for (nm in names(V)) {
    v <- V[[nm]]
    if (is.null(colnames(v))) stop("All matrices must have column (sample) names.")
    if (!setequal(colnames(v), sample_names))
      stop("All views must contain the same sample names.")
  }
  # reorder all views to same sample order
  V <- lapply(V, function(m) m[, sample_names, drop = FALSE])
  
  # optional scaling (feature-wise)
  if (scale_features) {
    V <- lapply(V, function(m) {
      m[] <- scale(t(m)) |> t()
      m[is.na(m)] <- 0
      m
    })
  }
  
  results <- list()
  
  # ----------------
  # FABIA (concat)
  # ----------------
  if ("FABIA" %in% methods) {
    .stop_pkg("fabia")
    if (verbose) message("Running FABIA ...")
    M_fabia <- .stack_views(V)   # features x samples
    p <- min(nrow(M_fabia), k)
    set.seed(21)#991489
    fit <- fabia(X = M_fabia, p = p, alpha = 0.01, cyc = 1000, spz = 0.5, lap = 1.0, nL = 2)
    # In fabia objects, loadings ~ L (features x k), factors/scores ~ Z (samples x k)
    
    L <- tryCatch(fit@L, error = function(e) NULL)
    Z <- tryCatch(fit@Z, error = function(e) NULL)
    
    # loadings: features x k
    if (!is.null(L)) {
      if (nrow(L) != nrow(M_fabia) && ncol(L) == nrow(M_fabia)) L <- t(L)
      rownames(L) <- rownames(M_fabia)
      if (ncol(L) > 0) colnames(L) <- paste0("F", seq_len(ncol(L)))
    } else {
      L <- matrix(0, nrow = nrow(M_fabia), ncol = 0, dimnames = list(rownames(M_fabia), NULL))
    }
    # scores: samples x k
    if (!is.null(Z)) {
      if (nrow(Z) != ncol(M_fabia) && ncol(Z) == ncol(M_fabia)) Z <- t(Z)
      rownames(Z) <- colnames(M_fabia)
      if (ncol(Z) > 0) colnames(Z) <- paste0("F", seq_len(ncol(Z)))
    } else {
      Z <- matrix(0, nrow = ncol(M_fabia), ncol = 0, dimnames = list(colnames(M_fabia), NULL))
    }
    results$FABIA <- list(scores = Z, loadings = L, meta = list(n_factors = ncol(Z)))
  }
  
  # -------------
  # MOFA2 (train with all features, let MOFA handle scaling)
  # -------------
  # if ("MOFA" %in% methods) {
  #   .stop_pkg("MOFA2")
  #   if (verbose) message("Running MOFA2 ...")
  #   library(MOFA2)
  #   
  #   # Ensure input matrices are numeric with rownames
  #   V_full <- lapply(V, function(m) {
  #     M <- as.matrix(m)
  #     storage.mode(M) <- "double"
  #     if (is.null(rownames(M))) rownames(M) <- paste0("f", seq_len(nrow(M)))
  #     M
  #   })
  #   
  #   # Build MOFA object
  #   mofa_obj <- MOFA2::create_mofa(V_full)
  #   
  #   # --- Data options: let MOFA do the scaling
  #   data_opts <- MOFA2::get_default_data_options(mofa_obj)
  #   data_opts$scale_views  <- TRUE
  #   data_opts$scale_groups <- TRUE
  #   
  #   # --- Model options
  #   model_opts <- MOFA2::get_default_model_options(mofa_obj)
  #   model_opts$num_factors <- k
  #   
  #   # --- Training options
  #   train_opts <- MOFA2::get_default_training_options(mofa_obj)
  #   #train_opts$verbose <- isTRUE(verbose)
  #   train_opts$seed    <- 21#12367
  #   train_opts$convergence_mode <- "slow"
  #   train_opts$maxiter <- 2000
  #   
  #   # --- Prepare and run MOFA
  #   mofa_obj <- MOFA2::prepare_mofa(
  #     object = mofa_obj,
  #     data_options     = data_opts,
  #     model_options    = model_opts,
  #     training_options = train_opts
  #   )
  # 
  #   outfile <- paste0(getwd(),"model_raw_pipeline.hdf5")
  #   #file.path(
  #   #   tempdir(),
  #   #   paste0("mofa_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".hdf5")
  #   # )
  #   
  #   fit_trained <- MOFA2::run_mofa(mofa_obj, outfile, use_basilisk = TRUE)
  #   
  #   # --- Load model (keep inactive factors)
  #   fit <- tryCatch({
  #     MOFA2::load_model(outfile, remove_inactive_factors = FALSE)
  #   }, error = function(e) {
  #     if (methods::is(fit_trained, "MOFA")) {
  #       warning("MOFA2: load_model() failed; using in-memory model. Inactive factors may be pruned.")
  #       fit_trained
  #     } else {
  #       stop("MOFA2: load_model() failed and run_mofa() did not return a MOFA object: ",
  #            conditionMessage(e))
  #     }
  #   })
  #   
  #   # --- Extract scores
  #   Z_list <- MOFA2::get_factors(fit, factors = "all", as.data.frame = FALSE)
  #   Z <- as.matrix(Z_list[[1]])
  #   if (!is.null(Z) && ncol(Z) > 0) {
  #     colnames(Z) <- paste0("F", seq_len(ncol(Z)))
  #   } else {
  #     Z <- matrix(0, nrow = length(sample_names), ncol = 0,
  #                 dimnames = list(sample_names, NULL))
  #   }
  #   
  #   # --- Extract loadings
  #   W_list <- MOFA2::get_weights(fit, views = "all", factors = "all", as.data.frame = FALSE)
  #   L <- .stack_views(W_list)
  #   if (ncol(L) > 0) colnames(L) <- paste0("F", seq_len(ncol(L)))
  #   
  #   results$MOFA <- list(
  #     scores   = Z,
  #     loadings = L,
  #     meta     = list(
  #       n_factors = ncol(Z),
  #       trained_hdf5 = outfile,
  #       scaling_in_mofa = TRUE
  #     )
  #   )
  # }
  if ("MOFA" %in% methods) {
    .stop_pkg("MOFA2")
    if (verbose) message("Running MOFA2 ...")
    library(MOFA2)
    
    # Ensure input matrices are numeric with rownames
    V_full <- lapply(V, function(m) {
      M <- as.matrix(m)
      storage.mode(M) <- "double"
      if (is.null(rownames(M))) rownames(M) <- paste0("f", seq_len(nrow(M)))
      M
    })
    
    # Build MOFA object
    mofa_obj <- MOFA2::create_mofa(V_full)
    
    # --- Data options: let MOFA do the scaling
    data_opts <- MOFA2::get_default_data_options(mofa_obj)
    data_opts$scale_views  <- TRUE
    data_opts$scale_groups <- TRUE
    
    # --- Model options
    model_opts <- MOFA2::get_default_model_options(mofa_obj)
    model_opts$num_factors <- k
    
    # --- Training options
    train_opts <- MOFA2::get_default_training_options(mofa_obj)
    train_opts$seed    <- 21
    train_opts$convergence_mode <- "slow"
    train_opts$maxiter <- 2000
    
    # --- Prepare and run MOFA
    mofa_obj <- MOFA2::prepare_mofa(
      object = mofa_obj,
      data_options     = data_opts,
      model_options    = model_opts,
      training_options = train_opts
    )
    
    # Safer path construction
    outfile <- file.path(getwd(), "model_raw_pipeline.hdf5")
    
    fit_trained <- MOFA2::run_mofa(mofa_obj, outfile, use_basilisk = TRUE)
    
    # --- Load model (keep inactive factors)
    fit <- tryCatch({
      MOFA2::load_model(outfile, remove_inactive_factors = FALSE)
    }, error = function(e) {
      if (methods::is(fit_trained, "MOFA")) {
        warning("MOFA2: load_model() failed; using in-memory model. Inactive factors may be pruned.")
        fit_trained
      } else {
        stop("MOFA2: load_model() failed and run_mofa() did not return a MOFA object: ",
             conditionMessage(e))
      }
    })
    
    # --- Extract scores
    Z_list <- MOFA2::get_factors(fit, factors = "all", as.data.frame = FALSE)
    Z <- as.matrix(Z_list[[1]])
    if (!is.null(Z) && ncol(Z) > 0) {
      colnames(Z) <- paste0("F", seq_len(ncol(Z)))
    } else {
      Z <- matrix(0, nrow = length(sample_names), ncol = 0,
                  dimnames = list(sample_names, NULL))
    }
    
    # --- Extract loadings
    W_list <- MOFA2::get_weights(fit, views = "all", factors = "all", as.data.frame = FALSE)
    L <- .stack_views(W_list)
    if (ncol(L) > 0) colnames(L) <- paste0("F", seq_len(ncol(L)))
    
    # --- New: variance explained (plot + table)
    p_varExp <- tryCatch(
      MOFA2::plot_variance_explained(fit, max_r2 = 15),
      error = function(e) { warning("plot_variance_explained failed: ", conditionMessage(e)); NULL }
    )
    varExp_tbl <- tryCatch(
      MOFA2::calculate_variance_explained(fit),
      error = function(e) NULL
    )
    
    # --- New: factor-by-group boxplots
    p_factorByGroup <- tryCatch({
      if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 not installed")
      # Get sample metadata from the fitted object if available
      smd <- tryCatch(MOFA2::samples_metadata(fit), error = function(e) NULL)
      # Build a minimal metadata if missing
      if (is.null(smd) || nrow(smd) == 0) {
        smd <- data.frame(sample = rownames(Z), Group = "All", row.names = rownames(Z))
      } else {
        # ensure rownames and a Group column
        if (is.null(rownames(smd)) && "sample" %in% names(smd)) {
          rownames(smd) <- smd$sample
        }
        smd <- smd[rownames(Z), , drop = FALSE]
        if (!("Group" %in% colnames(smd))) {
          # Try to pick a categorical column as Group; else fallback to "All"
          cat_cols <- names(which(vapply(smd, function(x) is.factor(x) || is.character(x), logical(1))))
          if (length(cat_cols) > 0) {
            smd$Group <- smd[[cat_cols[1]]]
          } else {
            smd$Group <- "All"
          }
        }
      }
      
      # Long table of factor scores
      df_long <- utils::stack(as.data.frame(Z))
      df_long$sample <- rep(rownames(Z), times = ncol(Z))
      names(df_long) <- c("score", "Factor", "sample")
      df_long$Factor <- as.character(df_long$Factor)
      df_long$Group  <- smd$Group[match(df_long$sample, rownames(smd))]
      
      ggplot2::ggplot(df_long, ggplot2::aes(x = Factor, y = score, fill = Group)) +
        ggplot2::geom_boxplot(outlier.size = 0.5) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "MOFA factor scores by group",
          x = "Factor",
          y = "Score",
          fill = "Group"
        )
    }, error = function(e) { warning("factor-by-group plot failed: ", conditionMessage(e)); NULL })
    
    results$MOFA <- list(
      scores   = Z,
      loadings = L,
      meta     = list(
        n_factors = ncol(Z),
        trained_hdf5 = outfile,
        scaling_in_mofa = TRUE,
        variance_explained = varExp_tbl
      ),
      plots = list(
        variance_by_factor = p_varExp,      # ggplot (may be NULL on error)
        factors_by_group   = p_factorByGroup # ggplot (may be NULL on error)
      )
    )
  }
  
  # -------------
  # GFA (multiview) — single run, fixed labeling & single-matrix W support
  # -------------
  if ("GFA" %in% methods) {
    .stop_pkg("GFA")
    if (verbose) message("Running GFA ...")
    library(GFA)
    
    ## Ensure views have stable names: 'rna' and 'prot' for 2 views, else view1, view2, ...
    if (is.null(names(V)) || any(!nzchar(names(V)))) {
      names(V) <- if (length(V) == 2) c("rna", "prot") else paste0("view", seq_along(V))
    }
    
    ## GFA expects samples x features; keep feature names
    V_gfa <- lapply(V, t)                 # each: samples x features
    names(V_gfa) <- names(V)              # keep labels ('rna','prot')
    
    ## Options (yours)
    model_option <- GFA::getDefaultOpts()
    model_option$iter.burnin <- 500
    model_option$iter.max    <- 1000
    
    ## Normalization (as you had)
    norm <- GFA::normalizeData(V_gfa, type = "center")
    
    ## Fit once
    fit <- GFA::gfa(norm$train, K = k, opts = model_option)
    
    ## ---------- SCORES (samples x k) ----------
    Z <- {
      Zcand <- if (is.list(fit$X) && length(fit$X) >= 1) fit$X[[1]] else fit$X
      if (is.null(Zcand)) {
        matrix(0, nrow = nrow(V_gfa[[1]]), ncol = 0,
               dimnames = list(rownames(V_gfa[[1]]), NULL))
      } else {
        Zcand <- as.matrix(Zcand)
        if (is.null(rownames(Zcand))) rownames(Zcand) <- rownames(V_gfa[[1]])
        if (NCOL(Zcand) > 0 && is.null(colnames(Zcand))) colnames(Zcand) <- paste0("F", seq_len(NCOL(Zcand)))
        Zcand
      }
    }
    
    ## ---------- LOADINGS (features x k), with 'rna::' / 'prot::' prefixes ----------
    make_empty_block <- function(i) {
      nf <- ncol(V_gfa[[i]])
      feats <- colnames(V_gfa[[i]]); if (is.null(feats)) feats <- paste0("feature_", seq_len(nf))
      M <- matrix(0, nrow = nf, ncol = 0)
      rownames(M) <- paste0(names(V_gfa)[i], "::", feats)
      M
    }
    
    W_list <- list()
    
    if (is.list(fit$W) && length(fit$W) == length(V_gfa)) {
      ## Standard case: one W per view
      W_list <- lapply(seq_along(fit$W), function(i) {
        Wi <- as.matrix(fit$W[[i]])
        if (is.null(dim(Wi))) Wi <- matrix(Wi, nrow = 1)
        
        n_feat <- ncol(V_gfa[[i]])               # features in this view
        feats  <- colnames(V_gfa[[i]]); if (is.null(feats)) feats <- paste0("feature_", seq_len(n_feat))
        
        ## Orient to (features x k)
        if (nrow(Wi) == n_feat) {
          # ok
        } else if (ncol(Wi) == n_feat) {
          Wi <- t(Wi)
        } else if (n_feat == 1 && (nrow(Wi) == 1 || ncol(Wi) == 1)) {
          if (nrow(Wi) != 1) Wi <- t(Wi)
        } else {
          # unreconcilable for this view → empty, correctly labeled
          return(make_empty_block(i))
        }
        
        rownames(Wi) <- paste0(names(V_gfa)[i], "::", feats)
        if (NCOL(Wi) > 0 && is.null(colnames(Wi))) colnames(Wi) <- paste0("F", seq_len(NCOL(Wi)))
        Wi
      })
      
    } else if (is.matrix(fit$W)) {
      ## All views merged into a single W matrix
      Wi_all <- as.matrix(fit$W)
      
      view_feat   <- vapply(V_gfa, ncol, integer(1))   # features per view
      total_feat  <- sum(view_feat)
      
      ## Orient to (total_features x k)
      if (nrow(Wi_all) == total_feat) {
        # ok
      } else if (ncol(Wi_all) == total_feat) {
        Wi_all <- t(Wi_all)
      } else if (nrow(Wi_all) == ncol(Wi_all) && nrow(Wi_all) < total_feat) {
        warning(sprintf("GFA returned square W (%d x %d); cannot map to %d features. Returning empty blocks.",
                        nrow(Wi_all), ncol(Wi_all), total_feat))
        W_list <- lapply(seq_along(V_gfa), make_empty_block)
      } else {
        warning(sprintf("Unexpected W shape: %d x %d (expected total features %d). Returning empty blocks.",
                        nrow(Wi_all), ncol(Wi_all), total_feat))
        W_list <- lapply(seq_along(V_gfa), make_empty_block)
      }
      
      ## If Wi_all is oriented correctly, split rows into per-view blocks
      if (length(W_list) == 0) {
        cuts <- c(0, cumsum(view_feat))
        W_list <- vector("list", length(V_gfa))
        for (i in seq_along(V_gfa)) {
          rows  <- (cuts[i] + 1):cuts[i + 1]
          Wi    <- Wi_all[rows, , drop = FALSE]            # features_i x k
          feats <- colnames(V_gfa[[i]]); if (is.null(feats)) feats <- paste0("feature_", seq_len(nrow(Wi)))
          rownames(Wi) <- paste0(names(V_gfa)[i], "::", feats)
          if (NCOL(Wi) > 0 && is.null(colnames(Wi))) colnames(Wi) <- paste0("F", seq_len(NCOL(Wi)))
          W_list[[i]] <- Wi
        }
      }
      
    } else {
      ## Unknown W type → labeled empty blocks
      W_list <- lapply(seq_along(V_gfa), make_empty_block)
    }
    
    L <- if (length(W_list)) do.call(rbind, W_list) else NULL
    
    results$GFA <- list(
      scores   = Z,
      loadings = L,
      meta     = list(
        n_factors   = NCOL(Z),
        iter.max    = model_option$iter.max,
        iter.burnin = model_option$iter.burnin,
        K           = k
      )
    )
  }
  
  # -------------
  # MFA (FactoMineR) — filter zero-variance, NULL/empty-safe; robust loadings from quanti.var$coord
  # -------------
  # if ("MFA" %in% methods) {
  #   .stop_pkg("FactoMineR")
  #   if (verbose) message("Running MFA ...")
  #   library(FactoMineR)
  #   # samples x features per view (optionally filter zero-variance columns)
  #   Vt <- lapply(V, t)
  #   nz_mask <- function(M) if (ncol(M)) apply(M, 2, function(v) sd(v[is.finite(v)]) > 0) else logical(0)
  #   Vt <- lapply(Vt, function(m) { keep <- nz_mask(m); m[, keep, drop = FALSE] })
  #   
  #   # bind blocks; compute (or accept) group sizes
  #   X_mfa <- do.call(cbind, Vt)                  # samples x all_features
  #   if (is.null(mfa_group_sizes)) mfa_group_sizes <- vapply(Vt, ncol, 1L)
  #   stopifnot(sum(mfa_group_sizes) == ncol(X_mfa))
  #   
  #   # data.frame is friendlier to FactoMineR
  #   X_df <- as.data.frame(X_mfa)
  #   if (!is.null(rownames(X_mfa))) rownames(X_df) <- rownames(X_mfa)
  #   
  #   fit <- FactoMineR::MFA(
  #     X_df,
  #     group = mfa_group_sizes,
  #     type  = rep("c", length(mfa_group_sizes)),
  #     ncp   = k,
  #     graph = FALSE
  #   )
  #   
  #   ## Scores (samples x q)
  #   Z <- fit$ind$coord
  #   if (is.null(Z)) {
  #     Z <- matrix(0, nrow = nrow(X_df), ncol = 0,
  #                 dimnames = list(rownames(X_df), NULL))
  #     q <- 0L
  #   } else {
  #     Z <- as.matrix(Z); q <- ncol(Z)
  #     if (q > 0) colnames(Z) <- paste0("F", seq_len(q))
  #   }
  #   
  #   ## Loadings (features x q)
  #   # 1) Try global var coords; 2) else use quanti.var$coord; 3) else empty
  #   L <- fit$var$coord
  #   if (is.null(L)) {
  #     Lq <- fit$quanti.var$coord
  #     if (!is.null(Lq)) {
  #       # keep only dimensions (Dim.1, Dim.2, ...) and coerce to matrix
  #       dim_cols <- grep("^Dim\\.", colnames(Lq))
  #       L <- as.matrix(Lq[, dim_cols, drop = FALSE])
  #       # rename columns to F1..Fq (q may be <= k)
  #       if (ncol(L) > 0) colnames(L) <- paste0("F", seq_len(ncol(L)))
  #       
  #       # Reorder rows to match our concatenated block order
  #       # Build desired order of base variable names (per group, in the same order as X_mfa)
  #       want_base <- unlist(lapply(seq_along(Vt), function(i) colnames(Vt[[i]])), use.names = FALSE)
  #       
  #       # Row names in L often look like "group.var" or just "var"
  #       row_base <- sub("^.*?\\.", "", rownames(L))  # strip optional "group."
  #       # Try to match per-group to avoid cross-group duplicates
  #       offs <- c(0, cumsum(mfa_group_sizes)); offs <- offs[-length(offs)]
  #       idx <- integer(0)
  #       for (i in seq_along(Vt)) {
  #         g_vars <- colnames(Vt[[i]])
  #         # prefer matches where original rowname already equals "group.var"
  #         # but fall back to base-name match
  #         rnames <- rownames(L)
  #         cand_group <- paste0(names(V)[i], ".", g_vars)
  #         m <- match(cand_group, rnames)
  #         m[is.na(m)] <- match(g_vars[is.na(m)], row_base)
  #         idx <- c(idx, m)
  #       }
  #       # guard: drop any NAs (shouldn’t happen unless names collided)
  #       keep_ok <- !is.na(idx)
  #       L <- L[idx[keep_ok], , drop = FALSE]
  #       # final nice rownames with view prefixes
  #       nice_rn <- unlist(lapply(seq_along(Vt), function(i) paste0(names(V)[i], "::", colnames(Vt[[i]]))),
  #                         use.names = FALSE)
  #       rownames(L) <- nice_rn[keep_ok]
  #     }
  #   } else {
  #     L <- as.matrix(L)
  #     if (ncol(L) > 0) {
  #       # trim to q if needed
  #       if (!is.null(q) && q > 0 && ncol(L) > q) L <- L[, seq_len(q), drop = FALSE]
  #       colnames(L) <- paste0("F", seq_len(ncol(L)))
  #     }
  #     # add view prefixes in our concatenated order
  #     if (is.null(colnames(X_mfa))) colnames(X_mfa) <- paste0("v", seq_len(ncol(X_mfa)))
  #     block_names <- rep(names(V), times = mfa_group_sizes)
  #     rownames(L) <- paste0(block_names, "::", colnames(X_mfa))
  #   }
  #   
  #   # If still NULL/empty, return a dimensionally correct empty matrix
  #   if (is.null(L)) {
  #     rn <- unlist(lapply(seq_along(Vt), function(i)
  #       paste0(names(V)[i], "::", colnames(Vt[[i]]))), use.names = FALSE)
  #     L <- matrix(0, nrow = length(rn), ncol = 0, dimnames = list(rn, NULL))
  #   }
  #   
  #   results$MFA <- list(scores = Z, loadings = L, meta = list(n_factors = ncol(Z)))
  # }
  if ("MFA" %in% methods) {
    .stop_pkg("FactoMineR")
    if (verbose) message("Running MFA ...")
    library(FactoMineR)
    
    # transpose to samples x features per view
    Vt <- lapply(V, t)
    
    # add tiny jitter to any zero-variance columns so nothing is dropped
    jitter_cols <- function(M, eps = 1e-8) {
      if (ncol(M) == 0) return(M)
      sds <- apply(M, 2, sd, na.rm = TRUE)
      bad <- which(!is.finite(sds) | sds == 0)
      if (length(bad)) {
        M[, bad] <- M[, bad, drop = FALSE] +
          matrix(rnorm(nrow(M) * length(bad), 0, eps), nrow(M))
      }
      M
    }
    Vt <- lapply(Vt, jitter_cols)
    
    # bind blocks; compute (or accept) group sizes
    X_mfa <- do.call(cbind, Vt)
    if (is.null(mfa_group_sizes)) mfa_group_sizes <- vapply(Vt, ncol, 1L)
    stopifnot(sum(mfa_group_sizes) == ncol(X_mfa))
    
    # data.frame is friendlier to FactoMineR
    X_df <- as.data.frame(X_mfa)
    if (!is.null(rownames(X_mfa))) rownames(X_df) <- rownames(X_mfa)
    
    fit <- FactoMineR::MFA(
      X_df,
      group = mfa_group_sizes,
      type  = rep("c", length(mfa_group_sizes)),
      ncp   = k,
      graph = FALSE
    )
    
    ## Scores
    Z <- fit$ind$coord
    if (is.null(Z)) {
      Z <- matrix(0, nrow = nrow(X_df), ncol = 0,
                  dimnames = list(rownames(X_df), NULL))
      q <- 0L
    } else {
      Z <- as.matrix(Z); q <- ncol(Z)
      if (q > 0) colnames(Z) <- paste0("F", seq_len(q))
    }
    
    ## Loadings
    L <- fit$var$coord
    if (is.null(L)) {
      Lq <- fit$quanti.var$coord
      if (!is.null(Lq)) {
        dim_cols <- grep("^Dim\\.", colnames(Lq))
        L <- as.matrix(Lq[, dim_cols, drop = FALSE])
        if (ncol(L) > 0) colnames(L) <- paste0("F", seq_len(ncol(L)))
        
        # reorder and prefix
        nice_rn <- unlist(lapply(seq_along(Vt),
                                 function(i) paste0(names(V)[i], "::", colnames(Vt[[i]]))),
                          use.names = FALSE)
        rownames(L) <- nice_rn
      }
    } else {
      L <- as.matrix(L)
      if (ncol(L) > 0 && !is.null(q) && q > 0 && ncol(L) > q)
        L <- L[, seq_len(q), drop = FALSE]
      colnames(L) <- paste0("F", seq_len(ncol(L)))
      block_names <- rep(names(V), times = mfa_group_sizes)
      rownames(L) <- paste0(block_names, "::", colnames(X_mfa))
    }
    
    if (is.null(L)) {
      rn <- unlist(lapply(seq_along(Vt),
                          function(i) paste0(names(V)[i], "::", colnames(Vt[[i]]))),
                   use.names = FALSE)
      L <- matrix(0, nrow = length(rn), ncol = 0, dimnames = list(rn, NULL))
    }
    
    results$MFA <- list(scores = Z, loadings = L,
                        meta = list(n_factors = ncol(Z)))
  }
  
  
  if (verbose) message("Done. Methods fit: ", paste(names(results), collapse = ", "))
  results
}

#===========================================
# Create Benchmark Datasets
#===========================================
sim_benchmark_res <- run_factor_pipeline(
  list(omic1 = as.matrix(t(dataset_omic1)), omic2 = as.matrix(t(dataset_omic2))),
  methods = c("fabia","mofa","gfa","mfa"),
  k = 2,
  scale_features = TRUE
)

#============================================
# Merge 'scores' with metadata by sample for each method in a benchmark object
#============================================
merge_scores_with_metadata <- function(bench, metadata, sample_col = "sample",
                                       methods = intersect(names(bench), c("FABIA","MOFA","GFA","MFA","GT")),
                                       join = c("left","inner")) {
  join <- match.arg(join)
  
  # ensure metadata has the join key and unique rows
  meta <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (!(sample_col %in% names(meta))) {
    if (!is.null(rownames(meta))) {
      meta[[sample_col]] <- rownames(meta)
    } else {
      stop("`metadata` must have a sample column named ", sample_col, " or rownames to use as sample IDs.")
    }
  }
  meta[[sample_col]] <- as.character(meta[[sample_col]])
  meta <- distinct(meta, .data[[sample_col]], .keep_all = TRUE)
  
  by_method <- list()
  for (m in methods) {
    if (is.null(bench[[m]]$scores)) next
    
    S <- bench[[m]]$scores
    # coerce to data.frame and add 'sample' column from rownames if needed
    Sdf <- as.data.frame(S, check.names = FALSE, stringsAsFactors = FALSE)
    if (!("sample" %in% names(Sdf))) {
      rn <- rownames(Sdf)
      if (is.null(rn)) stop("Scores for ", m, " have no rownames and no 'sample' column.")
      Sdf$sample <- rn
    }
    Sdf$sample <- as.character(Sdf$sample)
    
    # join: scores$sample -> metadata[[sample_col]]
    by_map <- setNames(sample_col, "sample")
    out <- if (join == "left") {
      left_join(Sdf, meta, by = by_map)
    } else {
      inner_join(Sdf, meta, by = by_map)
    }
    
    # useful ordering: sample, factors..., metadata...
    fac_cols <- grep("^F\\d+(\\s*\\(|$)", names(Sdf), value = TRUE)
    front <- c("sample", fac_cols)
    out <- out[, c(intersect(front, names(out)), setdiff(names(out), front)), drop = FALSE]
    out$method <- m
    by_method[[m]] <- out
  }
  
  all_long <- bind_rows(by_method)
  list(by_method = by_method, all_long = all_long)
}

library(dplyr)

# Attach metadata to scores for every method in a benchmark object
attach_scores_metadata <- function(bench,
                                   metadata,
                                   sample_col = "sample",
                                   methods = intersect(names(bench), c("FABIA","MOFA","GFA","MFA","GT")),
                                   store_name = "scores_with_meta",
                                   replace_scores = FALSE) {
  # --- prep metadata ---
  meta <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (!(sample_col %in% names(meta))) {
    if (!is.null(rownames(meta))) {
      meta[[sample_col]] <- rownames(meta)
    } else {
      stop("`metadata` must have a '", sample_col, "' column or rownames.")
    }
  }
  meta[[sample_col]] <- as.character(meta[[sample_col]])
  meta <- distinct(meta, .data[[sample_col]], .keep_all = TRUE)
  
  for (m in methods) {
    if (is.null(bench[[m]]$scores)) next
    
    # coerce scores to data frame; pull sample IDs from rownames if needed
    Sdf <- as.data.frame(bench[[m]]$scores, check.names = FALSE, stringsAsFactors = FALSE)
    if (!("sample" %in% names(Sdf))) {
      rn <- rownames(Sdf)
      if (is.null(rn)) stop("Scores for ", m, " have no rownames and no 'sample' column.")
      Sdf$sample <- rn
    }
    Sdf$sample <- as.character(Sdf$sample)
    
    # left-join by sample
    merged <- left_join(Sdf, meta, by = setNames(sample_col, "sample"))
    merged$method <- m
    
    if (isTRUE(replace_scores)) {
      # keep original numeric matrix under a backup name, replace scores with merged df
      bench[[m]]$scores_matrix <- bench[[m]]$scores
      bench[[m]]$scores <- merged
    } else {
      # store alongside the original scores matrix
      bench[[m]][[store_name]] <- merged
    }
  }
  bench
}

#merged <- attach_scores_metadata(benchmark_cll, CLL_metadata, sample_col = "sample")

#============================================
# Plot based on metadata
#============================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(grid)   # unit()

# ---------------- utilities ----------------
resolve_factor_col <- function(df, method, factor_idx = 1) {
  cand  <- paste0("F", factor_idx, " (", method, ")")
  short <- paste0("F", factor_idx)
  if (cand  %in% names(df)) return(cand)
  if (short %in% names(df)) return(short)
  hit <- grep(paste0("^F", factor_idx, "(\\b|\\s*\\()"), names(df), value = TRUE)
  if (length(hit)) hit[1] else stop("Could not find F", factor_idx, " for ", method)
}

# robust palette builder: user overrides defaults; input can be NULL, vector, or list
build_palette <- function(pal = NULL) {
  default <- c("0"="#E37449", "1"="#00366C", "NA"="#999999")#
  if (is.null(pal)) return(default)
  if (is.list(pal)) pal <- unlist(pal, use.names = TRUE)
  pal <- as.character(pal)
  if (is.null(names(pal))) stop("Custom palette must be a *named* vector with names '0','1','NA'.")
  default[names(pal)] <- pal
  default
}

# dataframe for one method
scores_df_for_method <- function(bench, method, factor_idx = 1, group_col = "IGHV") {
  x <- bench[[method]][["scores_with_meta"]]
  if (is.null(x)) stop("scores_with_meta missing for ", method, ". Run attach_scores_metadata() first.")
  score_col <- resolve_factor_col(x, method, factor_idx)
  
  grp_raw <- if (group_col %in% names(x)) x[[group_col]] else NA
  grp_chr <- dplyr::case_when(
    is.na(grp_raw) ~ "NA",
    as.character(grp_raw) %in% c("0","1") ~ as.character(grp_raw),
    TRUE ~ "NA"
  )
  
  tibble(
    sample      = if ("sample" %in% names(x)) x$sample else seq_len(nrow(x)),
    .sample_idx = seq_len(nrow(x)),
    .score      = as.numeric(x[[score_col]]),
    .group      = grp_chr,
    method      = method
  )
}

# ---------------- panel ----------------
plot_scores_panel <- function(
    df, y_lab, panel_tag = NULL,
    tag_position = c("tl_in","none"),
    palette = NULL
) {
  tag_position <- match.arg(tag_position)
  pal <- build_palette(palette)
  df$.group <- factor(df$.group, levels = c("0","1","NA"))
  
  p <- ggplot(df, aes(x = .sample_idx, y = .score, fill = .group, color = .group)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(shape = 21, size = 2.7, stroke = 0.8, alpha = 0.95) +
    # ⬇️ one legend only (from fill); hide color legend
    scale_fill_manual(values = pal, drop = FALSE, name = ".group") +
    scale_color_manual(values = pal, drop = FALSE, guide = "none") +
    labs(x = "Samples", y = y_lab) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title.x = element_text(size = 14, face = "plain", colour = "black"),#1f3765 #0072B2
      axis.title.y = element_text(size = 14, face = "plain", colour = "black"),
      plot.margin = margin(6, 6, 6, 6)
    )
  
  if (tag_position == "tl_in" && !is.null(panel_tag) && nzchar(panel_tag)) {
    p <- p + annotate("label", x = -Inf, y = Inf, label = panel_tag,
                      hjust = -0.1, vjust = 1.1, size = 4.5,
                      fill = "grey90", label.r = unit(0.15, "lines"))
  }
  p
}

# small plot used for OUTSIDE tag row
make_tag_plot <- function(tag) {
  if (is.null(tag) || !nzchar(tag)) return(patchwork::plot_spacer())
  ggplot() +
    annotate("label", x = 0, y = 0, label = tag,
             size = 4.5, fill = "grey90", label.r = unit(0.15, "lines")) +
    xlim(-1, 1) + ylim(-1, 1) + theme_void() +
    theme(plot.margin = margin(0, 6, 0, 6))
}
# # tag_position: "tl_in" (inside, top-left), "none" (no tag)
# plot_scores_panel <- function(
    #     df, y_lab, panel_tag = NULL,
#     tag_position = c("tl_in","none"),
#     palette = NULL
# ) {
#   tag_position <- match.arg(tag_position)
#   pal <- build_palette(palette)
#   df$.group <- factor(df$.group, levels = c("0","1","NA"))
#   
#   p <- ggplot(df, aes(x = .sample_idx, y = .score, fill = .group, color = .group)) +
#     geom_hline(yintercept = 0, linetype = "dashed") +
#     geom_point(shape = 21, size = 2.7, stroke = 0.8, alpha = 0.95) +
#     scale_fill_manual(values = pal, drop = FALSE) +
#     scale_color_manual(values = pal, drop = FALSE) +
#     labs(x = "Samples", y = y_lab) +
#     coord_cartesian(clip = "off") +
#     theme_minimal(base_size = 13) +
#     theme(
#       panel.grid.minor = element_blank(),
#       panel.grid.major.x = element_blank(),
#       axis.title.x = element_text(size = 14, face = "bold", colour = "#1f3765"),
#       axis.title.y = element_text(size = 14, face = "bold", colour = "#1f3765"),
#       plot.margin = margin(6, 6, 6, 6)
#     )
#   
#   if (tag_position == "tl_in" && !is.null(panel_tag) && nzchar(panel_tag)) {
#     p <- p + annotate(
#       "label", x = -Inf, y = Inf, label = panel_tag,
#       hjust = -0.1, vjust = 1.1, size = 4.5,
#       fill = "grey90", label.r = unit(0.15, "lines")
#     )
#   }
#   p
# }
# 
# # small plot used for OUTSIDE tag row
# make_tag_plot <- function(tag) {
#   if (is.null(tag) || !nzchar(tag)) return(patchwork::plot_spacer())
#   ggplot() +
#     annotate("label", x = 0, y = 0, label = tag,
#              size = 4.5, fill = "grey90", label.r = unit(0.15, "lines")) +
#     xlim(-1, 1) + ylim(-1, 1) + theme_void() +
#     theme(plot.margin = margin(0, 6, 0, 6))
# }

plot_scores_scatter_grid <- function(
    bench,
    factor_idx      = 1,
    group_col       = "IGHV",
    methods         = c("MOFA","FABIA","MFA","GFA"),
    panel_tags      = NULL,                 # e.g., c("A","B","C","D") or NULL
    tag_position    = c("tl_in","tl_out","none"),
    title_text      = "Factor scores (scatter plots)",
    palette         = NULL,                 # passed to plot_scores_panel()
    legend_position = "right",              # "right" or "bottom"
    title_height    = 0.15,                 # relative space for title row
    tag_row_height  = 0.10                  # relative space for outside tag row
) {
  stopifnot(length(methods) == 4)
  tag_position <- match.arg(tag_position)
  
  # local tag-plot helper (only used when tag_position == "tl_out")
  .make_tag_plot <- function(tag) {
    if (is.null(tag) || !nzchar(tag)) return(patchwork::plot_spacer())
    ggplot() +
      annotate("label", x = 0, y = 0, label = tag,
               size = 4.5, fill = "grey90", label.r = grid::unit(0.15, "lines")) +
      xlim(-1, 1) + ylim(-1, 1) + theme_void() +
      theme(plot.margin = margin(0, 6, 0, 6))
  }
  
  # data for the four panels (relies on your scores_df_for_method)
  dfs <- lapply(methods, function(m) scores_df_for_method(bench, m, factor_idx, group_col))
  names(dfs) <- methods
  
  # normalize panel tags
  if (is.null(panel_tags)) tags <- rep("", 4) else tags <- rep_len(panel_tags, 4)
  all_empty <- all(!nzchar(tags))
  
  # build panels (relies on your plot_scores_panel)
  pA <- plot_scores_panel(dfs[[methods[1]]], paste0("Score (", methods[1], ")"),
                          panel_tag = if (tag_position == "tl_in") tags[1] else NULL,
                          tag_position = if (tag_position == "tl_in") "tl_in" else "none",
                          palette = palette)
  pB <- plot_scores_panel(dfs[[methods[2]]], paste0("Score (", methods[2], ")"),
                          panel_tag = if (tag_position == "tl_in") tags[2] else NULL,
                          tag_position = if (tag_position == "tl_in") "tl_in" else "none",
                          palette = palette)
  pC <- plot_scores_panel(dfs[[methods[3]]], paste0("Score (", methods[3], ")"),
                          panel_tag = if (tag_position == "tl_in") tags[3] else NULL,
                          tag_position = if (tag_position == "tl_in") "tl_in" else "none",
                          palette = palette)
  pD <- plot_scores_panel(dfs[[methods[4]]], paste0("Score (", methods[4], ")"),
                          panel_tag = if (tag_position == "tl_in") tags[4] else NULL,
                          tag_position = if (tag_position == "tl_in") "tl_in" else "none",
                          palette = palette)
  
  # put the legend settings onto each panel (avoids the '&' operator)
  leg_theme <- theme(legend.position = legend_position,
                     legend.title = element_text(face = "bold"))
  pA <- pA + leg_theme; pB <- pB + leg_theme; pC <- pC + leg_theme; pD <- pD + leg_theme
  
  # title
  title_grob <- ggplot() +
    annotate("text", x = 0, y = 0, label = title_text,
             fontface = 2, size = 5.5, colour = "#6a1b9a") +
    theme_void() + theme(plot.margin = margin(0,0,6,0))
  
  # 2×2 grid of panels
  panel_grid <- (pA | pB) / (pC | pD)
  
  # assemble with optional outside tag row and collect guides → one legend
  if (tag_position == "tl_out" && !all_empty) {
    tag_row <- .make_tag_plot(tags[1]) | .make_tag_plot(tags[2]) |
      .make_tag_plot(tags[3]) | .make_tag_plot(tags[4])
    title_grob / tag_row / panel_grid +
      plot_layout(heights = c(title_height, tag_row_height, 1), guides = "collect")
  } else {
    title_grob / panel_grid +
      plot_layout(heights = c(title_height, 1), guides = "collect")
  }
}

# p_none <- plot_scores_scatter_grid(
#   benchmark_cll, factor_idx = 1, group_col = "IGHV",
#   methods = c("MOFA","FABIA","MFA","GFA"),
#   panel_tags = NULL,
#   tag_position = "none"
# )
# print(p_none)

#============================================
# Create scores dataframes
#============================================

# Scores
fabia_score_df <- data.frame(sim_benchmark_res$FABIA$scores); fabia_score_df$sample <- rownames(fabia_score_df)
# rename columns that start with "F"
colnames(fabia_score_df) <- sub("^F(\\d+)$", "F\\1 (FABIA)", colnames(fabia_score_df))

mofa_score_df <- data.frame(sim_benchmark_res$MOFA$scores); mofa_score_df$sample <- rownames(mofa_score_df)
# rename columns that start with "F"
colnames(mofa_score_df) <- sub("^F(\\d+)$", "F\\1 (MOFA)", colnames(mofa_score_df))

mfa_score_df <- data.frame(sim_benchmark_res$MFA$scores); mfa_score_df$sample <- rownames(mfa_score_df)
# rename columns that start with "F"
colnames(mfa_score_df) <- sub("^F(\\d+)$", "F\\1 (MFA)", colnames(mfa_score_df))

gfa_score_df <- data.frame(sim_benchmark_res$GFA$scores); gfa_score_df$sample <- rownames(gfa_score_df)
# rename columns that start with "F"
colnames(gfa_score_df) <- sub("^F(\\d+)$", "F\\1 (GFA)", colnames(gfa_score_df))


# Apply to each method
swap_F1_F2_names <- function(df) {
  nms <- names(df)
  
  f1_idx <- grep("F1", nms, fixed = TRUE)
  f2_idx <- grep("F2", nms, fixed = TRUE)
  
  if (length(f1_idx) != 1 || length(f2_idx) != 1) {
    stop("Could not uniquely detect columns with 'F1' and 'F2' in their names.")
  }
  
  # swap just the names
  tmp <- nms[f1_idx]
  nms[f1_idx] <- nms[f2_idx]
  nms[f2_idx] <- tmp
  
  names(df) <- nms
  df
}
fabia_score_df <- swap_F1_F2_names(fabia_score_df)
gfa_score_df <- swap_F1_F2_names(gfa_score_df)
names(fabia_score_df)
names(mofa_score_df)
names(mfa_score_df)
names(gfa_score_df)
names(simGTscores)
colnames(simGTscores) <- sub("^F(\\d+)$", "F\\1 (GT)", colnames(simGTscores))

library(dplyr)
library(ggplot2)
library(reshape2)

# SCORES VISUALIZATION
names(fabia_score_df)
names(mofa_score_df)
names(mfa_score_df)
names(gfa_score_df)

#### ============================================================
#### 1) Ground truth extraction (scores + per-omic loadings)
####    - Scores rows:  sample_1, sample_2, ...
####    - Loadings rows: "rna::omic1_feature_i", "prot::omic2_feature_i", ...
####    - Per-omic data.frames also carry a `feature` column for merges.
#### ============================================================

# ============================================================
# Robust benchmark for factor methods with differing #columns
# ============================================================
# ---------- Pairs-plot panels (base graphics) -------------  # NEW
.panel_hist <- function(x) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  h <- hist(x, plot = FALSE)
  y <- if (max(h$counts) == 0) h$counts else h$counts / max(h$counts)
  rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], y, col = "grey85", border = "white")
  dx <- try(density(x, na.rm = TRUE), silent = TRUE)
  if (!inherits(dx, "try-error")) lines(dx$x, dx$y / max(dx$y), lwd = 1)
}

.panel_cor <- function(x, y) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r  <- suppressWarnings(stats::cor(x, y, use = "pairwise.complete.obs"))
  txt <- if (is.finite(r)) sprintf("%.2f", abs(r)) else "NA"
  text(0.5, 0.5, txt, cex = 2.4, font = 2)
}

.panel_scatter <- function(x, y) {
  points(x, y, pch = 8, cex = 0.5)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) >= 3) abline(lm(y ~ x), col = "red", lwd = 1)
}

.make_pairs_plot <- function(mat, main = "", oma_top = 2) {
  old <- par(no.readonly = TRUE); on.exit(par(old))
  par(oma = c(0, 0, oma_top, 0))
  pairs(
    mat,
    diag.panel = .panel_hist,
    upper.panel = .panel_cor,
    lower.panel = .panel_scatter
  )
  if (nzchar(main)) mtext(main, outer = TRUE, line = 0.5, cex = 1.2, font = 2)
  recordPlot()
}

# Flip 'v' so that cor(ref, v) >= 0 (pairwise NA-safe). If cor is NA, do nothing.
.align_sign_vec <- function(ref, v) {
  r <- suppressWarnings(stats::cor(ref, v, use = "pairwise.complete.obs"))
  if (is.finite(r) && r < 0) -v else v
}

# ---------- Reference-based factor alignment ----------------  # NEW
# Map each method's factors to a reference method using Hungarian on |cor|.
.align_to_reference_scores <- function(scores_list, ref, match_factors_fun) {
  ref_mat <- scores_list[[ref]]
  ref_k   <- ncol(ref_mat)
  maps <- list()
  for (m in names(scores_list)) if (m != ref) {
    common <- intersect(rownames(ref_mat), rownames(scores_list[[m]]))
    A <- ref_mat[common, , drop = FALSE]
    B <- scores_list[[m]][common, , drop = FALSE]
    C <- stats::cor(A, B, use = "pairwise.complete.obs")
    ms <- match_factors_fun(C)
    map <- rep(NA_integer_, ref_k)
    if (!is.null(ms) && nrow(ms) > 0) map[ms$A] <- ms$B
    maps[[m]] <- map
  }
  maps
}

.align_to_reference_loadings <- function(loadings_list, ref, get_prefix, get_suffix, match_factors_fun) {
  ref_mat <- loadings_list[[ref]]
  ref_k   <- ncol(ref_mat)
  pref_ref <- get_prefix(rownames(ref_mat))
  ds <- unique(pref_ref)
  out <- setNames(vector("list", length(ds)), ds)
  for (om in ds) {
    idx_ref <- which(pref_ref == om)
    ids_ref <- get_suffix(rownames(ref_mat)[idx_ref])
    maps_om <- list()
    for (m in names(loadings_list)) if (m != ref) {
      pref_m <- get_prefix(rownames(loadings_list[[m]]))
      idx_m  <- which(pref_m == om)
      if (!length(idx_m)) { maps_om[[m]] <- rep(NA_integer_, ref_k); next }
      ids_m  <- get_suffix(rownames(loadings_list[[m]])[idx_m])
      common <- intersect(ids_ref, ids_m)
      if (!length(common)) { maps_om[[m]] <- rep(NA_integer_, ref_k); next }
      A <- ref_mat[idx_ref[match(common, ids_ref)], , drop = FALSE]
      B <- loadings_list[[m]][idx_m[match(common, ids_m)], , drop = FALSE]
      C <- stats::cor(A, B, use = "pairwise.complete.obs")
      ms <- match_factors_fun(C)
      map <- rep(NA_integer_, ref_k)
      if (!is.null(ms) && nrow(ms) > 0) map[ms$A] <- ms$B
      maps_om[[m]] <- map
    }
    out[[om]] <- maps_om
  }
  out
}

benchmark_factor_methods <- function(
    benchmark_res,
    methods = c("FABIA","MOFA","GFA","MFA","GT"),
    rename_cols = TRUE,          # label factor cols "F1 (METHOD)" etc.
    center_scale = TRUE,         # z-score each factor column
    feature_id_split = "::",     # used to align features (suffix after this)
    ground_truth_feature = NULL, # optional named 0/1 vector
    ground_truth_sample  = NULL, # optional named 0/1 vector
    plot_theme_base_size = 11){
  # ---- deps ----
  req <- c("ggplot2","reshape2")
  miss <- req[!vapply(req, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(miss)) stop("Please install required packages: ", paste(miss, collapse=", "))
  if (!requireNamespace("clue", quietly = TRUE)) stop("Please install package 'clue'")
  
  `%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
  
  # ---------- robust numeric ingestion ----------
  coerce_numeric_cols <- function(df) {
    df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
    f_like <- grepl("^F[0-9]+\\s*\\(", names(df))
    for (j in which(f_like)) if (!is.numeric(df[[j]])) df[[j]] <- suppressWarnings(as.numeric(df[[j]]))
    non_f <- which(!f_like)
    for (j in non_f) if (!is.numeric(df[[j]])) {
      tc <- suppressWarnings(type.convert(df[[j]], as.is = TRUE))
      if (is.numeric(tc)) df[[j]] <- tc
    }
    df
  }
  
  to_numeric_matrix <- function(obj, method_name, id_pref = c("sample","feature")) {
    df <- coerce_numeric_cols(obj)
    rn <- rownames(df)
    if (is.null(rn) || anyNA(rn) || any(rn == "")) {
      hit <- intersect(id_pref, names(df))
      if (length(hit)) rn <- as.character(df[[hit[1]]])
    }
    f_like <- grepl("^F[0-9]+\\s*\\(", names(df))
    keep <- f_like & vapply(df, is.numeric, logical(1))
    if (!any(keep)) {
      drop_names <- c("sample","feature", grep("^signal_", names(df), value = TRUE))
      keep <- vapply(df, is.numeric, logical(1)) & !names(df) %in% drop_names
    }
    if (!any(keep)) {
      stop("No numeric columns found for method ", method_name,
           ". Columns were: ", paste(names(df), collapse = ", "))
    }
    mat <- as.matrix(df[, keep, drop = FALSE])
    storage.mode(mat) <- "double"
    if (is.null(rn)) {
      rn <- if (identical(id_pref[1], "sample")) paste0("sample_", seq_len(nrow(mat)))
      else paste0("feature_", seq_len(nrow(mat)))
    }
    rownames(mat) <- rn
    mat
  }
  
  # robust scaling that never drops dimensions
  safe_scale <- function(X) {
    if (!center_scale) return(as.matrix(X))
    X <- as.matrix(X)
    if (ncol(X) == 0L) return(X)
    ok <- colSums(!is.na(X)) > 0
    if (!any(ok)) return(X)
    tmp <- scale(X[, ok, drop = FALSE])
    if (is.null(dim(tmp))) {
      tmp <- matrix(tmp, nrow = nrow(X), ncol = sum(ok),
                    dimnames = list(rownames(X), colnames(X)[ok]))
    }
    X[, ok] <- as.matrix(tmp)
    X
  }
  
  rename_factor_cols <- function(X, method_name) {
    k <- ncol(X)
    if (rename_cols) {
      colnames(X) <- paste0("F", seq_len(k), " (", method_name, ")")
    } else if (is.null(colnames(X))) {
      colnames(X) <- paste0("F", seq_len(k))
    }
    X
  }
  
  # ---- 0) Which methods are present ----
  available <- intersect(methods, names(benchmark_res))
  if (length(available) < 2)
    stop("Need at least two methods in 'benchmark_res'. Found: ", paste(available, collapse=", "))
  
  # ---- 1) Collect & standardize matrices ----
  get_scores <- function(m) {
    X <- to_numeric_matrix(benchmark_res[[m]]$scores, m, id_pref = "sample")
    X <- rename_factor_cols(X, m)
    X <- safe_scale(X)
    X
  }
  get_loadings <- function(m) {
    X <- to_numeric_matrix(benchmark_res[[m]]$loadings, m, id_pref = "feature")
    X <- rename_factor_cols(X, m)
    X <- safe_scale(X)
    X
  }
  
  scores_list   <- lapply(available, get_scores);   names(scores_list)   <- available
  loadings_list <- lapply(available, get_loadings); names(loadings_list) <- available
  
  # ---- 2) Helpers for per-dataset loadings --------------------------------
  get_prefix <- function(rn) sub("::.*$", "", rn)       # "rna", "prot", "omic3", ...
  get_suffix <- function(rn) sub("^.*::", "", rn)       # "omic1_feature_123", ...
  normalize_ids <- function(xnames) {
    if (is.null(xnames)) return(NULL)
    if (is.null(feature_id_split)) return(xnames)
    parts <- strsplit(xnames, split = feature_id_split, fixed = TRUE)
    vapply(parts, function(p) if (length(p)) tail(p, 1) else "", character(1))
  }
  
  # Per-dataset loading correlations
  cor_mat_loadings_peromic <- function(A, B) {
    pA <- get_prefix(rownames(A)); pB <- get_prefix(rownames(B))
    omics <- intersect(unique(pA), unique(pB))
    out <- vector("list", length(omics)); names(out) <- omics
    for (om in omics) {
      idxA <- which(pA == om); idxB <- which(pB == om)
      if (!length(idxA) || !length(idxB)) next
      idsA <- get_suffix(rownames(A)[idxA])
      idsB <- get_suffix(rownames(B)[idxB])
      common <- intersect(idsA, idsB)
      if (!length(common)) next
      AA <- A[idxA[match(common, idsA)], , drop = FALSE]
      BB <- B[idxB[match(common, idsB)], , drop = FALSE]
      out[[om]] <- stats::cor(AA, BB, use = "pairwise.complete.obs")
    }
    Filter(Negate(is.null), out)
  }
  
  # --- aligned vectors for scattering loadings within a dataset ---
  get_aligned_loading_vectors <- function(A, B, dataset, i, j) {
    pA <- get_prefix(rownames(A)); pB <- get_prefix(rownames(B))
    idxA <- which(pA == dataset); idxB <- which(pB == dataset)
    if (!length(idxA) || !length(idxB)) return(NULL)
    idsA <- get_suffix(rownames(A)[idxA]); idsB <- get_suffix(rownames(B)[idxB])
    common <- intersect(idsA, idsB)
    if (!length(common)) return(NULL)
    x <- A[idxA[match(common, idsA)], i]
    y <- B[idxB[match(common, idsB)], j]
    list(x = x, y = y)
  }
  
  # ---- plotting helpers (scatter with r) -----------------------------------
  make_scatter <- function(x, y, title, xlab, ylab) {
    df <- data.frame(x = as.numeric(x), y = as.numeric(y))
    r  <- suppressWarnings(stats::cor(df$x, df$y, use = "pairwise.complete.obs"))
    ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(alpha = 0.6, size = 1) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.4, linetype = "dashed") +
      ggplot2::annotate("text", x = Inf, y = Inf, label = sprintf("r = %.3f", r),
                        hjust = 1.05, vjust = 1.5, size = 4) +
      ggplot2::labs(title = title, x = xlab, y = ylab) +
      ggplot2::theme_minimal(base_size = plot_theme_base_size) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  }
  
  # ---- 3) Factor matching (Hungarian on |corr|), NA-safe -------------------
  match_factors <- function(C) {
    if (is.null(C) || any(dim(C) == 0)) return(NULL)
    C[!is.finite(C)] <- NA_real_
    absC <- abs(C); absC[is.na(absC)] <- 0
    ka <- nrow(absC); kb <- ncol(absC)
    maxv <- if (all(absC == 0)) 0 else max(absC, na.rm = TRUE)
    cost <- maxv - absC
    pad_val <- if (is.finite(suppressWarnings(max(cost, na.rm = TRUE))))
      suppressWarnings(max(cost, na.rm = TRUE)) + 1 else 1
    if (ka > kb) cost <- cbind(cost, matrix(pad_val, nrow = ka, ncol = ka - kb))
    if (kb > ka) cost <- rbind(cost, matrix(pad_val, nrow = kb - ka, ncol = kb))
    storage.mode(cost) <- "double"
    ass <- clue::solve_LSAP(cost)
    real_rows <- seq_len(ka)
    j_all     <- as.integer(ass[real_rows])
    keep      <- which(j_all <= kb)
    if (!length(keep)) return(NULL)
    i_keep <- real_rows[keep]; j_keep <- j_all[keep]
    data.frame(A = i_keep, B = j_keep, corr = mapply(function(i,j) C[i,j], i_keep, j_keep),
               stringsAsFactors = FALSE)
  }
  
  # ---- 4) Pairwise similarities (+ scatters) -------------------------------
  method_pairs <- t(combn(available, 2))
  pair_summaries_scores   <- list()
  pair_summaries_loadings <- list()
  pair_heatmaps_scores    <- list()
  pair_heatmaps_loadings  <- list()
  pair_scatter_scores     <- list()   # NEW
  pair_scatter_loadings   <- list()   # NEW
  
  for (r in seq_len(nrow(method_pairs))) {
    m1 <- method_pairs[r,1]; m2 <- method_pairs[r,2]
    
    # ----- SCORES -----
    common_samp <- intersect(rownames(scores_list[[m1]]), rownames(scores_list[[m2]]))
    A_s <- scores_list[[m1]][common_samp, , drop = FALSE]
    B_s <- scores_list[[m2]][common_samp, , drop = FALSE]
    Cs  <- stats::cor(A_s, B_s, use = "pairwise.complete.obs")
    ms  <- match_factors(Cs)
    
    # summary + heatmap
    pair_summaries_scores[[paste(m1,m2,sep="|")]] <- list(
      summary = data.frame(
        pair = paste(m1, m2, sep = " vs "),
        k_m1_scores = ncol(A_s),
        k_m2_scores = ncol(B_s),
        k_matched_scores = if (!is.null(ms)) nrow(ms) else 0L,
        mean_abs_corr_scores = if (!is.null(ms)) mean(abs(ms$corr), na.rm = TRUE) else NA_real_,
        median_abs_corr_scores = if (!is.null(ms)) stats::median(abs(ms$corr), na.rm = TRUE) else NA_real_,
        max_abs_corr_scores = if (!is.null(ms)) max(abs(ms$corr), na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
      ),
      scores_cor = Cs, scores_match = ms
    )
    
    dfh <- reshape2::melt(Cs,
                          varnames = c(paste0("Factor (", m1,")"),
                                       paste0("Factor (", m2,")")),
                          value.name = "Correlation")
    xcol <- names(dfh)[2]; ycol <- names(dfh)[1]
    pair_heatmaps_scores[[paste(m1,m2,sep="|")]] <-
      ggplot2::ggplot(dfh, ggplot2::aes(x = .data[[xcol]], y = .data[[ycol]], fill = .data[["Correlation"]])) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_gradient2(limits = c(-1,1), midpoint = 0, name = "Corr") +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = paste0("Scores correlation: ", m1, " vs ", m2), x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = plot_theme_base_size) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(face = "bold"))
    
    # score scatters for matched pairs
    ss_plots <- list()
    if (!is.null(ms) && nrow(ms) > 0) {
      for (k in seq_len(nrow(ms))) {
        i <- ms$A[k]; j <- ms$B[k]
        x <- A_s[, i]; y <- B_s[, j]
        title <- sprintf("Scores scatter: %s[%s] vs %s[%s]",
                         m1, colnames(A_s)[i], m2, colnames(B_s)[j])
        ss_plots[[k]] <- make_scatter(x, y, title,
                                      xlab = paste0(colnames(A_s)[i]),
                                      ylab = paste0(colnames(B_s)[j]))
      }
    }
    pair_scatter_scores[[paste(m1,m2,sep="|")]] <- ss_plots
    
    # ----- LOADINGS (per dataset) -----
    Cl_list <- cor_mat_loadings_peromic(loadings_list[[m1]], loadings_list[[m2]])
    ml_list <- lapply(Cl_list, match_factors)
    
    # per-dataset summary + heatmap + scatters
    rows <- list()
    scat_by_ds <- list()
    for (om in names(Cl_list)) {
      ml <- ml_list[[om]]
      rows[[length(rows)+1]] <- data.frame(
        pair = paste(m1, m2, sep = " vs "),
        dataset = om,
        k_m1_loadings = ncol(loadings_list[[m1]]),
        k_m2_loadings = ncol(loadings_list[[m2]]),
        k_matched_loadings = if (!is.null(ml)) nrow(ml) else 0L,
        mean_abs_corr_loadings = if (!is.null(ml)) mean(abs(ml$corr), na.rm = TRUE) else NA_real_,
        median_abs_corr_loadings = if (!is.null(ml)) stats::median(abs(ml$corr), na.rm = TRUE) else NA_real_,
        max_abs_corr_loadings = if (!is.null(ml)) max(abs(ml$corr), na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
      )
      
      # heatmap for this dataset
      dfhL <- reshape2::melt(Cl_list[[om]],
                             varnames = c(paste0("Factor (", m1,")"),
                                          paste0("Factor (", m2,")")),
                             value.name = "Correlation")
      xL <- names(dfhL)[2]; yL <- names(dfhL)[1]
      pair_heatmaps_loadings[[paste(m1,m2,om,sep="|")]] <-
        ggplot2::ggplot(dfhL, ggplot2::aes(x = .data[[xL]], y = .data[[yL]], fill = .data[["Correlation"]])) +
        ggplot2::geom_tile(color = "white", linewidth = 0.2) +
        ggplot2::scale_fill_gradient2(limits = c(-1,1), midpoint = 0, name = "Corr") +
        ggplot2::coord_fixed() +
        ggplot2::labs(title = paste0("Loadings correlation (", om, "): ", m1, " vs ", m2),
                      x = NULL, y = NULL) +
        ggplot2::theme_minimal(base_size = plot_theme_base_size) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                       panel.grid = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(face = "bold"))
      
      # loading scatters for matched pairs within this dataset
      ds_plots <- list()
      if (!is.null(ml) && nrow(ml) > 0) {
        for (k in seq_len(nrow(ml))) {
          i <- ml$A[k]; j <- ml$B[k]
          vecs <- get_aligned_loading_vectors(loadings_list[[m1]], loadings_list[[m2]], om, i, j)
          if (is.null(vecs)) next
          title <- sprintf("Loadings scatter (%s): %s[%s] vs %s[%s]",
                           om, m1, colnames(loadings_list[[m1]])[i], m2, colnames(loadings_list[[m2]])[j])
          ds_plots[[k]] <- make_scatter(vecs$x, vecs$y, title,
                                        xlab = paste0(colnames(loadings_list[[m1]])[i]),
                                        ylab = paste0(colnames(loadings_list[[m2]])[j]))
        }
      }
      scat_by_ds[[om]] <- ds_plots
    }
    
    pair_summaries_loadings[[paste(m1,m2,sep="|")]] <- list(
      cor   = Cl_list,
      match = ml_list,
      summary = if (length(rows)) do.call(rbind, rows) else NULL
    )
    pair_scatter_loadings[[paste(m1,m2,sep="|")]] <- scat_by_ds
  }
  
  # ======================== 5b) Pairs-style plots (NEW) =========================
  # Reference method for alignment
  ref_method <- if ("FABIA" %in% available) "FABIA" else available[1]
  
  # ---- Scores: build a pairs plot per reference factor -------------------------
  maps_scores <- .align_to_reference_scores(scores_list, ref_method, match_factors)
  
  # samples common to all selected methods
  common_samples_all <- Reduce(intersect, lapply(scores_list, rownames))
  pairs_scores_by_factor <- list()
  ref_k <- ncol(scores_list[[ref_method]])
  
  # for (i in seq_len(ref_k)) {
  #   cols <- list()
  #   for (m in available) {
  #     if (m == ref_method) {
  #       cols[[paste0("Score (", m, ")")]] <- scores_list[[m]][common_samples_all, i]
  #     } else {
  #       j <- maps_scores[[m]][i]
  #       if (is.na(j)) next
  #       cols[[paste0("Score (", m, ")")]] <- scores_list[[m]][common_samples_all, j]
  #     }
  #   }
  #   if (length(cols) >= 2) {
  #     mat <- as.data.frame(cols, check.names = FALSE)
  #     pairs_scores_by_factor[[paste0("F", i)]] <-
  #       .make_pairs_plot(mat, main = sprintf("Factor scores (cor.) — aligned to %s F%d", ref_method, i))
  #   }
  # }
  for (i in seq_len(ref_k)) {
    cols <- list()
    ref_vec <- scores_list[[ref_method]][common_samples_all, i]
    
    for (m in available) {
      if (m == ref_method) {
        cols[[paste0("Score (", m, ")")]] <- ref_vec
      } else {
        j <- maps_scores[[m]][i]
        if (is.na(j)) next
        v <- scores_list[[m]][common_samples_all, j]
        v <- .align_sign_vec(ref_vec, v)        # << align sign to reference
        cols[[paste0("Score (", m, ")")]] <- v
      }
    }
    if (length(cols) >= 2) {
      mat <- as.data.frame(cols, check.names = FALSE)
      pairs_scores_by_factor[[paste0("F", i)]] <-
        .make_pairs_plot(mat, main = sprintf("Factor scores (Pearson cor.)"))
    }
  }
  
  # ---- Per-omic loadings: pairs plot per dataset & factor ----------------------
  maps_loadings <- .align_to_reference_loadings(loadings_list, ref_method, get_prefix, get_suffix, match_factors)
  pairs_loadings_by_dataset_and_factor <- list()
  
  # datasets seen in reference loadings
  datasets_ref <- unique(sub("::.*$", "", rownames(loadings_list[[ref_method]])))
  
  for (om in datasets_ref) {
    pref_ref <- get_prefix(rownames(loadings_list[[ref_method]]))
    idx_ref  <- which(pref_ref == om)
    ids_ref  <- get_suffix(rownames(loadings_list[[ref_method]])[idx_ref])
    
    # universal common features across all methods for this dataset
    common_ids <- ids_ref
    for (m in available) if (m != ref_method) {
      pref_m <- get_prefix(rownames(loadings_list[[m]]))
      ids_m  <- get_suffix(rownames(loadings_list[[m]])[pref_m == om])
      common_ids <- intersect(common_ids, ids_m)
    }
    if (!length(common_ids)) next
    
    # row indices by method for the same features in the same order
    row_idx <- list()
    row_idx[[ref_method]] <- idx_ref[match(common_ids, ids_ref)]
    for (m in available) if (m != ref_method) {
      pref_m <- get_prefix(rownames(loadings_list[[m]]))
      idx_m  <- which(pref_m == om)
      ids_m  <- get_suffix(rownames(loadings_list[[m]])[idx_m])
      row_idx[[m]] <- idx_m[match(common_ids, ids_m)]
    }
    
    # per-factor pairs
    ref_k <- ncol(loadings_list[[ref_method]])
    #   for (i in seq_len(ref_k)) {
    #     cols <- list()
    #     for (m in available) {
    #       if (m == ref_method) {
    #         cols[[paste0(m)]] <- loadings_list[[m]][row_idx[[m]], i]
    #       } else {
    #         j <- maps_loadings[[om]][[m]][i]
    #         if (is.na(j)) next
    #         cols[[paste0(m)]] <- loadings_list[[m]][row_idx[[m]], j]
    #       }
    #     }
    #     if (length(cols) >= 2) {
    #       mat <- as.data.frame(cols, check.names = FALSE)
    #       lab <- paste0("Loadings (", om, ") cor. — aligned to ", ref_method, " F", i)
    #       pairs_loadings_by_dataset_and_factor[[om]][[paste0("F", i)]] <- .make_pairs_plot(mat, main = lab)
    #     }
    #   }
    
    for (i in seq_len(ref_k)) {
      cols <- list()
      ref_vec <- loadings_list[[ref_method]][row_idx[[ref_method]], i]  # same features/order
      
      for (m in available) {
        if (m == ref_method) {
          cols[[paste0(m)]] <- ref_vec
        } else {
          j <- maps_loadings[[om]][[m]][i]
          if (is.na(j)) next
          v <- loadings_list[[m]][row_idx[[m]], j]
          v <- .align_sign_vec(ref_vec, v)        # << align sign to reference
          cols[[paste0(m)]] <- v
        }
      }
      if (length(cols) >= 2) {
        mat <- as.data.frame(cols, check.names = FALSE)
        lab <- paste0("Loadings (", om, ") — Pearson cor. ")
        pairs_loadings_by_dataset_and_factor[[om]][[paste0("F", i)]] <- .make_pairs_plot(mat, main = lab)
      }
    }
  }
  
  # ---- 6) Global all-in-one correlation heatmap (scores) ----
  merge_scores <- function(lst){
    X <- lst[[1]]; colnames(X) <- make.unique(colnames(X))
    for (i in seq(2, length(lst))) {
      Y <- lst[[i]]; colnames(Y) <- make.unique(colnames(Y))
      common <- intersect(rownames(X), rownames(Y))
      X <- cbind(X[common,,drop=FALSE], Y[common,,drop=FALSE])
    }
    X
  }
  all_scores_mat <- merge_scores(scores_list)
  all_cor <- stats::cor(all_scores_mat, use = "pairwise.complete.obs")
  
  get_method <- function(nm) sub(".*\\(([^)]+)\\)\\s*$", "\\1", nm)
  get_fnum   <- function(nm) suppressWarnings(as.integer(sub("^\\s*F(\\d+).*", "\\1", nm)))
  mth <- vapply(colnames(all_cor), get_method, character(1))
  fno <- vapply(colnames(all_cor), get_fnum, numeric(1))
  lev <- unique(mth)
  ord <- order(match(mth, lev), fno, colnames(all_cor))
  all_cor <- all_cor[ord, ord, drop = FALSE]
  
  dfall <- reshape2::melt(all_cor, varnames = c("Factor1","Factor2"), value.name = "Correlation")
  dfall$Factor1 <- factor(dfall$Factor1, levels = rownames(all_cor))
  dfall$Factor2 <- factor(dfall$Factor2, levels = colnames(all_cor))
  
  tab <- as.data.frame(table(mth[ord]), stringsAsFactors = FALSE)
  tab$cum <- cumsum(tab$Freq); bounds <- tab$cum + 0.5
  
  p_all <- ggplot2::ggplot(dfall, ggplot2::aes(x = Factor2, y = Factor1, fill = Correlation)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(limits = c(-1,1), midpoint = 0, name = "Corr") +
    ggplot2::geom_vline(xintercept = bounds, linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = bounds, linewidth = 0.3) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = "All methods — factor correlation (scores)", x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = plot_theme_base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"))
  
  # ---- 7) Optional AUCs vs truth (unchanged) ----
  normalize_ids_suffix <- function(xnames) normalize_ids(xnames)
  auc_results <- NULL
  if (!is.null(ground_truth_feature) || !is.null(ground_truth_sample)) {
    if (!requireNamespace("pROC", quietly = TRUE)) stop("To compute AUC, please install 'pROC'.")
    calc_auc <- function(scores, truth) {
      truth <- as.numeric(truth)
      if (!length(scores) || !length(truth) || length(scores) != length(truth)) return(NA_real_)
      roc <- pROC::roc(response = truth, predictor = scores, quiet = TRUE, direction = "<")
      as.numeric(pROC::auc(roc))
    }
    rows_auc <- list()
    if (!is.null(ground_truth_feature)) {
      for (m in available) {
        L <- loadings_list[[m]]
        s <- apply(abs(L), 1, max, na.rm = TRUE)
        idsL <- normalize_ids_suffix(names(s))
        idsT <- names(ground_truth_feature)
        if (!is.null(idsT)) idsT <- normalize_ids_suffix(idsT)
        common <- intersect(idsL, idsT)
        if (length(common)) {
          rows_auc[[length(rows_auc)+1]] <- data.frame(
            target="feature", method=m,
            AUC=calc_auc(scores = s[match(common, idsL)],
                         truth  = ground_truth_feature[match(common, idsT)]),
            stringsAsFactors = FALSE
          )
        }
      }
    }
    if (!is.null(ground_truth_sample)) {
      for (m in available) {
        S <- scores_list[[m]]
        ss <- apply(abs(S), 1, max, na.rm = TRUE)
        common <- intersect(names(ss), names(ground_truth_sample))
        if (length(common)) {
          rows_auc[[length(rows_auc)+1]] <- data.frame(
            target="sample", method=m,
            AUC=calc_auc(scores = ss[common], truth = ground_truth_sample[common]),
            stringsAsFactors = FALSE
          )
        }
      }
    }
    auc_results <- if (length(rows_auc)) do.call(rbind, rows_auc) else NULL
  }
  
  summary_table_scores   <- do.call(rbind, lapply(pair_summaries_scores, `[[`, "summary"))
  rownames(summary_table_scores) <- NULL
  summary_table_loadings <- do.call(rbind, lapply(pair_summaries_loadings, `[[`, "summary"))
  rownames(summary_table_loadings) <- NULL
  
  list(
    available_methods        = available,
    scores                   = scores_list,
    loadings                 = loadings_list,
    pairwise_scores          = pair_summaries_scores,
    pairwise_loadings        = pair_summaries_loadings,  # per-dataset
    pair_heatmaps_scores     = pair_heatmaps_scores,
    pair_heatmaps_loadings   = pair_heatmaps_loadings,   # per dataset
    pair_scatter_scores      = pair_scatter_scores,      # << NEW
    pair_scatter_loadings    = pair_scatter_loadings,    # << NEW (list per dataset)
    all_scores_heatmap       = p_all,
    summary_table_scores     = summary_table_scores,
    summary_table_loadings   = summary_table_loadings,
    auc_results              = auc_results,
    pairs_scores_by_factor               = pairs_scores_by_factor,           # NEW
    pairs_loadings_by_dataset_and_factor = pairs_loadings_by_dataset_and_factor  # NEW
    
  )
}

sim_bench <- benchmark_factor_methods(sim_benchmark_res)


#============================================
# Create scores dataframes
#============================================

# --------- 06-nov-2025 -----
get_signal_gt <- function(simulatedData, iter = "iteration_1") {
  # pull iteration object
  it <- simulatedData[[iter]]
  
  # sample indices that truly belong to factor 1 (GT) for omic 1 and omic 2
  idx_score1 <- it$indices_samples[[1]]
  idx_score2 <- it$indices_samples[[2]]
  
  # factor score matrices
  fs_score1 <- it$factor_scores[[1]]
  fs_score2 <- it$factor_scores[[2]]
  
  # build binary columns: 1 if row index is in GT indices, 0 otherwise
  F1_GT_score1 <- as.integer(seq_along(fs_score1) %in% idx_score1)
  F2_GT_score2 <- as.integer(seq_along(fs_score2) %in% idx_score2)
  
  # attach to factor score matrices
  fs_score1 <- cbind(fs_score1, F1_GT = F1_GT_score1)
  fs_score2 <- cbind(fs_score2, F2_GT = F2_GT_score2)
  
  # return a clean list
  list(
    score1 = fs_score1,
    score2 = fs_score2
  )
}

## use it
simGTscores <- get_signal_gt(simulatedData)
signal_simGTscores <- data.frame(simGTscores$score1,simGTscores$score2)
signal_simGTscores$sample <- paste0("sample_", seq_len(NROW(signal_simGTscores)))
names(signal_simGTscores) <- c("fs_score1","F1 (GT)", "fs_score1", "F2 (GT)", "sample")
signal_simGTscores <- signal_simGTscores %>%
  select(`F1 (GT)`, `F2 (GT)`, sample)

names(signal_simGTscores)
names(fabia_score_df)
names(mofa_score_df)
names(mfa_score_df)
names(gfa_score_df)

# ----------------------------------- 10 Nov 2025 ------------------------------------
## ------------------------------------------------------------
## helpers to grab factor columns
## ------------------------------------------------------------
## --- helpers --------------------------------------------------

# get method predicted factor, e.g. "F1 (FABIA)"
get_method_factor <- function(df, k = 1L, method_tag) {
  colname <- sprintf("F%d (%s)", k, method_tag)
  if (!colname %in% names(df)) {
    stop("Could not find column ", colname, " in:\n", paste(names(df), collapse = ", "))
  }
  df[[colname]]
}

# get GT factor, e.g. "F1 (GT)"
get_gt_factor <- function(df, k = 1L) {
  colname <- sprintf("F%d (GT)", k)
  if (!colname %in% names(df)) {
    stop("Could not find GT column ", colname, " in:\n", paste(names(df), collapse = ", "))
  }
  df[[colname]]
}


## ------------------------------------------------------------
## threshold rules (same, but we will pass abs(pred) in)
## ------------------------------------------------------------
thr_zmad <- function(x, c = 3, center = c("median","mean")) {
  center <- match.arg(center)
  ax <- x   # x already abs
  mu <- if (center == "median") median(ax) else mean(ax)
  mad0 <- mad(ax, center = mu, constant = 1)
  tau <- mu + c * mad0
  list(mask = ax >= tau, tau = tau)
}

thr_percentile <- function(x, q = 0.90) {
  ax <- x   # x already abs
  tau <- as.numeric(quantile(ax, probs = q, names = FALSE, type = 8))
  list(mask = ax >= tau, tau = tau)
}

thr_fdr <- function(x, alpha = 0.05, central_frac = 0.8) {
  # for scores, fdr is optional, but we keep it
  x0 <- x - median(x)
  rad <- as.numeric(quantile(abs(x0), central_frac))
  sigma0 <- sd(x0[abs(x0) <= rad]); if (!is.finite(sigma0) || sigma0 == 0) sigma0 <- sd(x0)
  z <- x0 / sigma0
  p <- 2 * pnorm(-abs(z))
  qv <- p.adjust(p, method = "BH")
  list(mask = qv <= alpha, q = qv, sigma0 = sigma0)
}

roll_mean <- function(x, d = 31) {
  L <- length(x)
  half <- floor(d / 2)
  out <- numeric(L)
  for (i in seq_len(L)) {
    lo <- max(1, i - half)
    hi <- min(L, i + half)
    out[i] <- mean(x[lo:hi])
  }
  out
}

thr_phi6_rolling_mean <- function(x, delta = 31) {
  R <- roll_mean(x, d = delta)
  b <- x > R
  if (!any(b)) return(list(mask = rep(FALSE, length(x)), tau = NA))
  tau <- mean(x[b])
  list(mask = b & (x >= tau), tau = tau)
}

## φ1: max2noise ----
## estimate a noise level from the central part of the data,
## then keep values whose magnitude is >= snr * noise
thr_phi1_max_to_noise <- function(x, central_frac = 0.8, snr = 2) {
  ax <- x  # assume already abs() upstream
  # central band to estimate noise
  rad <- as.numeric(quantile(ax, central_frac))
  noise_vals <- ax[ax <= rad]
  noise_sd <- sd(noise_vals)
  if (!is.finite(noise_sd) || noise_sd == 0) noise_sd <- sd(ax)
  tau <- snr * noise_sd
  list(mask = ax >= tau, tau = tau, noise_sd = noise_sd)
}

## φ2: top-k mean ----
## take the k largest values, compute their mean, and use that as threshold
# thr_phi2_topk_mean <- function(x, k = 9) {
#   ax <- x
#   n <- length(ax)
#   k <- min(k, n)
#   topk <- sort(ax, decreasing = TRUE)[seq_len(k)]
#   tau <- mean(topk)
#   list(mask = ax >= tau, tau = tau, k = k)
# }
## φ2 (scaled): take the k largest values of the *scaled* vector,
## compute their mean as tau, and threshold the *scaled* vector by tau.
thr_phi2_topk_mean <- function(x,
                               sigma,
                               type = c("loading","score"),
                               k_loading = 100,
                               k_score   = 9,
                               use_abs = TRUE,
                               na_rm = TRUE) {
  type <- match.arg(type)
  if (!is.numeric(x)) stop("'x' must be numeric.")
  if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma <= 0)
    stop("'sigma' must be a single positive finite numeric (noise variance).")
  
  v <- as.numeric(x)
  if (na_rm) v <- v[!is.na(v)]
  if (!length(v)) return(list(mask = logical(0), tau = NA_real_, mean = NA_real_, sd = NA_real_))
  
  ax <- if (use_abs) abs(v) else v
  denom <- if (type == "loading") sqrt(sigma) else sigma
  ax_norm <- ax / denom
  
  k <- if (type == "loading") min(k_loading, length(ax_norm)) else min(k_score, length(ax_norm))
  topk <- sort(ax_norm, decreasing = TRUE)[seq_len(k)]
  tau <- mean(topk)
  
  mu <- mean(ax_norm)
  s  <- stats::sd(ax_norm)
  
  list(mask = ax_norm >= tau, tau = tau, mean = mu, sd = s)
}

## φ3: quantile(pi) ----
## general quantile-based cutoff (true prop φ3)
# φ3 (quantile, π from truth): threshold at (1 - π)-quantile
thr_phi3_quantile_pi <- function(x,
                                 truth_vec = NULL,
                                 pi = 0.20,
                                 use_abs = TRUE,
                                 na_rm = TRUE,
                                 quantile_type = 8) {
  if (!is.numeric(x)) stop("'x' must be numeric.")
  xv <- as.numeric(x)
  
  # Keep index of non-NA entries for alignment with truth_vec and masking
  idx <- if (na_rm) which(!is.na(xv)) else seq_along(xv)
  if (!length(idx)) return(list(mask = logical(0), tau = NA_real_, pi = NA_real_))
  
  ax_full <- if (use_abs) abs(xv) else xv
  ax <- ax_full[idx]
  
  # Determine π from truth_vec when provided
  if (!is.null(truth_vec)) {
    n <- length(ax)
    
    count_truth <- NA_integer_
    
    if (is.logical(truth_vec) && length(truth_vec) == length(xv)) {
      count_truth <- sum(truth_vec[idx], na.rm = TRUE)
      
    } else if (is.numeric(truth_vec) && length(truth_vec) == length(xv) &&
               all(is.na(truth_vec) | truth_vec %in% c(0, 1))) {
      count_truth <- sum(truth_vec[idx] == 1, na.rm = TRUE)
      
    } else if (is.numeric(truth_vec) && all(truth_vec %% 1 == 0) &&
               all(truth_vec >= 1, na.rm = TRUE) &&
               all(truth_vec <= length(xv), na.rm = TRUE)) {
      # indices into original x
      count_truth <- length(intersect(as.integer(truth_vec), idx))
      
    } else {
      # Fallback: treat as a collection of positives
      count_truth <- length(truth_vec)
    }
    
    # Proportion π (clamped to [0,1])
    pi <- max(0, min(1, count_truth / n))
  }
  
  # Compute tau at (1 - π)-quantile of the working vector
  tau <- as.numeric(quantile(ax, probs = 1 - pi, names = FALSE, type = quantile_type))
  
  # Build mask on the same working vector (non-NA subset)
  mask_sub <- ax >= tau
  
  # Expand mask back to original length (respecting NA handling)
  mask <- rep(FALSE, length(xv))
  mask[idx] <- mask_sub
  
  list(mask = mask, tau = tau, pi = pi)
}


## φ4: normspread ----
## simple location + spread rule: mean + 1*sd
thr_phi4_norm_spread <- function(x) {
  ax <- x
  mu <- mean(ax)
  s  <- sd(ax)
  tau <- mu + s
  list(mask = ax >= tau, tau = tau, mean = mu, sd = s)
}

## φ5: fixed Q80 ----
## you wrote "fixed Q80" but called it thr_phi5_fixed_q10(); we'll keep your name,
## but use q = 0.80 so it actually does Q80.
thr_phi5_fixed_q20 <- function(x, q = 0.80) {
  ax <- x
  tau <- as.numeric(quantile(ax, probs = q, names = FALSE, type = 8))
  list(mask = ax >= tau, tau = tau, q = q)
}
## ------------------------------------------------------------
## metrics (unchanged)
## ------------------------------------------------------------
binary_metrics <- function(pred, truth) {
  TP <- sum(pred & truth)
  FP <- sum(pred & !truth)
  FN <- sum(!pred & truth)
  TN <- sum(!pred & !truth)
  sens <- ifelse(TP + FN == 0, NA, TP / (TP + FN))
  spec <- ifelse(TN + FP == 0, NA, TN / (TN + FP))
  prec <- ifelse(TP + FP == 0, NA, TP / (TP + FP))
  rec  <- sens
  f1   <- ifelse(is.na(prec) || is.na(rec) || (prec + rec == 0),
                 NA, 2 * prec * rec / (prec + rec))
  jacc <- ifelse(TP + FP + FN == 0, NA, TP / (TP + FP + FN))
  c(TP = TP, FP = FP, FN = FN, TN = TN,
    Sensitivity = sens, Specificity = spec,
    Precision = prec, Recall = rec, F1 = f1, Jaccard = jacc)
}

## ------------------------------------------------------------
## main evaluator: uses abs(pred), GT as given (0/1)
## ------------------------------------------------------------
evaluate_score_factor <- function(method_name, method_df, gt_df,
                                  id_col = "sample", k = 1L) {
  
  # merge on sample
  merged <- merge(method_df, gt_df, by = id_col, all = FALSE)
  
  # predicted scores → abs
  pred_vec <- get_method_factor(merged, k = k, method_tag = method_name)
  pred_abs <- abs(pred_vec)
  
  # GT already 0/1
  truth_vec <- get_gt_factor(merged, k = k)
  truth_bin <- truth_vec == 1
  
  # thresholds on abs(pred)
  masks <- list(
    # #"φ1 max2noise"  = thr_phi1_max_to_noise(pred_abs)$mask,
    # # "φ2 top-k mean (sd)" = thr_phi2_topk_mean(pred_abs, sigma = sqrt(7), type = "score", k_score = 15)$mask,
    # # "φ2 top-k mean (var)" = thr_phi2_topk_mean(pred_abs, sigma = 7, type = "score", k_score = 15)$mask,
    # # "φ2 top-1 mean (sd)" = thr_phi2_topk_mean(pred_abs, sigma = sqrt(7), type = "score", k_score = 1)$mask,
    # # "φ2 top-1 mean (var)" = thr_phi2_topk_mean(pred_abs, sigma = 7, type = "score", k_score = 1)$mask,
    # "φ1 normspread" = thr_phi4_norm_spread(pred_abs)$mask,
    # #"φ2 fixed Q80"  = thr_phi5_fixed_q20(pred_abs)$mask,
    # "φ2 true quantile"   = thr_phi3_quantile_pi(pred_abs, truth_bin)$mask,
    # "φ3 Percentile 80%" = thr_percentile(pred_abs, q = 0.80)$mask,
    # "φ4 Z/MAD"          = thr_zmad(pred_abs, c = 3)$mask,
    # "φ5 FDR 0.05"       = thr_fdr(pred_abs, alpha = 0.05)$mask,
    # "φ6 rollmean"    = thr_phi6_rolling_mean(pred_abs, delta = 31)$mask
    "τ1" = thr_phi4_norm_spread(pred_abs)$mask,
    "τ2"   = thr_phi3_quantile_pi(pred_abs, truth_bin)$mask,
    "τ3" = thr_percentile(pred_abs, q = 0.80)$mask,
    "τ4"          = thr_zmad(pred_abs, c = 3)$mask,
    "τ5"       = thr_fdr(pred_abs, alpha = 0.05)$mask,
    "τ6"    = thr_phi6_rolling_mean(pred_abs, delta = 31)$mask
  )
  #x, k, scale = "sd",  na_rm = na_rm
  out <- lapply(names(masks), function(rule_name) {
    met <- binary_metrics(masks[[rule_name]], truth_bin)
    c(method = method_name, factor = paste0("F", k), rule = rule_name, met)
  })
  do.call(rbind, out)
}


## ------------------------------------------------------------
## run for all your score dfs
## ------------------------------------------------------------
fabia_score_df=swap_F1_F2_names(fabia_score_df)
score_methods <- list(
  FABIA = fabia_score_df,
  MOFA  = mofa_score_df,
  MFA   = mfa_score_df,
  GFA   = gfa_score_df
)

all_res <- lapply(names(score_methods), function(m) {
  dfm <- score_methods[[m]]
  r1 <- evaluate_score_factor(m, dfm, signal_simGTscores, id_col = "sample", k = 1L)
  r2 <- evaluate_score_factor(m, dfm, signal_simGTscores, id_col = "sample", k = 2L)
  rbind(r1, r2)
})
all_res <- do.call(rbind, all_res)

# convert numeric columns
all_res <- as.data.frame(all_res, stringsAsFactors = FALSE)
num_cols <- c("TP","FP","FN","TN","Sensitivity","Specificity",
              "Precision","Recall","F1","Jaccard")
all_res[num_cols] <- lapply(all_res[num_cols], as.numeric)

all_res



# ---------- Packages ----------
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(patchwork)  # if you prefer cowplot, swap accordingly
# install.packages(c("patchwork","viridis")) if needed
library(viridis)

# ---------- Prep ----------
num_cols <- c("TP","FP","FN","TN","Sensitivity","Specificity","Precision","Recall","F1","Jaccard")
all_res[num_cols] <- lapply(all_res[num_cols], as.numeric)

all_res <- all_res %>%
  mutate(
    method = factor(method, levels = c("FABIA","MOFA","MFA","GFA")),
    factor = factor(factor, levels = c("F1","F2")),
    # Order the rules as they appear conceptually
    rule   = factor(rule,
                    levels = c("τ1","τ2","τ3",
                               "τ4","τ5","τ6"))
  )


######################################################################################################################################
# Packages
library(dplyr); library(tidyr); library(ggplot2); library(forcats); library(viridis); library(scales)

# Make numeric & set factor orders
num_cols <- c("TP","FP","FN","TN","Sensitivity","Specificity","Precision","Recall","F1","Jaccard")
all_res[num_cols] <- lapply(all_res[num_cols], as.numeric)

all_res <- all_res %>%
  mutate(
    method = factor(method, levels = c("FABIA","MOFA","MFA","GFA")),
    factor = factor(factor, levels = c("F1","F2")),
    rule   = factor(rule, levels = c("τ1","τ2","τ3",
                                     "τ4","τ5","τ6"))
  )
# Clean paper theme
theme_pub <- function(base_size = 10, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_line(linewidth = 0.2),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
}

# --------------------------- Figure A — Performance heatmaps (with values on tiles) ------------------------
perf_long <- all_res %>%
  pivot_longer(
    cols = c("Precision","F1"),
    names_to = "metric", values_to = "value"
  ) %>%
  mutate(metric = factor(metric, levels = c("F1","Precision")))

p_heat <- ggplot(perf_long, aes(x = rule, y = method, fill = value)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  scale_fill_viridis(option = "C", limits = c(0,1), breaks = seq(0,1,0.2),
                     name = "Score", oob = scales::squish) +
  facet_grid(factor ~ metric, scales = "free_x", space = "free_x") +
  labs(x = "Rule", y = "Method",
       title = "") + # Performance across methods, rules, and factors
  theme_pub()+
  theme(
    legend.position = "bottom",            # move legend below
    legend.direction = "horizontal"        # horizontal guide
  )

# Outline non-perfect tiles (for F1/Jaccard)
issues <- perf_long %>% 
  dplyr::filter(metric %in% c("F1","Precision") & value < 1) %>% 
  dplyr::distinct(factor, metric, method, rule)

if (nrow(issues) > 0) {
  p_heat <- p_heat +
    geom_tile(data = issues, aes(x = rule, y = method),
              fill = NA, color = "black", linewidth = 0.6, inherit.aes = FALSE)
}

p_heat


# Save
if (!dir.exists("figures_benchmark")) dir.create("figures_benchmark")
ggsave("figures_benchmark/FigureA_Heatmaps.pdf", p_heat, width = 12, height = 6, device = cairo_pdf)
ggsave("figures_benchmark/FigureA_Heatmaps.png", p_heat, width = 12, height = 6, dpi = 600)

#
# -------------------- Figure B — Jaccard lollipop (per factor) ----------------------------------------------
# Bright, colorblind-safe palette (Okabe–Ito)
pal_okabe <- c(
  FABIA = "#0072B2",  # blue
  MOFA  = "#D55E00",  # vermillion
  MFA   = "#009E73",  # bluish green
  GFA   = "#CC79A7"   # reddish purple
)

# consistent horizontal dodge for stems & points
pos <- position_dodge(width = 0.55)

p_lollipop <- ggplot(all_res, aes(x = rule, y = Jaccard)) +
  geom_segment(aes(x = rule, xend = rule, y = 0, yend = Jaccard, color = method),
               position = pos, linewidth = 0.9, alpha = 0.7, lineend = "round") +
  geom_point(aes(fill = method), position = pos, size = 3.3,
             shape = 21, color = "black", stroke = 0.4) +
  scale_color_manual(values = pal_okabe, name = "Method") +
  scale_fill_manual(values = pal_okabe,  name = "Method") +
  facet_wrap(~ factor, nrow = 1) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     labels = scales::number_format(accuracy = 0.01)) +
  labs(x = "Rule", y = "Jaccard",
       title = "") +
  theme_pub() +
  theme(legend.position = "right")
p_lollipop

ggsave("figures_benchmark/FigureB_JaccardLollipop.pdf", p_lollipop, width = 12, height = 6, device = cairo_pdf)
ggsave("figures_benchmark/FigureB_JaccardLollipop.png", p_lollipop, width = 10, height = 4, dpi = 600)

# -------------------- Figure C — Confusion composition (TP/FP/FN/TN as proportions) ----------------------
counts_long <- all_res %>%
  pivot_longer(cols = c(TP, FP, FN, TN), names_to = "type", values_to = "count") %>%
  mutate(type = factor(type, levels = c("TP","FP","FN","TN")))

p_counts <- ggplot(counts_long, aes(x = rule, y = count, fill = type)) +
  geom_bar(stat = "identity", position = "fill", color = "white", linewidth = 0.2) +
  facet_grid(factor ~ method) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_viridis(discrete = TRUE, option = "C", name = "Confusion") +
  labs(x = "Rule", y = "Proportion",
       title = "Confusion composition by method and rule",
       subtitle = "TP/FP/FN/TN normalized within factor × method × rule") +
  theme_pub() +
  theme(legend.position = "right")

p_counts

ggsave("figures_benchmark/FigureC_ConfusionComposition.pdf", p_counts, width = 11, height = 6, device = cairo_pdf)
ggsave("figures_benchmark/FigureC_ConfusionComposition.png", p_counts, width = 11, height = 6, dpi = 600)

# ------------------- Figure Sensitity specificity -----------
library(dplyr); library(tidyr); library(ggplot2); library(forcats); library(viridis); library(scales)

# factor ordering (adjust if you already set this)
all_res <- all_res %>%
  mutate(
    method = factor(method, levels = c("FABIA","MOFA","MFA","GFA")),
    factor = factor(factor, levels = c("F1","F2")),
    rule   = factor(rule, levels = c("τ1","τ2","τ3",
                                     "τ4","τ5","τ6"))
  )

theme_pub <- function() {
  theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 1),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
}

rates_long <- all_res %>%
  select(method, factor, rule, Sensitivity, Specificity) %>%
  pivot_longer(c(Sensitivity, Specificity), names_to = "metric", values_to = "value")

p_rates <- ggplot(rates_long, aes(x = rule, y = value, color = method, shape = method, group = method)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(size = 2) +
  geom_line(alpha = 0.4) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.2)) +
  scale_color_viridis_d(option = "C", end = 0.9) +
  labs(x = "Rule", y = "Rate", color = "Method", shape = "Method",
       title = "") + #Sensitivity (TPR) and Specificity (TNR)
  facet_grid(metric ~ factor) +
  theme_pub()

p_rates

# ---- Inputs you already have ----
# signal_simGTscores, fabia_score_df, mofa_score_df, mfa_score_df, gfa_score_df
# each with columns like "F1 (GT)", "F2 (GT)", ..., and a "sample" column

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# ---------- 1) Join everything by sample ----------
scores_all <- list(signal_simGTscores, fabia_score_df, mofa_score_df, mfa_score_df, gfa_score_df) %>%
  reduce(inner_join, by = "sample")

# ---------- 2) Helper to extract and rename a factor’s columns ----------
pick_factor <- function(df, factor_lab = c("F1","F2")) {
  factor_lab <- match.arg(factor_lab)
  # target column names
  wanted <- c(
    sprintf("%s (GT)",  factor_lab),
    sprintf("%s (FABIA)",factor_lab),
    sprintf("%s (MOFA)", factor_lab),
    sprintf("%s (MFA)",  factor_lab),
    sprintf("%s (GFA)",  factor_lab)
  )
  out <- df[, wanted, drop = FALSE]
  names(out) <- c("GT","FABIA","MOFA","MFA","GFA")
  out
}

X_F1 <- pick_factor(scores_all, "F1")
X_F2 <- pick_factor(scores_all, "F2")

# ---------- 3) Compute correlation matrices ----------
# method = "pearson" or "spearman"
corr_mat <- function(X, method = c("pearson","spearman")) {
  method <- match.arg(method)
  stats::cor(X, use = "pairwise.complete.obs", method = method)
}

# choose your metric here:
method_corr <- "pearson"  # or "spearman"

C1 <- corr_mat(X_F1, method_corr)
C2 <- corr_mat(X_F2, method_corr)

# ---------- 4) Long data for faceted plotting ----------
to_long <- function(M, facet_label) {
  as.data.frame(as.table(M)) |>
    rename(row = Var1, col = Var2, value = Freq) |>
    mutate(Factor = facet_label)
}
df_corr <- bind_rows(to_long(C1, "F1"), to_long(C2, "F2"))

# ---------- 5) Plot (side-by-side panels) ----------
theme_pub <- function() {
  theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal"
    )
}

p_corr <- ggplot(df_corr, aes(x = col, y = row, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", # 
                       limits = c(-1, 1), midpoint = 0,
                       name = ifelse(method_corr == "pearson","Pearson r","Spearman ρ")) +
  facet_wrap(~ Factor, nrow = 1) +
  labs(
    x = NULL, y = NULL,
    title = sprintf("%s correlation between GT and methods", 
                    ifelse(method_corr == "pearson","Pearson","Spearman"))
  ) +
  theme_pub() +
  guides(fill = guide_colorbar(title.position = "top", barwidth = unit(8, "cm"), barheight = unit(0.35, "cm")))

p_corr

# ---------- 6) Save (optional) ----------
if (!dir.exists("figures_benchmark")) dir.create("figures_benchmark")
ggsave(sprintf("figures_benchmark/Corr_%s_F1_F2.pdf", method_corr), p_corr, width = 9, height = 4.5, device = cairo_pdf)
ggsave(sprintf("figures_benchmark/Corr_%s_F1_F2.png", method_corr), p_corr, width = 9, height = 4.5, dpi = 600)

# -------------- 2. LOADING - PER OMIC -------------------

library(dplyr)
library(ggplot2)
library(reshape2)

# GT LOADINDS
get_loadsignal_gt <- function(simulatedData, iter = "iteration_1") {
  # pull iteration object
  it <- simulatedData[[iter]]
  
  # sample indices that truly belong to factor 1 (GT) for omic 1 and omic 2
  idx_omic1_loadfactor1 <- it$indices_features.1[[1]]
  idx_omic1_loadfactor2 <- it$indices_features.1[[2]]
  
  idx_omic2_loadfactor1 <- it$indices_features.2[[1]]

  # factor score matrices
  ls_factor1_omic1 <- it$loading_omic1$beta1
  ls_factor2_omic1 <- it$loading_omic1$beta1
  ls_factor1_omic2 <- it$loading_omic2$delta1
  
  # build binary columns: 1 if row index is in GT indices, 0 otherwise
  F1_GT_factor1_omic1 <- as.integer(seq_along(ls_factor1_omic1) %in% idx_omic1_loadfactor1)
  F2_GT_factor2_omic1 <- as.integer(seq_along(ls_factor2_omic1) %in% idx_omic1_loadfactor2)

  F1_GT_factor1_omic2 <- as.integer(seq_along(ls_factor1_omic2) %in% idx_omic2_loadfactor1)
  F2_GT_factor2_omic2 <- rep(0L, length(ls_factor1_omic2))
  
  # attach to factor score matrices
  ls_factor1_omic1 <- cbind(ls_factor1_omic1, F1_GT = F1_GT_factor1_omic1)
  ls_factor2_omic1 <- cbind(ls_factor2_omic1, F2_GT = F2_GT_factor2_omic1)

  ls_factor1_omic2 <- cbind(ls_factor1_omic2, F1_GT = F1_GT_factor1_omic2)
  ls_factor2_omic2 <- data.frame(F2_GT_factor2_omic2)
  colnames(ls_factor2_omic2) <- "F2_GT"
  
  # return a clean list
  list(
    factor1_omic1 = ls_factor1_omic1,
    factor2_omic1 = ls_factor2_omic1,
    factor1_omic2 = ls_factor1_omic2,
    factor2_omic2 = ls_factor2_omic2
  )
}

simGTloadings <- get_loadsignal_gt(simulatedData)
signal_simGTomi1loadings <- data.frame(simGTloadings$factor1_omic1,simGTloadings$factor2_omic1)
signal_simGTomi2loadings <- data.frame(simGTloadings$factor1_omic2,simGTloadings$factor2_omic2)

signal_simGTomi1loadings$feature <- paste0("feature_", seq_len(NROW(signal_simGTomi1loadings)))
signal_simGTomi2loadings$feature <- paste0("feature_", seq_len(NROW(signal_simGTomi2loadings)))

names(signal_simGTomi1loadings) <- c("ls_factor1_omic1","F1 (GT)", "ls_factor2_omic1", "F2 (GT)", "feature")
signal_simGTomi1loadings <- signal_simGTomi1loadings %>%
  select(`F1 (GT)`, `F2 (GT)`, feature)
names(signal_simGTomi2loadings) <- c("ls_factor1_omic2","F1 (GT)", "F2 (GT)", "feature")
signal_simGTomi2loadings <- signal_simGTomi2loadings %>%
  select(`F1 (GT)`, `F2 (GT)`, feature)

# LOADINGS
fabia_loading_df <- data.frame(sim_benchmark_res$FABIA$loading); fabia_loading_df$feature <- rownames(fabia_loading_df)
# rename columns that start with "F"
colnames(fabia_loading_df) <- sub("^F(\\d+)$", "F\\1 (FABIA)", colnames(fabia_loading_df))

mofa_loading_df <- data.frame(sim_benchmark_res$MOFA$loading); mofa_loading_df$feature <- rownames(mofa_loading_df)
# rename columns that start with "F"
colnames(mofa_loading_df) <- sub("^F(\\d+)$", "F\\1 (MOFA)", colnames(mofa_loading_df))

mfa_loading_df <- data.frame(sim_benchmark_res$MFA$loading); mfa_loading_df$feature <- rownames(mfa_loading_df)
# rename columns that start with "F"
colnames(mfa_loading_df) <- sub("^F(\\d+)$", "F\\1 (MFA)", colnames(mfa_loading_df))

gfa_loading_df <- data.frame(sim_benchmark_res$GFA$loading); gfa_loading_df$feature <- rownames(gfa_loading_df)
# rename columns that start with "F"
colnames(gfa_loading_df) <- sub("^F(\\d+)$", "F\\1 (GFA)", colnames(gfa_loading_df))

names(fabia_loading_df)
names(mofa_loading_df)
names(mfa_loading_df)
names(gfa_loading_df)
names(signal_simGTomi1loadings)
names(signal_simGTomi2loadings)

# Helper function: split loading df by prefix
split_by_prefix <- function(df, prefix_omic1 = "^omic1", prefix_omic2 = "^omic2") {
  df_omic1 <- df[grepl(prefix_omic1, df$feature), , drop = FALSE]
  df_omic2 <- df[grepl(prefix_omic2, df$feature), , drop = FALSE]
  list(omic1 = df_omic1, omic2 = df_omic2)
}


# if swap is needed
fabia_loading_df2 <- swap_F1_F2_names(fabia_loading_df)
gfa_loading_df2 <- swap_F1_F2_names(gfa_loading_df)

fabia_split <- split_by_prefix(fabia_loading_df2)
mofa_split  <- split_by_prefix(mofa_loading_df)
mfa_split   <- split_by_prefix(mfa_loading_df)
gfa_split   <- split_by_prefix(gfa_loading_df2)

# Keep only "feature_<number>" at the end of any string
.keep_feature_id <- function(x) sub(".*(feature_\\d+)$", "\\1", trimws(x))

# Clean a single data.frame/tibble if it has a 'feature' column
.clean_feature_col <- function(df, col = "feature") {
  if (is.data.frame(df) && col %in% names(df)) {
    df[[col]] <- .keep_feature_id(df[[col]])
  }
  df
}

# Clean an entire list (recurses if any elements are lists of dfs)
clean_feature_in_list <- function(x, col = "feature") {
  lapply(x, function(el) {
    if (is.list(el) && !is.data.frame(el)) {
      clean_feature_in_list(el, col)  # recurse for nested lists
    } else {
      .clean_feature_col(el, col)
    }
  })
}

loading_omic1_methods <- list(
  FABIA = fabia_split$omic1,
  MOFA  = mofa_split$omic1,
  MFA   = mfa_split$omic1,
  GFA   = gfa_split$omic1
)
loading_omic1_methods <- clean_feature_in_list(loading_omic1_methods)

loading_omic2_methods <- list(
  FABIA = fabia_split$omic2,
  MOFA  = mofa_split$omic2,
  MFA   = mfa_split$omic2,
  GFA   = gfa_split$omic2
)
loading_omic2_methods <- clean_feature_in_list(loading_omic2_methods)

## ------------------------------------------------------------
## main evaluator: uses abs(pred), GT as given (0/1)
## ------------------------------------------------------------
evaluate_loading <- function(method_name, method_df, gt_df,
                                  id_col = "feature", k = 1L) {
  
  # merge on sample
  merged <- merge(method_df, gt_df, by = id_col, all = FALSE)
  
  # predicted scores → abs
  pred_vec <- get_method_factor(merged, k = k, method_tag = method_name)
  pred_abs <- abs(pred_vec)
  
  # GT already 0/1
  truth_vec <- get_gt_factor(merged, k = k)
  truth_bin <- truth_vec == 1
  
  # thresholds on abs(pred)
  masks <- list(
    "τ1" = thr_phi4_norm_spread(pred_abs)$mask,
    "τ2"   = thr_phi3_quantile_pi(pred_abs, truth_bin)$mask,
    "τ3" = thr_percentile(pred_abs, q = 0.80)$mask,
    "τ4"          = thr_zmad(pred_abs, c = 3)$mask,
    "τ5"       = thr_fdr(pred_abs, alpha = 0.05)$mask,
    "τ6"    = thr_phi6_rolling_mean(pred_abs, delta = 31)$mask
  )
  out <- lapply(names(masks), function(rule_name) {
    met <- binary_metrics(masks[[rule_name]], truth_bin)
    c(method = method_name, factor = paste0("F", k), rule = rule_name, met)
  })
  do.call(rbind, out)
}

all_omic1_res <- lapply(names(loading_omic1_methods), function(m) {
  dfm <- loading_omic1_methods[[m]]
  r1 <- evaluate_loading(m, dfm, signal_simGTomi1loadings, id_col = "feature", k = 1L)
  r2 <- evaluate_loading(m, dfm, signal_simGTomi1loadings, id_col = "feature", k = 2L)
  rbind(r1, r2)
})
all_omic1_res <- do.call(rbind, all_omic1_res)

all_omic2_res <- lapply(names(loading_omic2_methods), function(m) {
  dfm <- loading_omic2_methods[[m]]
  r1 <- evaluate_loading(m, dfm, signal_simGTomi2loadings, id_col = "feature", k = 1L)
  r2 <- evaluate_loading(m, dfm, signal_simGTomi2loadings, id_col = "feature", k = 2L)
  rbind(r1, r2)
})
all_omic2_res <- do.call(rbind, all_omic2_res)

# convert numeric columns
all_omic1_res <- as.data.frame(all_omic1_res, stringsAsFactors = FALSE)
all_omic2_res <- as.data.frame(all_omic2_res, stringsAsFactors = FALSE)

num_cols <- c("TP","FP","FN","TN","Sensitivity","Specificity",
              "Precision","Recall","F1","Jaccard")
all_omic1_res[num_cols] <- lapply(all_omic1_res[num_cols], as.numeric)
all_omic2_res[num_cols] <- lapply(all_omic2_res[num_cols], as.numeric)

all_omic1_res
all_omic2_res



# ---------- Packages ----------
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(patchwork)  # if you prefer cowplot, swap accordingly
# install.packages(c("patchwork","viridis")) if needed
library(viridis)

# ---------- Prep ----------
num_cols <- c("TP","FP","FN","TN","Sensitivity","Specificity","Precision","Recall","F1","Jaccard")
all_omic2_res[num_cols] <- lapply(all_omic2_res[num_cols], as.numeric)

all_omic2_res <- all_omic2_res %>%
  mutate(
    method = factor(method, levels = c("FABIA","MOFA","MFA","GFA")),
    factor = factor(factor, levels = c("F1","F2")),
    # Order the rules as they appear conceptually
    rule   = factor(rule,
                    levels = c("τ1","τ2","τ3",
                               "τ4","τ5","τ6"))
  )


######################################################################################################################################
# Packages
library(dplyr); library(tidyr); library(ggplot2); library(forcats); library(viridis); library(scales)

# Make numeric & set factor orders
num_cols <- c("TP","FP","FN","TN","Sensitivity","Specificity","Precision","Recall","F1","Jaccard")
all_omic2_res[num_cols] <- lapply(all_omic2_res[num_cols], as.numeric)

all_omic2_res <- all_omic2_res %>%
  mutate(
    method = factor(method, levels = c("FABIA","MOFA","MFA","GFA")),
    factor = factor(factor, levels = c("F1","F2")),
    rule   = factor(rule, levels = c("τ1","τ2","τ3",
                                     "τ4","τ5","τ6"))
  )

# Clean paper theme
theme_pub <- function(base_size = 10, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_line(linewidth = 0.2),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
}

# --------------------------- Figure A — Performance heatmaps (with values on tiles) ------------------------
perf_long_omic2 <- all_omic2_res %>%
  pivot_longer(
    cols = c("Precision","F1"),
    names_to = "metric", values_to = "value"
  ) %>%
  mutate(metric = factor(metric, levels = c("F1","Precision")))

p_heat <- ggplot(perf_long_omic2, aes(x = rule, y = method, fill = value)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  scale_fill_viridis(option = "C", limits = c(0,1), breaks = seq(0,1,0.2),
                     name = "Loading", oob = scales::squish) +
  facet_grid(factor ~ metric, scales = "free_x", space = "free_x") +
  labs(x = "Rule", y = "Method",
       title = "") + # Performance across methods, rules, and factors
  theme_pub()+
  theme(
    legend.position = "bottom",            # move legend below
    legend.direction = "horizontal"        # horizontal guide
  )

# Outline non-perfect tiles (for F1/Jaccard)
issues <- perf_long_omic2 %>% 
  dplyr::filter(metric %in% c("F1","Precision") & value < 1) %>% 
  dplyr::distinct(factor, metric, method, rule)

if (nrow(issues) > 0) {
  p_heat <- p_heat +
    geom_tile(data = issues, aes(x = rule, y = method),
              fill = NA, color = "black", linewidth = 0.6, inherit.aes = FALSE)
}

p_heat


# Save
if (!dir.exists("figures_benchmark")) dir.create("figures_benchmark")
ggsave("figures_benchmark/FigureA_Heatmaps.pdf", p_heat, width = 12, height = 6, device = cairo_pdf)
ggsave("figures_benchmark/FigureA_Heatmaps.png", p_heat, width = 12, height = 6, dpi = 600)

#
# -------------------- Figure B — Jaccard lollipop (per factor) ----------------------------------------------
# Bright, colorblind-safe palette (Okabe–Ito)
pal_okabe <- c(
  FABIA = "#0072B2",  # blue
  MOFA  = "#D55E00",  # vermillion
  MFA   = "#009E73",  # bluish green
  GFA   = "#CC79A7"   # reddish purple
)

# consistent horizontal dodge for stems & points
pos <- position_dodge(width = 0.55)

p_lollipop <- ggplot(all_omic2_res, aes(x = rule, y = Jaccard)) +
  geom_segment(aes(x = rule, xend = rule, y = 0, yend = Jaccard, color = method),
               position = pos, linewidth = 0.9, alpha = 0.7, lineend = "round") +
  geom_point(aes(fill = method), position = pos, size = 3.3,
             shape = 21, color = "black", stroke = 0.4) +
  scale_color_manual(values = pal_okabe, name = "Method") +
  scale_fill_manual(values = pal_okabe,  name = "Method") +
  facet_wrap(~ factor, nrow = 1) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     labels = scales::number_format(accuracy = 0.01)) +
  labs(x = "Rule", y = "Jaccard",
       title = " ") +
  theme_pub() +
  theme(legend.position = "right")
p_lollipop

ggsave("figures_benchmark/FigureB_JaccardLollipop.pdf", p_lollipop, width = 12, height = 6, device = cairo_pdf)
ggsave("figures_benchmark/FigureB_JaccardLollipop.png", p_lollipop, width = 10, height = 4, dpi = 600)

# -------------------- Figure C — Confusion composition (TP/FP/FN/TN as proportions) ----------------------
counts_long <- all_omic1_res %>%
  pivot_longer(cols = c(TP, FP, FN, TN), names_to = "type", values_to = "count") %>%
  mutate(type = factor(type, levels = c("TP","FP","FN","TN")))

p_counts <- ggplot(counts_long, aes(x = rule, y = count, fill = type)) +
  geom_bar(stat = "identity", position = "fill", color = "white", linewidth = 0.2) +
  facet_grid(factor ~ method) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_viridis(discrete = TRUE, option = "C", name = "Confusion") +
  labs(x = "Rule", y = "Proportion",
       title = "Confusion composition by method and rule",
       subtitle = "TP/FP/FN/TN normalized within factor × method × rule") +
  theme_pub() +
  theme(legend.position = "right")

p_counts

ggsave("figures_benchmark/FigureC_ConfusionComposition.pdf", p_counts, width = 11, height = 6, device = cairo_pdf)
ggsave("figures_benchmark/FigureC_ConfusionComposition.png", p_counts, width = 11, height = 6, dpi = 600)

# ------------------- Figure Sensitity specificity -----------
library(dplyr); library(tidyr); library(ggplot2); library(forcats); library(viridis); library(scales)

# factor ordering (adjust if you already set this)
all_omic2_res1 <- all_omic2_res %>%
  mutate(
    method = factor(method, levels = c("FABIA","MOFA","MFA","GFA")),
    factor = factor(factor, levels = c("F1","F2")),
    rule   = factor(rule, levels = c("τ1","τ2","τ3",
                                     "τ4","τ5","τ6"))
  )

theme_pub <- function() {
  theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 1),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
}

rates_long <- all_omic2_res1 %>%
  select(method, factor, rule, Sensitivity, Specificity) %>%
  pivot_longer(c(Sensitivity, Specificity), names_to = "metric", values_to = "value")

p_rates <- ggplot(rates_long, aes(x = rule, y = value, color = method, shape = method, group = method)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(size = 2) +
  geom_line(alpha = 0.4) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.2)) +
  scale_color_viridis_d(option = "C", end = 0.9) +
  labs(x = "Rule", y = "Rate", color = "Method", shape = "Method",
       title = "") +
  facet_grid(metric ~ factor) +
  theme_pub()

p_rates

# ---- Inputs you already have ----
# signal_simGTscores, fabia_score_df, mofa_score_df, mfa_score_df, gfa_score_df
# each with columns like "F1 (GT)", "F2 (GT)", ..., and a "sample" column

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# ---------- 1) Join everything by sample ----------
loadings_omic2_all <- list(signal_simGTomi2loadings, loading_omic2_methods$FABIA, loading_omic2_methods$MOFA, loading_omic2_methods$MFA, loading_omic2_methods$GFA) %>%
  reduce(inner_join, by = "feature")

# ---------- 2) Helper to extract and rename a factor’s columns ----------
pick_factor <- function(df, factor_lab = c("F1","F2")) {
  factor_lab <- match.arg(factor_lab)
  # target column names
  wanted <- c(
    sprintf("%s (GT)",  factor_lab),
    sprintf("%s (FABIA)",factor_lab),
    sprintf("%s (MOFA)", factor_lab),
    sprintf("%s (MFA)",  factor_lab),
    sprintf("%s (GFA)",  factor_lab)
  )
  out <- df[, wanted, drop = FALSE]
  names(out) <- c("GT","FABIA","MOFA","MFA","GFA")
  out
}

X_F1 <- pick_factor(loadings_omic2_all, "F1")
X_F2 <- pick_factor(loadings_omic2_all, "F2")

# ---------- 3) Compute correlation matrices ----------
# method = "pearson" or "spearman"
corr_mat <- function(X, method = c("pearson","spearman")) {
  method <- match.arg(method)
  stats::cor(X, use = "pairwise.complete.obs", method = method)
}

# choose your metric here:
method_corr <- "pearson"  # or "spearman"

C1 <- corr_mat(X_F1, method_corr)
C2 <- corr_mat(X_F2, method_corr)

# ---------- 4) Long data for faceted plotting ----------
to_long <- function(M, facet_label) {
  as.data.frame(as.table(M)) |>
    rename(row = Var1, col = Var2, value = Freq) |>
    mutate(Factor = facet_label)
}
df_corr <- bind_rows(to_long(C1, "F1"), to_long(C2, "F2"))

# ---------- 5) Plot (side-by-side panels) ----------
theme_pub <- function() {
  theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 1),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal"
    )
}

p_corr <- ggplot(df_corr, aes(x = col, y = row, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", # 
                       limits = c(-1, 1), midpoint = 0,
                       name = ifelse(method_corr == "pearson","Pearson r","Spearman ρ")) +
  facet_wrap(~ Factor, nrow = 1) +
  labs(
    x = NULL, y = NULL,
    title = sprintf("%s correlation between GT and methods", 
                    ifelse(method_corr == "pearson","Pearson","Spearman"))
  ) +
  theme_pub() +
  guides(fill = guide_colorbar(title.position = "top", barwidth = unit(8, "cm"), barheight = unit(0.35, "cm")))

p_corr

# ---------- 6) Save (optional) ----------
if (!dir.exists("figures_benchmark")) dir.create("figures_benchmark")
ggsave(sprintf("figures_benchmark/Corr_%s_F1_F2.pdf", method_corr), p_corr, width = 9, height = 4.5, device = cairo_pdf)
ggsave(sprintf("figures_benchmark/Corr_%s_F1_F2.png", method_corr), p_corr, width = 9, height = 4.5, dpi = 600)
