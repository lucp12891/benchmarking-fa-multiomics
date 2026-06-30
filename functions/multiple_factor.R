# #  --------------------------------- MultiFactor: Simulation function ----------------------------
# multiple_factor <- function(n_features_one, n_features_two, n_samples, sigmas, iterations, n_factors){
#   set.seed(694)
#   
#   # Initialize lists to store omic data for all iterations and sigmas
#   all_omic_data <- list()
#   
#   for (iter in 1:iterations) {  # Iterate through the specified number of iterations
#     omic.one <- list()
#     omic.two <- list()
#     
#     for (k in 1:length(sigmas)) {  # Iterate through the list of sigmas
#       n_s <- n_samples
#       valid_data <- FALSE  # Flag to check if data passes FABIA condition
#       
#       #num_factor = num_factors
#       omic.one <- list()
#       omic.two <- list()
#       
#       # Generate first omic data
#       s_sig_s1 <- ceiling(n_s / 5.3)   # Start index for first range
#       s_sig_e1 <- ceiling(n_s / 2.8)   # End index for first range
#       s_sig_s2 <- 50                   # Start index for second range
#       s_sig_e2 <- 70                   # End index for second range
#       
#       alpha <- rnorm(n_s, 0, 0.05)
#       
#       # Assign indices for the first range
#       assigned_indices_samples1 <- sample(s_sig_s1:s_sig_e1, length(s_sig_s1:s_sig_e1))
#       
#       # Assign indices for the second range
#       assigned_indices_samples2 <- sample(s_sig_s2:s_sig_e2, length(s_sig_s2:s_sig_e2))
#       
#       # Combine the two sets of indices
#       assigned_indices_samples <- list(assigned_indices_samples1, assigned_indices_samples2)
#       
#       max_factors <- length(assigned_indices_samples)# Same as n_factors
#       
#       
#       # Generate random alpha values based on the max_factors
#       list_alphas <- list()
#       for (i in 1:max_factors) {
#         list_alphas[[paste0("alpha", i)]] <- rnorm(n_samples, 0, 0.05)
#       }
#       
#       # Assign corresponding values to alpha variables based on assigned_indices_samples
#       for (i in seq_along(assigned_indices_samples)) {
#         indices <- assigned_indices_samples[[i]]
#         list_alphas[[i]][indices] <- rnorm(length(indices), (3 + 0.5*i), 0.05)  # Adjust values dynamically
#       }
#       
#       list_gammas <- list_alphas
#       
#       # Add the separated list of the indices selected for omics of multiple features
#       # Initial vector
#       vector <- c(1:max_factors)
#       
#       # Specify the number of elements to sample for the shared vector
#       num_shared <- ceiling(runif(1, min = 1, max = length(vector) - 1))
#       
#       # Select the 'shared' elements
#       if (num_shared > 0) {
#         shared <- 1
#         # Remove the selected shared elements from the vector
#         remain_vector <- vector[!vector %in% shared]
#       }
#       
#       # Shuffle the remaining vector
#       if(length(remain_vector == 1)){
#         shuffled_remain = remain_vector
#       } else{
#         shuffled_remain <- sample(remain_vector)
#       }
#       # Split the shuffled elements equally or nearly equally between omic_one_unique and omic_two_unique
#       num_elements <- length(shuffled_remain)
#       split_point <- sample(1:num_elements, 1) # Randomly select a split point
#       
#       omic_one_unique <- shuffled_remain[0:split_point] # First part goes to omic_one_unique
#       omic_two_unique <- shuffled_remain[!shuffled_remain %in% omic_one_unique] # Second part goes to omic_two_unique
#       
#       # Print the results
#       list(
#         shared = shared,
#         omic_one_unique = omic_one_unique,
#         omic_two_unique = omic_two_unique
#       )
#       
#       # End of indices assignment
#       
#       # Assignment of factors
#       list_omic_one_factors = c(shared, omic_one_unique)
#       list_omic_two_factors = c(shared, omic_two_unique)
#       #all_indices = c(shared, omic_one_unique, omic_two_unique)
#       
#       factor_xtics <- list(
#         shared = shared,
#         omic_one_unique = omic_one_unique,
#         omic_two_unique = omic_two_unique
#       )
#       
#       # Generate first omic data
#       
#       # Define the first range of feature indices
#       f_sig_s1 <- 3650      # Start index for first range
#       f_sig_e1 <- 4000    # End index for first range
#       
#       # Define the second range of feature indices
#       f_sig_s2 <- 1600   # Start index for second range
#       f_sig_e2 <- 1950   # End index for second range
#       
#       # Assign indices for the first range
#       assigned_indices_features1 <- sample(f_sig_s1:f_sig_e1, length(f_sig_s1:f_sig_e1))
#       
#       # Assign indices for the second range
#       assigned_indices_features2 <- sample(f_sig_s2:f_sig_e2, length(f_sig_s2:f_sig_e2))
#       
#       # Combine the two sets of indices
#       assigned_indices_features <- list(assigned_indices_features1, assigned_indices_features2)
#       
#       # Generate random beta values based on the max_factors
#       list_betas <- list()
#       
#       # Shuffle the list
#       (shuffled_assigned_indices_features <- assigned_indices_features)# sample(assigned_indices_features)) # Do not reshuffle first
#       
#       # Create only one vector in list_betas using the random index
#       for (i in seq_along(assigned_indices_features)) {
#         list_betas[[i]] <- rnorm(n_features_one, 0, 0.05)
#       }
#       
#       # Loop through the ordered indices in list_omic_one_factors
#       for (i in seq_along(assigned_indices_features)) {
#         indices_ns <- assigned_indices_features[[i]]  # Get the corresponding indices from assigned_indices_features
#         if (length(indices_ns) > 0) {
#           # Initialize list_betas[[i]] if it's not already
#           if (is.null(list_betas[[i]])) {
#             list_betas[[i]] <- numeric(max(indices_ns))
#           }
#           # Assign values to the appropriate indices in list_betas
#           list_betas[[i]][indices_ns] <- rnorm(length(indices_ns), mean = (4.0 + 0.5 * i), sd = 0.05)
#         }
#       }
#       
#       # Renaming list_betas to desired factor number
#       for (i in seq_along(list_omic_one_factors)) {
#         names(list_betas)[i] <- paste0("beta", list_omic_one_factors[[i]])#list_omic_one_factors[[i]]
#       }
#       
#       pattern_alpha <-"^alpha\\d+"
#       matches_alpha <- grepl(pattern_alpha, names(list_alphas), perl = TRUE)
#       n_alpha = length(matches_alpha[matches_alpha == TRUE])
#       
#       pattern_beta <- "^beta\\d+"
#       matches_beta <- grepl(pattern_beta, names(list_betas), perl = TRUE)
#       n_beta = length(matches_beta[matches_beta == TRUE])
#       
#       # Initialize the data list to store the results
#       data_list_i <- list()
#       
#       # Sort the list by names in ascending order
#       list_betas <- list_betas[order(names(list_betas))]
#       list_alphas <- list_alphas[order(names(list_alphas))]
#       
#       # Extract the last numeric values from the names of the vectors
#       alphas_names <- as.numeric(gsub("[^0-9]", "", names(list_alphas)))
#       betas_names <- as.numeric(gsub("[^0-9]", "", names(list_betas)))
#       
#       # Find common numbers between the last values of the vector names
#       common_names <- intersect(alphas_names, betas_names)
#       
#       # Filter the vectors in each list based on the common names
#       list_alphas <- list_alphas[paste0("alpha", common_names)]
#       list_betas <- list_betas[paste0("beta", common_names)]
#       
#       # Loop through each alpha and beta combination
#       for (i in 1:min(length(list_alphas), length(list_betas))) {
#         data_i <- list_alphas[[i]] %*% t(list_betas[[i]][1:n_features_one])
#         data_list_i[[paste0("data.", i)]] <- data_i
#       }
#       
#       # Combine the results into a single data variable
#       data.1 <- Reduce(`+`, data_list_i)
#       #data.1 <- do.call(rbind, data_list_i)
#       
#       eps1 <- rnorm(n_samples * n_features_one, 0, sigmas) # noise
#       omic1_data <- matrix(data.1, n_samples, n_features_one) + matrix(eps1, n_samples, n_features_one) # signal + noise
#       colnames(omic1_data) <- paste0('omic1_feature_', seq_len(n_features_one))
#       rownames(omic1_data) <- paste0('sample_', seq_len(n_samples))
#       
#       
#       names(list_gammas) <- gsub("alpha", "gamma", names(list_gammas))
#       
#       # Extract the numeric parts from the names of list_gammas
#       numeric_part <- as.numeric(gsub("\\D", "", names(list_gammas)))
#       
#       # Retain elements where the numeric part of the name matches values in list_omic_two_factors
#       list_gammas <- list_gammas[numeric_part %in% list_omic_two_factors]
#       
#       # Empty list
#       list_deltas <- list()
#       
#       # Define the range for the feature indices
#       f_sig2_s1 <- 1   # Start index
#       f_sig2_e1 <- 300   # End index
#       
#       # Sample exactly 300 indices from the defined range
#       assigned_indices_features2 <- sample(f_sig2_s1:f_sig2_e1, 300)
#       assigned_indices_features_omic.two <- list(assigned_indices_features2)
#       
#       # Create only one vector in list_betas using the random index
#       for (i in seq_along(assigned_indices_features_omic.two)) {
#         list_deltas[[i]] <- rnorm(n_features_two, 0, 0.05)
#       }
#       
#       
#       for (i in seq_along(assigned_indices_features_omic.two)) {
#         indices_ns <- assigned_indices_features_omic.two[[i]]  # Get the corresponding indices from assigned_indices_features
#         
#         if (length(indices_ns) > 0) {
#           # Ensure list_deltas has a slot for the current index
#           if (length(list_deltas) < i) {
#             list_deltas[[i]] <- numeric(0)  # Initialize an empty numeric vector if it doesn't exist
#           }
#           
#           # Ensure list_deltas[[i]] has enough space for indices_ns
#           if (length(list_deltas[[i]]) < max(indices_ns)) {
#             length(list_deltas[[i]]) <- max(indices_ns)
#           }
#           
#           # Assign values to the appropriate indices in list_deltas
#           list_deltas[[i]][indices_ns] <- rnorm(length(indices_ns), mean = (5.0 + 0.0 * i), sd = 0.05)
#         }
#       }
#       
#       
#       # Renaming list_deltas to desired factor number
#       for (i in seq_along(list_omic_two_factors)) {
#         names(list_deltas)[i] <- paste0("delta", list_omic_two_factors[[i]])#list_omic_one_factors[[i]]
#       }
#       
#       pattern_gamma <-"^gamma\\d+"
#       matches_gamma <- grepl(pattern_gamma, names(list_gammas), perl = TRUE)
#       n_gamma = length(matches_alpha[matches_gamma == TRUE])
#       
#       pattern_delta <- "^delta\\d+"
#       matches_delta <- grepl(pattern_delta, names(list_deltas), perl = TRUE)
#       n_delta = length(matches_delta[matches_delta == TRUE])
#       
#       # Initialize the data list to store the results
#       data_list_j <- list()
#       
#       # Sort the list by names in ascending order
#       list_gammas <- list_gammas[order(names(list_gammas))]
#       list_deltas <- list_deltas[order(names(list_deltas))]
#       
#       # Extract the last numeric values from the names of the vectors
#       gammas_names <- as.numeric(gsub("[^0-9]", "", names(list_gammas)))
#       deltas_names <- as.numeric(gsub("[^0-9]", "", names(list_deltas)))
#       
#       # Find common numbers between the last values of the vector names
#       common_names <- intersect(gammas_names, deltas_names)
#       
#       # Filter the vectors in each list based on the common names
#       list_gammas <- list_gammas[paste0("gamma", common_names)]
#       list_deltas <- list_deltas[paste0("delta", common_names)]
#       
#       # Loop through each alpha and beta combination
#       for (j in 1:min(length(list_gammas), length(list_deltas))) {
#         data_j <- list_gammas[[j]] %*% t(list_deltas[[j]][1:n_features_two])
#         data_list_j[[paste0("data.", j)]] <- data_j
#       }
#       
#       # Combine the results into a single data variable
#       data.2 <- Reduce(`+`, data_list_j)
#       #data.2 <- do.call(rbind, data_list_j)
#       
#       eps2 <- rnorm(n_samples * n_features_two, 0, sigmas)
#       omic2_data <- matrix(data.2, n_samples, n_features_two) + matrix(eps2, n_samples, n_features_two) # signal + noise
#       colnames(omic2_data) <- paste0('omic2_feature_', seq_len(n_features_two))
#       rownames(omic2_data) <- paste0('sample_', seq_len(n_samples))
#       
#       # Concatenate datasets
#       concatenated_data <- cbind(omic2_data, omic1_data)
#       
#       # Save the validated data
#       omic.one[[k]] <- omic1_data
#       omic.two[[k]] <- omic2_data
#     }
#     
#     # Combine datasets for this iteration
#     simulated_datasets <- list(object.one = omic.one, object.two = omic.two)
#     concatenated_datasets <- list()
#     for (i in seq_along(simulated_datasets$object.one)) {
#       concatenated_datasets[[i]] <- cbind(simulated_datasets$object.one[[i]], simulated_datasets$object.two[[i]])
#     }
#     
#     # Save data for the current iteration
#     all_omic_data[[paste0("iteration_", iter)]] <- list(
#       concatenated_datasets = concatenated_datasets,
#       indices_features.1 = assigned_indices_features,
#       indices_samples = assigned_indices_samples,
#       indices_features.2 = assigned_indices_features_omic.two,
#       sample_sig_start1 = s_sig_s1,
#       sample_sig_end1 = s_sig_e1,
#       sample_sig_start2 = s_sig_s2,
#       sample_sig_end2 = s_sig_e2,
#       feature_sig_start = f_sig_s1,
#       feature_sig_end = f_sig_e1,
#       feature_sig_start = f_sig_s2,
#       feature_sig_end = f_sig_e2,
#       feature_sig2_start = f_sig2_s1,
#       feature_sig2_end = f_sig2_e1,
#       list_alphas = list_alphas,
#       list_gammas = list_gammas,
#       list_betas = list_betas,
#       list_deltas = list_deltas
#     )
#   }
#   
#   return(all_omic_data)
# }
# 
# #simulatedData <- multiple_factor(n_features_one = 4000, n_features_two = 3000, n_samples = 100, sigmas = 6, n_factors = 2, iterations = 1)
#    
# #dataset <- simulatedData$iteration_1$concatenated_datasets[[1]]
# #set.seed(657)
# #image(c(1:dim(dataset)[1]), c(1:dim(dataset)[2]), dataset, ylab = "Features", xlab = "Samples", font.lab = 2)
# #abline(h = 4000, col = "brown", lty = "dashed")  # lty = "dashed" sets the line type
# # # Add text annotations for the matrices
# # # Adjust x and y coordinates based on the location of the dark matrices
# # text(x = 28, y = 2550, labels = "shared factor", col = "black", cex = 0.8, font = 2)
# # text(x = 60.5, y = 5540, labels = "unique factor", col = "black", cex = 0.80, font = 2)
# # persp(c(1:dim(t(dataset))[2]), c(1:dim(t(dataset))[1]), dataset, theta = 325, phi = 15, col = "yellow", xlab = "", ylab = "Samples", zlab = " ")
# --------------------------------- MultiFactor: Simulation function ----------------------------
multiple_factor <- function(n_features_one, n_features_two, n_samples, sigmas, iterations, n_factors){
  set.seed(694)
  
  # Initialize lists to store omic data for all iterations and sigmas
  all_omic_data <- list()
  
  for (iter in 1:iterations) {  # Iterate through the specified number of iterations
    omic.one <- list()
    omic.two <- list()
    
    for (k in 1:length(sigmas)) {  # Iterate through the list of sigmas
      n_s <- n_samples
      
      #num_factor = num_factors
      omic.one <- list()
      omic.two <- list()
      
      # ----------------------- Sample indices and alphas -----------------------
      s_sig_s1 <- ceiling(n_s / 5.3)
      s_sig_e1 <- ceiling(n_s / 2.8)
      s_sig_s2 <- 50
      s_sig_e2 <- 70
      
      # Assign indices
      assigned_indices_samples1 <- sample(s_sig_s1:s_sig_e1, length(s_sig_s1:s_sig_e1))
      assigned_indices_samples2 <- sample(s_sig_s2:s_sig_e2, length(s_sig_s2:s_sig_e2))
      assigned_indices_samples <- list(assigned_indices_samples1, assigned_indices_samples2)
      
      max_factors <- length(assigned_indices_samples)
      
      # Generate alphas
      list_alphas <- list()
      for (i in 1:max_factors) {
        list_alphas[[paste0("alpha", i)]] <- rnorm(n_samples, 0, 0.05)
      }
      for (i in seq_along(assigned_indices_samples)) {
        indices <- assigned_indices_samples[[i]]
        list_alphas[[i]][indices] <- rnorm(length(indices), (3 + 0.5*i), 0.05)
      }
      list_gammas <- list_alphas
      
      # Shared vs unique factors
      vector <- c(1:max_factors)
      shared <- 1
      remain_vector <- vector[!vector %in% shared]
      shuffled_remain <- if(length(remain_vector) == 1) remain_vector else sample(remain_vector)
      split_point <- ifelse(length(shuffled_remain) > 0, sample(1:length(shuffled_remain), 1), 0)
      omic_one_unique <- if (split_point > 0) shuffled_remain[1:split_point] else integer(0)
      omic_two_unique <- shuffled_remain[!shuffled_remain %in% omic_one_unique]
      
      list_omic_one_factors <- c(shared, omic_one_unique)
      list_omic_two_factors <- c(shared, omic_two_unique)
      
      # ----------------------- Features and betas -----------------------
      f_sig_s1 <- 3650
      f_sig_e1 <- 4000
      f_sig_s2 <- 1600
      f_sig_e2 <- 1950
      assigned_indices_features1 <- sample(f_sig_s1:f_sig_e1, length(f_sig_s1:f_sig_e1))
      assigned_indices_features2 <- sample(f_sig_s2:f_sig_e2, length(f_sig_s2:f_sig_e2))
      assigned_indices_features <- list(assigned_indices_features1, assigned_indices_features2)
      
      list_betas <- list()
      for (i in seq_along(assigned_indices_features)) {
        list_betas[[i]] <- rnorm(n_features_one, 0, 0.05)
      }
      for (i in seq_along(assigned_indices_features)) {
        indices_ns <- assigned_indices_features[[i]]
        if (length(indices_ns) > 0) {
          list_betas[[i]][indices_ns] <- rnorm(length(indices_ns), mean = (4.0 + 0.5 * i), sd = 0.05)
        }
      }
      for (i in seq_along(list_omic_one_factors)) {
        names(list_betas)[i] <- paste0("beta", list_omic_one_factors[[i]])
      }
      
      # ----------------------- Build OMIC1 data -----------------------
      list_alphas <- list_alphas[order(names(list_alphas))]
      list_betas  <- list_betas[order(names(list_betas))]
      common_names <- intersect(
        as.numeric(gsub("[^0-9]", "", names(list_alphas))),
        as.numeric(gsub("[^0-9]", "", names(list_betas)))
      )
      list_alphas <- list_alphas[paste0("alpha", common_names)]
      list_betas  <- list_betas[paste0("beta",  common_names)]
      
      data_list_i <- list()
      for (i in 1:min(length(list_alphas), length(list_betas))) {
        vec_alpha <- as.numeric(list_alphas[[i]])
        vec_beta  <- as.numeric(list_betas[[i]][1:n_features_one])
        data_list_i[[paste0("data.", i)]] <- outer(vec_alpha, vec_beta)   # ✅ FIX
      }
      data.1 <- Reduce(`+`, data_list_i)
      
      eps1 <- rnorm(n_samples * n_features_one, 0, sigmas[k])
      omic1_data <- matrix(data.1, n_samples, n_features_one) + matrix(eps1, n_samples, n_features_one)
      colnames(omic1_data) <- paste0('omic1_feature_', seq_len(n_features_one))
      rownames(omic1_data) <- paste0('sample_', seq_len(n_samples))
      
      # ----------------------- Build OMIC2 data -----------------------
      names(list_gammas) <- gsub("alpha", "gamma", names(list_gammas))
      list_gammas <- list_gammas[order(names(list_gammas))]
      
      assigned_indices_features_omic.two <- list(sample(1:300, 300))
      list_deltas <- list()
      for (i in seq_along(assigned_indices_features_omic.two)) {
        vec <- rnorm(n_features_two, 0, 0.05)
        indices_ns <- assigned_indices_features_omic.two[[i]]
        vec[indices_ns] <- rnorm(length(indices_ns), mean = (5.0 + 0.0 * i), sd = 0.05)
        list_deltas[[i]] <- vec
      }
      for (i in seq_along(list_omic_two_factors)) {
        names(list_deltas)[i] <- paste0("delta", list_omic_two_factors[[i]])
      }
      list_deltas <- list_deltas[order(names(list_deltas))]
      common_names <- intersect(
        as.numeric(gsub("[^0-9]", "", names(list_gammas))),
        as.numeric(gsub("[^0-9]", "", names(list_deltas)))
      )
      list_gammas <- list_gammas[paste0("gamma", common_names)]
      list_deltas <- list_deltas[paste0("delta", common_names)]
      
      data_list_j <- list()
      for (j in 1:min(length(list_gammas), length(list_deltas))) {
        vec_gamma <- as.numeric(list_gammas[[j]])
        vec_delta <- as.numeric(list_deltas[[j]][1:n_features_two])
        data_list_j[[paste0("data.", j)]] <- outer(vec_gamma, vec_delta)   # ✅ FIX
      }
      data.2 <- Reduce(`+`, data_list_j)
      
      eps2 <- rnorm(n_samples * n_features_two, 0, sigmas[k])
      omic2_data <- matrix(data.2, n_samples, n_features_two) + matrix(eps2, n_samples, n_features_two)
      colnames(omic2_data) <- paste0('omic2_feature_', seq_len(n_features_two))
      rownames(omic2_data) <- paste0('sample_', seq_len(n_samples))
      
      # Concatenate datasets
      concatenated_data <- cbind(omic2_data, omic1_data)
      
      omic.one[[k]] <- omic1_data
      omic.two[[k]] <- omic2_data
    }
    
    concatenated_datasets <- lapply(seq_along(omic.one), function(i) cbind(omic.one[[i]], omic.two[[i]]))
    
    all_omic_data[[paste0("iteration_", iter)]] <- list(
      concatenated_datasets = concatenated_datasets,
      indices_features.1 = assigned_indices_features,
      indices_samples = assigned_indices_samples,
      indices_features.2 = assigned_indices_features_omic.two,
      list_alphas = list_alphas,
      list_gammas = list_gammas,
      list_betas = list_betas,
      list_deltas = list_deltas
    )
  }
  
  return(all_omic_data)
}
