#-------------------------------------------------------------------------------
library(dplyr)

# Step 1: List all relevant .rds files
file_paths <- list.files(path = "C:/Users/bosangir/Downloads/test_sim/Run_12_05_2025/simulation_results_hpc", 
                         pattern = "sim_omicResults_single\\.factor_sigma_\\d+\\.rds", 
                         full.names = TRUE)

# Initialize a list to store data for each varphi
varphi_data <- vector("list", 12)  # Assuming 15 possible varphi values

# Process each file
for (file_path in file_paths) {
  # Read the .rds file
  sim_data <- readRDS(file_path)
  
  # Extract the file identifier (X) from the current file path
  X <- as.numeric(gsub(".*factor_sigma_(\\d+)\\.rds", "\\1", basename(file_path)))
  
  # Extract `jaccard_results`
  jaccard_results <- sim_data[["jaccard_results"]]
  
  # If `jaccard_results` is NULL, skip this file
  if (is.null(jaccard_results)) {
    warning(paste("No valid `jaccard_results` found in:", file_path))
    next
  }
  
  # Process each key in `jaccard_results`
  for (key in names(jaccard_results)) {
    # Extract varphi index and iteration
    varphi_idx <- as.integer(gsub(".*_varphi_(\\d+).*", "\\1", key))
    iteration <- as.integer(gsub(".*iteration_(\\d+)_.*", "\\1", key))
    
    # Access the nested data for the current key
    varphi_content <- jaccard_results[[key]]
    
    # Create a data frame for the current key
    df <- data.frame(
      X = X,
      iteration = iteration,
      varphi = varphi_idx,
      ji_true_fabia = varphi_content$ji_true[1],
      ji_true_mofa = varphi_content$ji_true[2],
      ji_true_mfa = varphi_content$ji_true[3],
      ji_true_gfa = varphi_content$ji_true[4],
      ji_omic_one_fabia = varphi_content$ji_omic_one[1],
      ji_omic_one_mofa = varphi_content$ji_omic_one[2],
      ji_omic_one_mfa = varphi_content$ji_omic_one[3],
      ji_omic_one_gfa = varphi_content$ji_omic_one[4],
      ji_omic_two_fabia = varphi_content$ji_omic_two[1],
      ji_omic_two_mofa = varphi_content$ji_omic_two[2],
      ji_omic_two_mfa = varphi_content$ji_omic_two[3],
      ji_omic_two_gfa = varphi_content$ji_omic_two[4],
      ji_smp_fabia = varphi_content$ji_samples[1],
      ji_smp_mofa = varphi_content$ji_samples[2],
      ji_smp_mfa = varphi_content$ji_samples[3],
      ji_smp_gfa = varphi_content$ji_samples[4],
      ji_matrix_fabia = varphi_content$ji_matrix[1],
      ji_matrix_mofa = varphi_content$ji_matrix[2],
      ji_matrix_mfa = varphi_content$ji_matrix[3],
      ji_matrix_gfa = varphi_content$ji_matrix[4]
    )
    
    # Store the data in the respective varphi index
    if (is.null(varphi_data[[varphi_idx]])) {
      varphi_data[[varphi_idx]] <- df
    } else {
      varphi_data[[varphi_idx]] <- bind_rows(varphi_data[[varphi_idx]], df)
    }
  }
}

varphi_one <- varphi_data[[1]]
varphi_two <- varphi_data[[2]]
varphi_three <- varphi_data[[3]]
varphi_four <- varphi_data[[4]]
varphi_five <- varphi_data[[5]]
varphi_six <- varphi_data[[6]]
varphi_seven <- varphi_data[[7]]
varphi_eight <- varphi_data[[8]]
varphi_nine <- varphi_data[[9]]
varphi_ten <- varphi_data[[10]]
varphi_eleven <- varphi_data[[11]]
varphi_twelve <- varphi_data[[12]]

# ------------------------------------ JACCARD COMPARISONS RESULTS ------------------------------------

library(dplyr)

# Step 1: List all relevant .rds files #
file_paths <- list.files(path = "C:/Users/bosangir/Downloads/test_sim/Run_12_05_2025/simulation_results_hpc", 
                         pattern = "sim_omicResults_single\\.factor_sigma_\\d+\\.rds", 
                         full.names = TRUE)

# Initialize a list to store data for each varphi
varphi_data_comp <- vector("list", 12)  # Assuming 12 possible varphi values

# Process each file
for (file_path in file_paths) {
  # Read the .rds file
  sim_data <- readRDS(file_path)
  
  # Extract the file identifier (X) from the current file path
  X <- as.numeric(gsub(".*factor_sigma_(\\d+)\\.rds", "\\1", basename(file_path)))
  
  # Extract `jaccard_results`
  jaccard_comparison_results <- sim_data[["jaccard_comparison_results"]]
  
  # If `jaccard_results` is NULL, skip this file
  if (is.null(jaccard_comparison_results)) {
    warning(paste("No valid `jaccard_comparisons_results` found in:", file_path))
    next
  }
  
  # Process each key in `jaccard_results`
  for (key in names(jaccard_comparison_results)) {
    # Extract varphi index and iteration
    varphi_idx <- as.integer(gsub(".*_varphi_(\\d+).*", "\\1", key))
    iteration <- as.integer(gsub(".*iteration_(\\d+)_.*", "\\1", key))
    
    # Access the nested data for the current key
    varphi_content <- jaccard_comparison_results[[key]]
    
    # Create a data frame for the current key
    df <- data.frame(
      X = X,
      iteration = iteration,
      varphi = varphi_idx,
      ji_cp_fabia_mofa = varphi_content$ji_true_pairs$fabia_mofa_ji_phi,
      ji_cp_fabia_mfa = varphi_content$ji_true_pairs$fabia_mfa_ji_phi,
      ji_cp_fabia_gfa = varphi_content$ji_true_pairs$fabia_gfa_ji_phi,
      ji_cp_mofa_mfa = varphi_content$ji_true_pairs$mofa_mfa_ji_phi,
      ji_cp_mofa_gfa = varphi_content$ji_true_pairs$mofa_gfa_ji_phi,
      ji_cp_mfa_gfa = varphi_content$ji_true_pairs$mfa_gfa_ji_phi,
      
      ji_omic_one_fabia_mofa = varphi_content$ji_omic_one_pairs$fabia_mofa_ji_phi,
      ji_omic_one_fabia_mfa = varphi_content$ji_omic_one_pairs$fabia_mfa_ji_phi,
      ji_omic_one_fabia_gfa = varphi_content$ji_omic_one_pairs$fabia_gfa_ji_phi,
      ji_omic_one_mofa_mfa = varphi_content$ji_omic_one_pairs$mofa_mfa_ji_phi,
      ji_omic_one_mofa_gfa = varphi_content$ji_omic_one_pairs$mofa_gfa_ji_phi,
      ji_omic_one_mfa_gfa = varphi_content$ji_omic_one_pairs$mfa_gfa_ji_phi,
      
      ji_omic_two_fabia_mofa = varphi_content$ji_omic_two_pairs$fabia_mofa_ji_phi,
      ji_omic_two_fabia_mfa = varphi_content$ji_omic_two_pairs$fabia_mfa_ji_phi,
      ji_omic_two_fabia_gfa = varphi_content$ji_omic_two_pairs$fabia_gfa_ji_phi,
      ji_omic_two_mofa_mfa = varphi_content$ji_omic_two_pairs$mofa_mfa_ji_phi,
      ji_omic_two_mofa_gfa = varphi_content$ji_omic_two_pairs$mofa_gfa_ji_phi,
      ji_omic_two_mfa_gfa = varphi_content$ji_omic_two_pairs$mfa_gfa_ji_phi,
      
      ji_smp_fabia_mofa = varphi_content$ji_smp_pairs$fabia_mofa_ji_phi,
      ji_smp_fabia_mfa = varphi_content$ji_smp_pairs$fabia_mfa_ji_phi,
      ji_smp_fabia_gfa = varphi_content$ji_smp_pairs$fabia_gfa_ji_phi,
      ji_smp_mofa_mfa = varphi_content$ji_smp_pairs$mofa_mfa_ji_phi,
      ji_smp_mofa_gfa = varphi_content$ji_smp_pairs$mofa_gfa_ji_phi,
      ji_smp_mfa_gfa = varphi_content$ji_smp_pairs$mfa_gfa_ji_phi,      
      
      ji_matrix_fabia_mofa = varphi_content$ji_matrix_pairs[1],
      ji_matrix_fabia_mfa = varphi_content$ji_matrix_pairs[2],
      ji_matrix_fabia_gfa = varphi_content$ji_matrix_pairs[3],
      ji_matrix_mofa_mfa = varphi_content$ji_matrix_pairs[4],
      ji_matrix_mofa_gfa = varphi_content$ji_matrix_pairs[5],
      ji_matrix_mfa_gfa = varphi_content$ji_matrix_pairs[6]
    )
    
    # Store the data in the respective varphi index
    if (is.null(varphi_data_comp[[varphi_idx]])) {
      varphi_data_comp[[varphi_idx]] <- df
    } else {
      varphi_data_comp[[varphi_idx]] <- bind_rows(varphi_data_comp[[varphi_idx]], df)
    }
  }
}

varphi_one_cp <- varphi_data_comp[[1]]
varphi_two_cp <- varphi_data_comp[[2]]
varphi_three_cp <- varphi_data_comp[[3]]
varphi_four_cp <- varphi_data_comp[[4]]
varphi_five_cp <- varphi_data_comp[[5]]
varphi_six_cp <- varphi_data_comp[[6]]
varphi_seven_cp <- varphi_data_comp[[7]]
varphi_eight_cp <- varphi_data_comp[[8]]
varphi_nine_cp <- varphi_data_comp[[9]]
varphi_ten_cp <- varphi_data_comp[[10]]
varphi_eleven_cp <- varphi_data_comp[[11]]
varphi_twelve_cp <- varphi_data_comp[[12]]

# PERFORMANCE MEASURES

# Step 1: List all relevant .rds files
file_paths <- list.files(path = "C:/Users/bosangir/Downloads/test_sim/Run_12_05_2025/simulation_results_hpc", 
                         pattern = "sim_omicResults_single\\.factor_sigma_\\d+\\.rds", 
                         full.names = TRUE)

# Initialize a list to store data for each varphi
varphi_metric_data <- vector("list", 12) # Assuming 12 varphi values

# Process each file
for (file_path in file_paths) {
  # Read the .rds file
  sim_data <- readRDS(file_path)
  
  # Extract the file identifier X
  #X <- gsub(".*sim_omicResults_single.factor_sigma_(\\d+)\\.rds", "\\1", basename(file_path))
  X <- as.numeric(gsub(".*factor_sigma_(\\d+)\\.rds", "\\1", basename(file_path)))
  
  # Construct the path dynamically using paste
  path_to_metrics1 <- paste0("sim_omicResults_single.factor_sigma_", X)

  # Navigate to `jaccard_results`
  metrics_results1 <- sim_data$per_measures_results#sim_data[[path_to_metrics1]]$per_measures_results
  
  # Initialize an empty list to store combined metrics results
  combined_metrics_results <- list()
  
  # Check and merge `metrics_results` from both paths
  if (!path_to_metrics1 %in% names(sim_data)) {
    warning(paste("Missing expected object:", path_to_metrics1, "in file:", file_path))
    next
  }

  # If no valid results are found, skip to the next file
  if (length(combined_metrics_results) == 0) {
    warning(paste("No valid metrics_results found for file:", file_path))
    next
  }
  
  # Process each top-level key, e.g., `variance_9_iteration_1_varphi_1`
  for (key in names(combined_metrics_results)) {# Process each top-level key, e.g., `variance_9_iteration_1_varphi_1`
    # Extract the varphi index from the key
    varphi_idx <- as.integer(gsub(".*_varphi_(\\d+)", "\\1", key))
    
    # Extract the iteration index from the key
    iteration <- gsub(".*iteration_(\\d+)_.*", "\\1", key)
    
    # Access the nested data
    varphi_content <- combined_metrics_results[[key]]
    
    # Skip if varphi_content or required fields are NULL
    if (is.null(varphi_content)) {
      warning(paste("Skipping key due to NULL varphi_content:", key))
      next
    }
    
    # Helper function to safely extract values or return NA
    safe_extract <- function(data, path, index) {
      if (!is.null(data[[path]]) && length(data[[path]]) >= index) {
        return(data[[path]][index])
      }
      return(NA) # Return NA if the field doesn't exist
    }
    
    # Create a data frame with required columns
    df <- data.frame(
      X = X,
      iteration = iteration,
      varphi = varphi_idx,
      #varphi_content$per_measures_true$per_measures_true[7],#,$MOFA.Sensitivity true_sens_mofa = 
      #varphi_content$per_measures_true$per_measures_true[17],#$FABIA.Sensitivity true_sens_fabia = 
      #varphi_content$per_measures_true$per_measures_true[27],#$MFA.Sensitivity true_sens_mfa = 
      #varphi_content$per_measures_true$per_measures_true[37],#$GFA.Sensitivity true_sens_gfa = 
      # 
      omic_one_sens_mofa = varphi_content$per_measures_omic_one$per_measures_true[7],#$MOFA.Sensitivity
      omic_one_sens_fabia = varphi_content$per_measures_omic_one$per_measures_true[8],#$FABIA.Sensitivity
      omic_one_sens_mfa = varphi_content$per_measures_omic_one$per_measures_true[9],#$MFA.Sensitivity
      omic_one_sens_gfa = varphi_content$per_measures_omic_one$per_measures_true[10],#$GFA.Sensitivity
      # 
      omic_two_sens_mofa = varphi_content$per_measures_omic_two$per_measures_true[7],#$MOFA.Sensitivity
      omic_two_sens_fabia = varphi_content$per_measures_omic_two$per_measures_true[8],#$FABIA.Sensitivity
      omic_two_sens_mfa = varphi_content$per_measures_omic_two$per_measures_true[9],#$MFA.Sensitivity
      omic_two_sens_gfa = varphi_content$per_measures_omic_two$per_measures_true[10],#$GFA.Sensitivity,
      # 
      smp_sens_mofa = varphi_content$per_measures_smp[7],#$MOFA.Sensitivity
      smp_sens_fabia = varphi_content$per_measures_smp[8],#$FABIA.Sensitivity
      smp_sens_mfa = varphi_content$per_measures_smp[9],#$MFA.Sensitivity
      smp_sens_gfa = varphi_content$per_measures_smp[10],#$GFA.Sensitivity
      #  
      #varphi_content$per_measures_true$per_measures_true[8],#
      #varphi_content$per_measures_true$per_measures_true[18],#
      #varphi_content$per_measures_true$per_measures_true[28],#
      #varphi_content$per_measures_true$per_measures_true[38],#
      # 
      omic_one_spec_mofa = varphi_content$per_measures_omic_one$per_measures_true[7],#$MOFA.Specificity
      omic_one_spec_fabia = varphi_content$per_measures_omic_one$per_measures_true[8],#$FABIA.Specificity
      omic_one_spec_mfa = varphi_content$per_measures_omic_one$per_measures_true[9],#$MFA.Specificity
      omic_one_spec_gfa = varphi_content$per_measures_omic_one$per_measures_true[10],#$GFA.Specificity
      #
      omic_two_spec_mofa = varphi_content$per_measures_omic_two$per_measures_true[7],#$MOFA.Specificity
      omic_two_spec_fabia = varphi_content$per_measures_omic_two$per_measures_true[8],#$FABIA.Specificity
      omic_two_spec_mfa = varphi_content$per_measures_omic_two$per_measures_true[9],#$MFA.Specificity
      omic_two_spec_gfa = varphi_content$per_measures_omic_two$per_measures_true[10],#$GFA.Specificity
      #
      smp_spec_mofa = varphi_content$per_measures_smp$MOFA.Specificity[7],
      smp_spec_fabia = varphi_content$per_measures_smp$FABIA.Specificity[8],
      smp_spec_mfa = varphi_content$per_measures_smp$MFA.Specificity[9],
      smp_spec_gfa = varphi_content$per_measures_smp$GFA.Specificity[10],
      # 
      #varphi_content$per_measures_true$per_measures_true[9],
      #varphi_content$per_measures_true$per_measures_true[19],
      #varphi_content$per_measures_true$per_measures_true[29],
      #varphi_content$per_measures_true$per_measures_true[39],
      #
      omic_one_acc_mofa = varphi_content$per_measures_omic_one$per_measures_true$MOFA.Accuracy[7],
      omic_one_acc_fabia = varphi_content$per_measures_omic_one$per_measures_true$FABIA.Accuracy[8],
      omic_one_acc_mfa = varphi_content$per_measures_omic_one$per_measures_true$MFA.Accuracy[9],
      omic_one_acc_gfa = varphi_content$per_measures_omic_one$per_measures_true$GFA.Accuracy[10],
      #
      omic_two_acc_mofa = varphi_content$per_measures_omic_two$per_measures_true$MOFA.Accuracy[7],
      omic_two_acc_fabia = varphi_content$per_measures_omic_two$per_measures_true$FABIA.Accuracy[8],
      omic_two_acc_mfa = varphi_content$per_measures_omic_two$per_measures_true$MFA.Accuracy[9],
      omic_two_acc_gfa = varphi_content$per_measures_omic_two$per_measures_true$GFA.Accuracy[10],
      #
      smp_acc_mofa = varphi_content$per_measures_smp$MOFA.Accuracy[7],
      smp_acc_fabia = varphi_content$per_measures_smp$FABIA.Accuracy[8],
      smp_acc_mfa = varphi_content$per_measures_smp$MFA.Accuracy[9],
      smp_acc_gfa = varphi_content$per_measures_smp$GFA.Accuracy[10],
      # 
      # varphi_content$per_measures_true$per_measures_true[10],
      # varphi_content$per_measures_true$per_measures_true[20],
      # varphi_content$per_measures_true$per_measures_true[30],
      # varphi_content$per_measures_true$per_measures_true[40],
      #
      omic_one_auc_mofa = varphi_content$per_measures_omic_one$per_measures_true$MOFA.AUC[7],
      omic_one_auc_fabia = varphi_content$per_measures_omic_one$per_measures_true$FABIA.AUC[8],
      omic_one_auc_mfa = varphi_content$per_measures_omic_one$per_measures_true$MFA.AUC[3],
      omic_one_auc_gfa = varphi_content$per_measures_omic_one$per_measures_true$GFA.AUC[4],
      #
      omic_two_auc_mofa = varphi_content$per_measures_omic_two$per_measures_true$MOFA.AUC[7],
      omic_two_auc_fabia = varphi_content$per_measures_omic_two$per_measures_true$FABIA.AUC[8],
      omic_two_auc_mfa = varphi_content$per_measures_omic_two$per_measures_true$MFA.AUC[3],
      omic_two_auc_gfa = varphi_content$per_measures_omic_two$per_measures_true$GFA.AUC[4],
      #
      smp_auc_mofa = varphi_content$per_measures_smp$MOFA.AUC[7],
      smp_auc_fabia = varphi_content$per_measures_smp$FABIA.AUC[8],
      smp_auc_mfa = varphi_content$per_measures_smp$MFA.AUC[3],
      smp_auc_gfa = varphi_content$per_measures_smp$GFA.AUC[4]
    )
    
    # Skip if the data frame has zero rows
    if (nrow(df) == 0) {
      warning(paste("Skipping key due to empty data frame:", key))
      next
    }
    
    # Store the data in the respective varphi index
    if (is.null(varphi_metric_data[[varphi_idx]])) {
      varphi_metric_data[[varphi_idx]] <- df
    } else {
      varphi_metric_data[[varphi_idx]] <- bind_rows(varphi_metric_data[[varphi_idx]], df)
    }
  }
  
}

# Initialize a list to store data for each varphi
varphi_metric_data <- vector("list", 12) # Assuming 12 varphi values

# Process each file
for (file_path in file_paths) {
  # Read the .rds file
  sim_data <- readRDS(file_path)
  
  # Extract the sigma value from filename
  X <- as.numeric(gsub(".*factor_sigma_(\\d+)\\.rds", "\\1", basename(file_path)))
  
  # Access the per_measures_results list directly
  metrics_results1 <- sim_data$per_measures_results
  
  if (is.null(metrics_results1) || length(metrics_results1) == 0) {
    warning(paste("No per_measures_results found for file:", file_path))
    next
  }
  
  # Process each top-level key, e.g., variance_9_iteration_1_varphi_1
  for (key in names(metrics_results1)) {
    varphi_idx <- as.integer(gsub(".*_varphi_(\\d+)", "\\1", key))
    iteration <- gsub(".*iteration_(\\d+)_.*", "\\1", key)
    
    varphi_content <- metrics_results1[[key]]
    
    if (is.null(varphi_content)) {
      warning(paste("Skipping key due to NULL varphi_content:", key))
      next
    }
    
    # Create a data frame with required columns
    df <- data.frame(
      X = X,
      iteration = iteration,
      varphi = varphi_idx,
      
      # Sensitivities
      omic_one_sens_mofa = varphi_content$per_measures_omic_one$per_measures_true$Sensitivity[[2]],
      omic_one_sens_fabia = varphi_content$per_measures_omic_one$per_measures_true$Sensitivity[[1]],
      omic_one_sens_mfa   = varphi_content$per_measures_omic_one$per_measures_true$Sensitivity[[3]],
      omic_one_sens_gfa   = varphi_content$per_measures_omic_one$per_measures_true$Sensitivity[[4]],
      
      omic_two_sens_mofa = varphi_content$per_measures_omic_two$per_measures_true$Sensitivity[[2]],
      omic_two_sens_fabia = varphi_content$per_measures_omic_two$per_measures_true$Sensitivity[[1]],
      omic_two_sens_mfa   = varphi_content$per_measures_omic_two$per_measures_true$Sensitivity[[3]],
      omic_two_sens_gfa   = varphi_content$per_measures_omic_two$per_measures_true$Sensitivity[[4]],
      
      smp_sens_mofa = varphi_content$per_measures_smp$per_measures_true_smp$Sensitivity[[2]],
      smp_sens_fabia = varphi_content$per_measures_smp$per_measures_true_smp$Sensitivity[[1]],
      smp_sens_mfa   = varphi_content$per_measures_smp$per_measures_true_smp$Sensitivity[[3]],
      smp_sens_gfa   = varphi_content$per_measures_smp$per_measures_true_smp$Sensitivity[[4]],
      
      # Specificities
      omic_one_spec_mofa = varphi_content$per_measures_omic_one$per_measures_true$Specificity[[2]],
      omic_one_spec_fabia = varphi_content$per_measures_omic_one$per_measures_true$Specificity[[1]],
      omic_one_spec_mfa   = varphi_content$per_measures_omic_one$per_measures_true$Specificity[[3]],
      omic_one_spec_gfa   = varphi_content$per_measures_omic_one$per_measures_true$Specificity[[4]],
      
      omic_two_spec_mofa = varphi_content$per_measures_omic_two$per_measures_true$Specificity[[2]],
      omic_two_spec_fabia = varphi_content$per_measures_omic_two$per_measures_true$Specificity[[1]],
      omic_two_spec_mfa   = varphi_content$per_measures_omic_two$per_measures_true$Specificity[[3]],
      omic_two_spec_gfa   = varphi_content$per_measures_omic_two$per_measures_true$Specificity[[4]],
      
      smp_spec_mofa = varphi_content$per_measures_smp$per_measures_true_smp$Specificity[[2]],
      smp_spec_fabia = varphi_content$per_measures_smp$per_measures_true_smp$Specificity[[1]],
      smp_spec_mfa   = varphi_content$per_measures_smp$per_measures_true_smp$Specificity[[3]],
      smp_spec_gfa   = varphi_content$per_measures_smp$per_measures_true_smp$Specificity[[4]],
      
      # Accuracy
      omic_one_acc_mofa = varphi_content$per_measures_omic_one$per_measures_true$Accuracy[[2]],
      omic_one_acc_fabia = varphi_content$per_measures_omic_one$per_measures_true$Accuracy[[1]],
      omic_one_acc_mfa   = varphi_content$per_measures_omic_one$per_measures_true$Accuracy[[3]],
      omic_one_acc_gfa   = varphi_content$per_measures_omic_one$per_measures_true$Accuracy[[4]],
      
      omic_two_acc_mofa = varphi_content$per_measures_omic_two$per_measures_true$Accuracy[[2]],
      omic_two_acc_fabia = varphi_content$per_measures_omic_two$per_measures_true$Accuracy[[1]],
      omic_two_acc_mfa   = varphi_content$per_measures_omic_two$per_measures_true$Accuracy[[3]],
      omic_two_acc_gfa   = varphi_content$per_measures_omic_two$per_measures_true$Accuracy[[4]],
      
      smp_acc_mofa = varphi_content$per_measures_smp$per_measures_true_smp$Accuracy[[2]],
      smp_acc_fabia = varphi_content$per_measures_smp$per_measures_true_smp$Accuracy[[1]],
      smp_acc_mfa   = varphi_content$per_measures_smp$per_measures_true_smp$Accuracy[[3]],
      smp_acc_gfa   = varphi_content$per_measures_smp$per_measures_true_smp$Accuracy[[4]],
      
      # AUC
      omic_one_auc_mofa = varphi_content$per_measures_omic_one$per_measures_true$AUC[[2]],
      omic_one_auc_fabia = varphi_content$per_measures_omic_one$per_measures_true$AUC[[1]],
      omic_one_auc_mfa   = varphi_content$per_measures_omic_one$per_measures_true$AUC[[3]],
      omic_one_auc_gfa   = varphi_content$per_measures_omic_one$per_measures_true$AUC[[4]],
      
      omic_two_auc_mofa = varphi_content$per_measures_omic_two$per_measures_true$AUC[[2]],
      omic_two_auc_fabia = varphi_content$per_measures_omic_two$per_measures_true$AUC[[1]],
      omic_two_auc_mfa   = varphi_content$per_measures_omic_two$per_measures_true$AUC[[3]],
      omic_two_auc_gfa   = varphi_content$per_measures_omic_two$per_measures_true$AUC[[4]],
      
      smp_auc_mofa = varphi_content$per_measures_smp$per_measures_true_smp$AUC[[2]],
      smp_auc_fabia = varphi_content$per_measures_smp$per_measures_true_smp$AUC[[1]],
      smp_auc_mfa   = varphi_content$per_measures_smp$per_measures_true_smp$AUC[[3]],
      smp_auc_gfa   = varphi_content$per_measures_smp$per_measures_true_smp$AUC[[4]]
    )
    
    # Store in varphi list
    if (is.null(varphi_metric_data[[varphi_idx]])) {
      varphi_metric_data[[varphi_idx]] <- df
    } else {
      varphi_metric_data[[varphi_idx]] <- dplyr::bind_rows(varphi_metric_data[[varphi_idx]], df)
    }
  }
}

varphi_one_pm <- varphi_metric_data[[1]]
varphi_two_pm <- varphi_metric_data[[2]]
varphi_three_pm <- varphi_metric_data[[3]]
varphi_four_pm <- varphi_metric_data[[4]]
varphi_five_pm <- varphi_metric_data[[5]]
varphi_six_pm <- varphi_metric_data[[6]]
varphi_seven_pm <- varphi_metric_data[[7]]
varphi_eight_pm <- varphi_metric_data[[8]]
varphi_nine_pm <- varphi_metric_data[[9]]
varphi_ten_pm <- varphi_metric_data[[10]]
varphi_eleven_pm <- varphi_metric_data[[11]]
varphi_twelve_pm <- varphi_metric_data[[12]]

# PLOTTING THE DATA
generate_plot <- function(data, variable1, variable2, variable3, variable4, varphix){
  library(dplyr); library(tidyr); library(ggplot2)
  
  output_folder <- #"C:/Users/Lenovo/Downloads/results/single_factor/output"
  
  # Convert X to numeric or ordered factor for correct sorting
  data <- data %>% mutate(X = as.numeric(as.character(X)))
  
  # Apply the custom function to calculate confidence intervals for each variable
  calculate_ci <- function(x) {
    mean_val <- mean(x)
    ci <- 1.96 * sd(x) / sqrt(length(x))  # You can adjust the multiplier (1.96) for a different confidence level
    return(c(y = mean_val, ymin = mean_val - ci, ymax = mean_val + ci))
  }
  
  # Calculate confidence intervals for each variable
  ci_data_variable1 <- data %>%
    group_by(X) %>%
    summarize(conf_interval = list(calculate_ci(!!rlang::sym(variable1))))
  
  ci_data_variable2 <- data %>%
    group_by(X) %>%
    summarize(conf_interval = list(calculate_ci(!!rlang::sym(variable2))))
  
  ci_data_variable3 <- data %>%
    group_by(X) %>%
    summarize(conf_interval = list(calculate_ci(!!rlang::sym(variable3))))
  
  ci_data_variable4 <- data %>%
    group_by(X) %>%
    summarize(conf_interval = list(calculate_ci(!!rlang::sym(variable4))))
  
  # Unnest the list columns
  ci_data_variable1 <- unnest(ci_data_variable1, cols = conf_interval)
  ci_data_variable2 <- unnest(ci_data_variable2, cols = conf_interval)
  ci_data_variable3 <- unnest(ci_data_variable3, cols = conf_interval)
  ci_data_variable4 <- unnest(ci_data_variable4, cols = conf_interval)
  
  ci_data_variable1$type <- rep(c("mean", "lcl", "ucl"), times = nrow(ci_data_variable1) / 3)
  # Convert long to wide format with mean, lcl, and ucl as separate columns)
  wide_ci_data1 <- ci_data_variable1 %>% pivot_wider(names_from = type, values_from = conf_interval)
  wide_ci_data1 <- wide_ci_data1 %>%mutate(ucl = ifelse(ucl>1,1,ucl))
  
  ci_data_variable2$type <- rep(c("mean", "lcl", "ucl"), times = nrow(ci_data_variable2) / 3)
  # Convert long to wide format with mean, lcl, and ucl as separate columns)
  wide_ci_data2 <- ci_data_variable2 %>% pivot_wider(names_from = type, values_from = conf_interval)
  wide_ci_data2 <- wide_ci_data2 %>%mutate(ucl = ifelse(ucl>1,1,ucl))
  
  ci_data_variable3$type <- rep(c("mean", "lcl", "ucl"), times = nrow(ci_data_variable3) / 3)
  # Convert long to wide format with mean, lcl, and ucl as separate columns)
  wide_ci_data3 <- ci_data_variable3 %>% pivot_wider(names_from = type, values_from = conf_interval)
  wide_ci_data3 <- wide_ci_data3 %>%mutate(ucl = ifelse(ucl>1,1,ucl))
  
  ci_data_variable4$type <- rep(c("mean", "lcl", "ucl"), times = nrow(ci_data_variable4) / 3)
  # Convert long to wide format with mean, lcl, and ucl as separate columns)
  wide_ci_data4 <- ci_data_variable4 %>% pivot_wider(names_from = type, values_from = conf_interval)
  wide_ci_data4 <- wide_ci_data4 %>%mutate(ucl = ifelse(ucl>1,1,ucl))
  
  # Plot using ggplot
  plot <- ggplot() +
    geom_line(data = data, aes(x = X, y = !!rlang::sym(variable1), group = iteration, color = as.factor(iteration)), show.legend = TRUE) +
    geom_line(data = data, aes(x = X, y = !!rlang::sym(variable2), group = iteration, color = as.factor(iteration)), show.legend = TRUE) +
    geom_line(data = data, aes(x = X, y = !!rlang::sym(variable3), group = iteration, color = as.factor(iteration)), show.legend = TRUE) +
    geom_line(data = data, aes(x = X, y = !!rlang::sym(variable4), group = iteration, color = as.factor(iteration)), show.legend = TRUE) +
    
    geom_line(data = wide_ci_data1, aes(x = X, y = mean), color = "black", size = 0.8) +
    geom_line(data = wide_ci_data2, aes(x = X, y = mean), color = "blue", size = 0.8) +
    geom_line(data = wide_ci_data3, aes(x = X, y = mean), color = "red", size = 0.8) +
    geom_line(data = wide_ci_data4, aes(x = X, y = mean), color = "green", size = 0.8) +
    
    geom_errorbar(data = wide_ci_data1, aes(x = X, ymin = lcl, ymax = ucl), color = "black", width = 0.1, linetype = "solid") +
    geom_errorbar(data = wide_ci_data2, aes(x = X, ymin = lcl, ymax = ucl), color = "blue", width = 0.1, linetype = "solid") +
    geom_errorbar(data = wide_ci_data3, aes(x = X, ymin = lcl, ymax = ucl), color = "red", width = 0.1, linetype = "solid") +
    geom_errorbar(data = wide_ci_data4, aes(x = X, ymin = lcl, ymax = ucl), color = "green", width = 0.1, linetype = "solid") +
    
    # Add title with subscript for cutoff
    labs(
      title = bquote("Cutoff (" * varphi[.(varphix)] * ")"),#expression("Cutoff (" * varphi[varphix] * ")"),
      x = "Sigma (Noise)",
      y = "Jaccard Index"
    ) +
    
    # Adding a legend for Methods
    scale_color_manual(
      name = "Methods",
      values = c("MOFA" = "black", "FABIA" = "blue", "MFA" = "red", "GFA" = "green"),
      labels = c("MOFA", "FABIA", "GFA", "MFA")
    ) +
    
    theme_bw() +
    theme(
      legend.position = "top",  # Place legend at the top for clarity
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(color = "black")
    ) +
    scale_y_continuous(limits = c(0.0, 1.0))
  
  # Save the plot
  plot_name <- paste(variable1, variable2, variable3,variable4, sep = "_") 
  ggsave(file.path(output_folder, paste0(plot_name, "_plot.png")), plot, width = 8, height = 8, dpi = 300)
  
  return(plot)
}

# TRUE JI
generate_plot(varphi_one, "ji_true_mofa", "ji_true_fabia", "ji_true_mfa", "ji_true_gfa", varphix = 1)
generate_plot(varphi_two, "ji_true_mofa", "ji_true_fabia", "ji_true_mfa", "ji_true_gfa", varphix = 2)
generate_plot(varphi_three, "ji_true_mofa", "ji_true_fabia", "ji_true_mfa", "ji_true_gfa", varphix = 3)
generate_plot(varphi_four, "ji_true_mofa", "ji_true_fabia", "ji_true_mfa", "ji_true_gfa", varphix = 4)
generate_plot(varphi_five, "ji_true_mofa", "ji_true_fabia", "ji_true_mfa", "ji_true_gfa", varphix = 5)
generate_plot(varphi_six, "ji_true_mofa", "ji_true_fabia", "ji_true_mfa", "ji_true_gfa", varphix = 6)
generate_plot(varphi_seven, "ji_true_mofa", "ji_true_fabia", "ji_true_mfa", "ji_true_gfa", varphix = 7)
generate_plot(varphi_nine, "ji_true_mofa", "ji_true_fabia", "ji_true_mfa", "ji_true_gfa", varphix = 8)
generate_plot(varphi_eight, "ji_true_mofa", "ji_true_fabia", "ji_true_mfa", "ji_true_gfa", varphix = 9)
generate_plot(varphi_ten, "ji_true_mofa", "ji_true_fabia", "ji_true_mfa", "ji_true_gfa", varphix = 10)
generate_plot(varphi_eleven, "ji_true_mofa", "ji_true_fabia", "ji_true_mfa", "ji_true_gfa", varphix = 11)
generate_plot(varphi_twelve, "ji_true_mofa", "ji_true_fabia", "ji_true_mfa", "ji_true_gfa", varphix = 12)
generate_plot(varphi_thirteen, "ji_true_mofa", "ji_true_fabia", "ji_true_mfa", "ji_true_gfa", varphix = 13)

# OMIC 1 JI
generate_plot(varphi_one, "ji_omic_one_mofa", "ji_omic_one_fabia", "ji_omic_one_mfa", "ji_omic_one_gfa", varphix = 1)
generate_plot(varphi_two, "ji_omic_one_mofa", "ji_omic_one_fabia", "ji_omic_one_mfa", "ji_omic_one_gfa", varphix = 2)
generate_plot(varphi_three, "ji_omic_one_mofa", "ji_omic_one_fabia", "ji_omic_one_mfa", "ji_omic_one_gfa", varphix = 3)
generate_plot(varphi_four, "ji_omic_one_mofa", "ji_omic_one_fabia", "ji_omic_one_mfa", "ji_omic_one_gfa", varphix = 4)
generate_plot(varphi_five, "ji_omic_one_mofa", "ji_omic_one_fabia", "ji_omic_one_mfa", "ji_omic_one_gfa", varphix = 5)
generate_plot(varphi_six, "ji_omic_one_mofa", "ji_omic_one_fabia", "ji_omic_one_mfa", "ji_omic_one_gfa", varphix = 6)
generate_plot(varphi_seven, "ji_omic_one_mofa", "ji_omic_one_fabia", "ji_omic_one_mfa", "ji_omic_one_gfa", varphix = 7)
generate_plot(varphi_eight, "ji_omic_one_mofa", "ji_omic_one_fabia", "ji_omic_one_mfa", "ji_omic_one_gfa", varphix = 8)
generate_plot(varphi_nine, "ji_omic_one_mofa", "ji_omic_one_fabia", "ji_omic_one_mfa", "ji_omic_one_gfa", varphix = 9)
generate_plot(varphi_ten, "ji_omic_one_mofa", "ji_omic_one_fabia", "ji_omic_one_mfa", "ji_omic_one_gfa", varphix = 10)
generate_plot(varphi_eleven, "ji_omic_one_mofa", "ji_omic_one_fabia", "ji_omic_one_mfa", "ji_omic_one_gfa", varphix = 11)
generate_plot(varphi_twelve, "ji_omic_one_mofa", "ji_omic_one_fabia", "ji_omic_one_mfa", "ji_omic_one_gfa", varphix = 12)
generate_plot(varphi_thirteen, "ji_omic_one_mofa", "ji_omic_one_fabia", "ji_omic_one_mfa", "ji_omic_one_gfa", varphix = 13)

# OMIC 2 JI
generate_plot(varphi_one, "ji_omic_two_mofa", "ji_omic_two_fabia", "ji_omic_two_mfa", "ji_omic_two_gfa", varphix = 1)
generate_plot(varphi_two,  "ji_omic_two_mofa", "ji_omic_two_fabia", "ji_omic_two_mfa", "ji_omic_two_gfa", varphix = 2)
generate_plot(varphi_three,  "ji_omic_two_mofa", "ji_omic_two_fabia", "ji_omic_two_mfa", "ji_omic_two_gfa", varphix = 3)
generate_plot(varphi_four,  "ji_omic_two_mofa", "ji_omic_two_fabia", "ji_omic_two_mfa", "ji_omic_two_gfa", varphix = 4)
generate_plot(varphi_five,  "ji_omic_two_mofa", "ji_omic_two_fabia", "ji_omic_two_mfa", "ji_omic_two_gfa", varphix = 5)
generate_plot(varphi_six,  "ji_omic_two_mofa", "ji_omic_two_fabia", "ji_omic_two_mfa", "ji_omic_two_gfa", varphix = 6)
generate_plot(varphi_seven,  "ji_omic_two_mofa", "ji_omic_two_fabia", "ji_omic_two_mfa", "ji_omic_two_gfa", varphix = 7)
generate_plot(varphi_eight,  "ji_omic_two_mofa", "ji_omic_two_fabia", "ji_omic_two_mfa", "ji_omic_two_gfa", varphix = 8)
generate_plot(varphi_nine,  "ji_omic_two_mofa", "ji_omic_two_fabia", "ji_omic_two_mfa", "ji_omic_two_gfa", varphix = 9)
generate_plot(varphi_ten,  "ji_omic_two_mofa", "ji_omic_two_fabia", "ji_omic_two_mfa", "ji_omic_two_gfa", varphix = 10)
generate_plot(varphi_eleven,  "ji_omic_two_mofa", "ji_omic_two_fabia", "ji_omic_two_mfa", "ji_omic_two_gfa", varphix = 11)
generate_plot(varphi_twelve,  "ji_omic_two_mofa", "ji_omic_two_fabia", "ji_omic_two_mfa", "ji_omic_two_gfa", varphix = 12)
generate_plot(varphi_thirteen,  "ji_omic_two_mofa", "ji_omic_two_fabia", "ji_omic_two_mfa", "ji_omic_two_gfa", varphix = 13)

# MATRIX JI
generate_plot(varphi_one, "ji_matrix_mofa", "ji_matrix_fabia", "ji_matrix_mfa", "ji_matrix_gfa", varphix = 1)
generate_plot(varphi_two, "ji_matrix_mofa", "ji_matrix_fabia", "ji_matrix_mfa", "ji_matrix_gfa", varphix = 2)
generate_plot(varphi_three, "ji_matrix_mofa", "ji_matrix_fabia", "ji_matrix_mfa", "ji_matrix_gfa", varphix = 3)
generate_plot(varphi_four, "ji_matrix_mofa", "ji_matrix_fabia", "ji_matrix_mfa", "ji_matrix_gfa", varphix = 4)
generate_plot(varphi_five, "ji_matrix_mofa", "ji_matrix_fabia", "ji_matrix_mfa", "ji_matrix_gfa", varphix = 5)
generate_plot(varphi_six, "ji_matrix_mofa", "ji_matrix_fabia", "ji_matrix_mfa", "ji_matrix_gfa", varphix = 6)
generate_plot(varphi_seven, "ji_matrix_mofa", "ji_matrix_fabia", "ji_matrix_mfa", "ji_matrix_gfa", varphix = 7)
generate_plot(varphi_eight, "ji_matrix_mofa", "ji_matrix_fabia", "ji_matrix_mfa", "ji_matrix_gfa", varphix = 8)
generate_plot(varphi_nine, "ji_matrix_mofa", "ji_matrix_fabia", "ji_matrix_mfa", "ji_matrix_gfa", varphix = 9)
generate_plot(varphi_ten, "ji_matrix_mofa", "ji_matrix_fabia", "ji_matrix_mfa", "ji_matrix_gfa", varphix = 10)
generate_plot(varphi_eleven, "ji_matrix_mofa", "ji_matrix_fabia", "ji_matrix_mfa", "ji_matrix_gfa", varphix = 11)
generate_plot(varphi_twelve, "ji_matrix_mofa", "ji_matrix_fabia", "ji_matrix_mfa", "ji_matrix_gfa", varphix = 12)
generate_plot(varphi_thirteen, "ji_matrix_mofa", "ji_matrix_fabia", "ji_matrix_mfa", "ji_matrix_gfa", varphix = 13)
