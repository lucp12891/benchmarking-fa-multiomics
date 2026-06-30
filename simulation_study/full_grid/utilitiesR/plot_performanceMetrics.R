generate_perMetric_plot <- function(data, variables, varphi, metric, path) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  
  output_folder <- path
  
  data <- data %>%
    mutate(X = as.numeric(as.character(X)))
  
  long_data <- data %>%
    pivot_longer(
      cols = all_of(variables), 
      names_to = "method", 
      values_to = "value"
    ) %>%
    mutate(
      method = case_when(
        str_detect(method, "fabia") ~ "FABIA",
        str_detect(method, "gfa") ~ "GFA",
        str_detect(method, "mfa") ~ "MFA",
        str_detect(method, "mofa") ~ "MOFA",
        TRUE ~ method
      )
    )
  
  # Function to calculate CI
  calculate_ci <- function(x) {
    mean_val <- mean(x, na.rm = TRUE)
    ci <- 1.96 * sd(x, na.rm = TRUE) / sqrt(length(x))
    tibble(mean = mean_val, lcl = mean_val - ci, ucl = mean_val + ci)
  }
  
  ci_data <- long_data %>%
    group_by(X, method) %>%
    summarize(calculate_ci(value), .groups = "drop")
  
  mean_colors <- c(
    "FABIA" = "blue", "GFA" = "gray20", "MFA" = "limegreen", "MOFA" = "#990000"
  )
  
  plot <- ggplot(ci_data, aes(x = X, y = mean, color = method)) +
    geom_line(linewidth = 1.2) +
    geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.1) +
    labs(
      title = bquote("Cutoff (" * varphi[.(varphi)] * ")"),
      x = "Sigma (Noise)",
      y = metric,
      color = "Methods"
    ) +
    scale_color_manual(values = mean_colors) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 14),  # Increase X-axis tick size
      axis.text.y = element_text(size = 14),  # Increase Y-axis tick size
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    ) +
    scale_y_continuous(limits = c(0.0, 1.0)) +
    scale_x_continuous(breaks = seq(0, 30, by = 5))
  
  # Save
  plot_name <- paste0("plot_varphi_", varphi)
  ggsave(file.path(output_folder, paste0(plot_name, ".png")), plot, width = 5.00, height = 3.00, dpi = 300, bg = "white")
  ggsave(file.path(output_folder, paste0(plot_name, ".eps")), plot, width = 6.56, height = 4.5, device = cairo_ps, dpi = 300, bg = "white")
  
  return(plot)
}

# Define the datasets
#root_dir = "C:/Users/bosangir/Downloads/test_sim/Run_24_05_2025/plots_results/single_plots"
root_dir = "C:/Users/bosangir/Downloads/test_sim/Run_24_05_2025/plots_results/single_plots"
datasets <- varphi_metric_data[c(1, 3, 5, 6, 9, 12)]


# Generate all plots
generate_all_pm_plots <- function(datasets) {
  
  # 1) Define a root directory for the output
  root_dir <- root_dir
  
  # 2) Map each metric to its respective subfolder name
  subdir_map <- c(
    "Sensitivity"              = "sensitivity",
    "Specificity"              = "specificity",
    "Accuracy"                 = "accuracy",
    "Area Under the Curve (AUC)" = "auc"
  )
  
  # 3) For convenience, define a nested list of variables for each combination
  #    of data_group ("omic1", "omic2", "smp") and metric.
  #    Adjust names as needed if yours differ.
  var_map <- list(
    omic1 = list(
      "Sensitivity" = c("omic_one_sens_fabia", "omic_one_sens_mofa", 
                        "omic_one_sens_mfa",   "omic_one_sens_gfa"),
      "Specificity" = c("omic_one_spec_fabia", "omic_one_sens_mofa", 
                        "omic_one_sens_mfa",   "omic_one_sens_gfa"),
      "Accuracy"    = c("omic_one_acc_fabia", "omic_one_acc_mofa", 
                        "omic_one_acc_mfa",   "omic_one_acc_gfa"),
      "Area Under the Curve (AUC)" = c("omic_one_auc_fabia", "omic_one_auc_mofa",
                                       "omic_one_auc_mfa",   "omic_one_auc_gfa")
    ),
    omic2 = list(
      "Sensitivity" = c("omic_two_sens_fabia", "omic_two_sens_mofa", 
                        "omic_two_sens_mfa",   "omic_two_sens_gfa"),
      "Specificity" = c("omic_two_spec_fabia", "omic_two_spec_mofa",
                        "omic_two_spec_mfa",   "omic_two_spec_gfa"),
      "Accuracy"    = c("omic_two_acc_fabia", "omic_two_acc_mofa",
                        "omic_two_acc_mfa",   "omic_two_acc_gfa"),
      "Area Under the Curve (AUC)" = c("omic_two_auc_fabia", "omic_two_auc_mofa",
                                       "omic_two_auc_mfa",   "omic_two_auc_gfa")
    ),
    smp = list(
      "Sensitivity" = c("smp_sens_fabia", "smp_sens_mofa",
                        "smp_sens_mfa",   "smp_sens_gfa"),
      "Specificity" = c("smp_spec_fabia", "smp_spec_mofa",
                        "smp_spec_mfa",   "smp_spec_gfa"),
      "Accuracy"    = c("smp_acc_fabia", "smp_acc_mofa", 
                        "smp_acc_mfa",   "smp_acc_gfa"),
      "Area Under the Curve (AUC)" = c("smp_auc_fabia", "smp_auc_gfa",
                                       "smp_auc_mfa",   "smp_auc_mofa")
    )
  )
  
  # 4) Define the data groups and metrics you want to iterate over
  data_groups <- c("omic1", "omic2", "smp")
  metrics     <- c("Sensitivity", "Specificity", "Accuracy", "Area Under the Curve (AUC)")
  
  # 5) Loop over each data group, each metric, and each dataset
  for (dg in data_groups) {
    for (m in metrics) {
      
      # Build the path for saving
      metric_subdir <- subdir_map[m]  # e.g., "sensitivity", "specificity", etc.
      save_path     <- file.path(root_dir, dg, metric_subdir)
      
      # Get the appropriate variables for this combination
      vars <- var_map[[dg]][[m]]
      
      # Loop over your dataset list
      for (i in seq_along(datasets)) {
        
        # Call your plotting function with the needed parameters
        generate_perMetric_plot(
          data   = datasets[[i]],
          variables   = vars,
          varphi = i,
          metric = m,
          path   = save_path
        )
      }
    }
  }
  
  message("All plots have been generated and saved.")
}

# Run the function
generate_all_pm_plots(datasets)

