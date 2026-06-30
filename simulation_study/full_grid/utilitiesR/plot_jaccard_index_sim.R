generate_jaccard_plot <- function(data, variables, varphi, path) {
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
  
  calculate_ci <- function(x) {
    mean_val <- mean(x, na.rm = TRUE)
    ci <- 1.96 * sd(x, na.rm = TRUE) / sqrt(length(x))
    tibble(mean = mean_val, lcl = mean_val - ci, ucl = mean_val + ci)
  }
  
  ci_data <- long_data %>%
    group_by(X, method) %>%
    summarize(calculate_ci(value), .groups = "drop")
  
  mean_colors <- c("FABIA" = "blue", "GFA" = "gray20", "MFA" = "limegreen", "MOFA" = "#990000")
  
  plot <- ggplot(ci_data, aes(x = X, y = mean, color = method)) +
    geom_line(linewidth = 1.2) +
    geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.1) +
    labs(
      title = bquote("Cutoff (" * varphi[.(varphi)] * ")"),
      x = "Sigma (Noise)",
      y = "Jaccard Index",
      color = "Methods"
    ) +
    scale_color_manual(values = mean_colors) +
    theme_minimal() +
    theme(
      legend.position = "top",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    ) +
    scale_y_continuous(limits = c(0.0, 1.0)) +
    scale_x_continuous(breaks = seq(0, 30, by = 5))
  
  # Save plots
  plot_name <- paste0("plot_varphi_", varphi)
  ggsave(file.path(output_folder, paste0(plot_name, ".png")), plot, width = 7.35, height = 5.04, dpi = 300, bg = "white")
  ggsave(file.path(output_folder, paste0(plot_name, ".eps")), plot, width = 6.56, height = 4.5, device = cairo_ps, dpi = 300, bg = "white")
  
  return(plot)
}


generate_jaccard_plot(
  data = datasets[[4]], 
  variables = jaccard_var_map$omic1, 
  varphi = 4, 
  path = file.path(root_dir, "omic1", "jaccard")
)

generate_jaccard_plot(
  data = datasets[[4]], 
  variables = jaccard_var_map$omic2, 
  varphi = 4, 
  path = file.path(root_dir, "omic2", "jaccard")
)

generate_jaccard_plot(
  data = datasets[[4]], 
  variables = jaccard_var_map$smp, 
  varphi = 4, 
  path = file.path(root_dir, "smp", "jaccard")
)

generate_jaccard_plot(
  data = datasets[[4]], 
  variables = jaccard_var_map$matrix, 
  varphi = 4, 
  path = file.path(root_dir, "matrix", "jaccard")
)
