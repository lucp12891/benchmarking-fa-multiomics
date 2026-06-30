# R/utils/libs.R
# Centralized package loading. Sourced at the top of analysis scripts.

.required_packages <- c(
  # Factor analysis methods
  "fabia", "FactoMineR", "GFA", "MOFA2",
  # Data wrangling
  "dplyr", "tidyr", "tibble", "readr", "readxl", "stringr", "data.table",
  # Modelling / metrics utilities
  "pROC", "mclust", "zoo",
  # Plotting
  "ggplot2", "ggupset", "VennDiagram", "pheatmap",
  # Simulation
  "SUMO"
)

#' Load all required packages, quietly
#'
#' @param packages Character vector of package names (defaults to the full list).
#' @return Invisibly, a logical vector of which packages loaded successfully.
load_libraries <- function(packages = .required_packages) {
  ok <- vapply(packages, function(pkg) {
    suppressPackageStartupMessages(
      suppressWarnings(
        requireNamespace(pkg, quietly = TRUE) &&
          library(pkg, character.only = TRUE, logical.return = TRUE)
      )
    )
  }, logical(1))

  if (any(!ok)) {
    warning("The following packages could not be loaded: ",
            paste(names(ok)[!ok], collapse = ", "),
            "\nSee README.md for installation instructions.")
  }
  invisible(ok)
}
