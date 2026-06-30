# -------------------------- Load Required Libraries --------------------------
libs <- c(
  "tidyr", "fabia", "FactoMineR", "readr", "readxl", 
  "tidyverse", "dplyr", "stringr", "basilisk", "MOFA2", 
  "data.table", "GFA", "pROC", "mclust", "zoo"
)

# Load each library, suppress messages and warnings
invisible(lapply(libs, function(pkg) {
  suppressPackageStartupMessages(
    suppressWarnings(
      library(pkg, character.only = TRUE)
    )
  )
}))

