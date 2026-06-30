# scripts/02_run_methods.R
# Run all four factor-analysis methods (MOFA, MFA, GFA, FABIA) on each prepared
# real-world dataset, and save the harmonized fits.
#
# Input : data/CLL/cll_omics.rds, data/radiation/radiation_omics.rds (from 01)
# Output: results/<dataset>_fits.rds -- named list per method, each
#         list(scores = k x n, loadings = named list (features_m x k), model = raw fit)

source("R/utils/libs.R")
source("R/methods/run_mofa.R")
source("R/methods/run_mfa.R")
source("R/methods/run_gfa.R")
source("R/methods/run_fabia.R")

dir.create("results", recursive = TRUE, showWarnings = FALSE)

#' Fit every method on one dataset
#'
#' @param omics_list Named list of (features x samples) matrices.
#' @param n_factors  Number of factors to extract.
#' @param methods    Subset of c("MOFA","MFA","GFA","FABIA").
#' @return Named list of fits (NULL entries for methods that errored).
run_all_methods <- function(omics_list, n_factors = 10,
                            methods = c("MOFA", "MFA", "GFA", "FABIA")) {
  fits <- list()
  if ("MOFA"  %in% methods)
    fits$MOFA  <- tryCatch(run_mofa_multiomics(omics_list,  n_factors), error = function(e) { message("MOFA: ",  conditionMessage(e)); NULL })
  if ("MFA"   %in% methods)
    fits$MFA   <- tryCatch(run_mfa_multiomics(omics_list,   n_factors), error = function(e) { message("MFA: ",   conditionMessage(e)); NULL })
  if ("GFA"   %in% methods)
    fits$GFA   <- tryCatch(run_gfa_multiomics(omics_list,   n_factors), error = function(e) { message("GFA: ",   conditionMessage(e)); NULL })
  if ("FABIA" %in% methods)
    fits$FABIA <- tryCatch(run_fabia_multiomics(omics_list, n_factors), error = function(e) { message("FABIA: ", conditionMessage(e)); NULL })
  Filter(Negate(is.null), fits)
}

run_dataset <- function(rds_in, rds_out, n_factors = 10) {
  if (!file.exists(rds_in)) {
    warning("Missing input: ", rds_in, " (run scripts/01_prepare_data.R first).")
    return(invisible(NULL))
  }
  omics <- readRDS(rds_in)
  message("\n=== Fitting methods on ", rds_in, " ===")
  fits <- run_all_methods(omics, n_factors = n_factors)
  saveRDS(fits, rds_out)
  message("Saved fits -> ", rds_out)
  invisible(fits)
}

if (!interactive()) {
  run_dataset("data/CLL/cll_omics.rds",             "results/cll_fits.rds",       n_factors = 10)
  run_dataset("data/radiation/radiation_omics.rds", "results/radiation_fits.rds", n_factors = 5)
}
