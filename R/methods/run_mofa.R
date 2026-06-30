# R/methods/run_mofa.R
# MOFA / MOFA+ wrapper for multi-omics factor analysis
# Argelaguet et al. (2018); Argelaguet et al. (2020, MOFA+)
#
# Returns the common interface used across all four methods:
#   list(scores = k x n, loadings = named list of (features_m x k), model = raw fit)

#' Run MOFA+ on a multi-omics dataset
#'
#' @param omics_list  Named list of matrices, each (features x samples).
#'                    All matrices must share the same columns (samples).
#' @param n_factors   Integer. Number of latent factors to extract.
#' @param scale_views Logical. Let MOFA scale each view to unit variance. Default TRUE.
#' @param maxiter     Integer. Maximum training iterations. Default 1000.
#' @param convergence Character. MOFA convergence mode ("fast", "medium", "slow"). Default "fast".
#' @param seed        Integer. Random seed for reproducibility. Default 60667.
#' @return List with elements: scores (k x n), loadings (named list), model (MOFA object)
run_mofa_multiomics <- function(omics_list,
                                n_factors   = 5,
                                scale_views = TRUE,
                                maxiter     = 1000,
                                convergence = "fast",
                                seed        = 60667) {
  if (!requireNamespace("MOFA2", quietly = TRUE)) {
    stop("Package 'MOFA2' is required. Install via BiocManager::install('MOFA2').")
  }
  set.seed(seed)

  n_cols <- vapply(omics_list, ncol, integer(1))
  if (length(unique(n_cols)) != 1) {
    stop("All omics matrices must have the same number of samples (columns).")
  }

  # MOFA expects a named list of (features x samples) matrices
  mofa_input <- lapply(omics_list, function(m) {
    M <- as.matrix(m); storage.mode(M) <- "double"
    if (is.null(rownames(M))) rownames(M) <- paste0("feat", seq_len(nrow(M)))
    M
  })

  mofa_obj   <- MOFA2::create_mofa(mofa_input)

  data_opts  <- MOFA2::get_default_data_options(mofa_obj)
  data_opts$scale_views  <- scale_views
  data_opts$scale_groups <- TRUE

  model_opts <- MOFA2::get_default_model_options(mofa_obj)
  model_opts$num_factors <- n_factors

  train_opts <- MOFA2::get_default_training_options(mofa_obj)
  train_opts$maxiter          <- maxiter
  train_opts$convergence_mode <- convergence
  train_opts$seed             <- seed

  mofa_obj <- MOFA2::prepare_mofa(mofa_obj, data_opts, model_opts, train_opts)
  outfile  <- file.path(tempdir(), "mofa_model.hdf5")
  fit <- suppressWarnings(MOFA2::run_mofa(mofa_obj, outfile, use_basilisk = TRUE))

  # Scores: MOFA returns (samples x k) -> transpose to (k x samples)
  Z <- MOFA2::get_factors(fit, factors = "all")[[1]]
  scores <- t(as.matrix(Z))

  # Loadings: one (features_m x k) matrix per view
  W <- MOFA2::get_weights(fit, views = "all", factors = "all")
  loadings <- lapply(W, as.matrix)
  names(loadings) <- names(omics_list)

  list(scores = scores, loadings = loadings, model = fit)
}
