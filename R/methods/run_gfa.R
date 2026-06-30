# R/methods/run_gfa.R
# Group Factor Analysis (GFA) wrapper for multi-omics data
# Klami et al. (2014); implemented via the GFA package
#
# GFA is a Bayesian group factor model with automatic relevance determination
# (ARD) at the dataset (group) level, so factors can be active in a subset of views.
#
# Returns the common interface:
#   list(scores = k x n, loadings = named list of (features_m x k), model = raw fit)

#' Run GFA on a multi-omics dataset
#'
#' @param omics_list Named list of matrices, each (features x samples).
#'                   All matrices must share the same columns (samples).
#' @param n_factors  Integer. Number of latent factors (upper bound; ARD prunes).
#' @param iter_max   Integer. Maximum Gibbs iterations. Default 1000.
#' @param iter_burn  Integer. Burn-in iterations. Default 100.
#' @param seed       Integer. Random seed. Default 60667.
#' @return List with elements: scores (k x n), loadings (named list), model (GFA object)
run_gfa_multiomics <- function(omics_list,
                               n_factors = 5,
                               iter_max  = 1000,
                               iter_burn = 100,
                               seed      = 60667) {
  if (!requireNamespace("GFA", quietly = TRUE)) {
    stop("Package 'GFA' is required. install.packages('GFA').")
  }
  set.seed(seed)

  n_cols <- vapply(omics_list, ncol, integer(1))
  if (length(unique(n_cols)) != 1) {
    stop("All omics matrices must have the same number of samples (columns).")
  }

  # GFA expects a list of (samples x features) views sharing the same rows (samples)
  views <- lapply(omics_list, function(m) t(as.matrix(m)))
  views <- GFA::normalizeData(views, type = "scaleFeatures")

  opts <- GFA::getDefaultOpts()
  opts$iter.max    <- iter_max
  opts$iter.burnin <- iter_burn

  fit <- GFA::gfa(views$train, K = n_factors, opts = opts)

  # Scores: latent variables X (samples x k) -> (k x samples)
  scores <- t(fit$X[, seq_len(min(n_factors, ncol(fit$X))), drop = FALSE])

  # Loadings: projection matrix W (features x k), split into per-view blocks by group
  group_sizes <- vapply(views$train, ncol, integer(1))
  split_at    <- cumsum(c(0, group_sizes))
  loadings <- lapply(seq_along(omics_list), function(m) {
    rows <- (split_at[m] + 1):split_at[m + 1]
    fit$W[rows, seq_len(min(n_factors, ncol(fit$W))), drop = FALSE]
  })
  names(loadings) <- names(omics_list)

  list(scores = scores, loadings = loadings, model = fit)
}
