# R/methods/run_mfa.R
# Multiple Factor Analysis (MFA) wrapper for multi-omics data
# Abdi et al. (2007); implemented via FactoMineR::MFA
#
# MFA balances the influence of each omics block by dividing each block by its
# first singular value, then performs a global PCA on the concatenated, weighted data.
#
# Returns the common interface:
#   list(scores = k x n, loadings = named list of (features_m x k), model = raw fit)

#' Run MFA on a multi-omics dataset
#'
#' @param omics_list Named list of matrices, each (features x samples).
#'                   All matrices must share the same columns (samples).
#' @param n_factors  Integer. Number of dimensions to retain.
#' @param scale      Logical. Standardize features within each block. Default TRUE.
#' @param seed       Integer. Random seed. Default 60667.
#' @return List with elements: scores (k x n), loadings (named list), model (MFA object)
run_mfa_multiomics <- function(omics_list,
                               n_factors = 5,
                               scale     = TRUE,
                               seed      = 60667) {
  if (!requireNamespace("FactoMineR", quietly = TRUE)) {
    stop("Package 'FactoMineR' is required. install.packages('FactoMineR').")
  }
  set.seed(seed)

  n_cols <- vapply(omics_list, ncol, integer(1))
  if (length(unique(n_cols)) != 1) {
    stop("All omics matrices must have the same number of samples (columns).")
  }

  # FactoMineR works on (individuals x variables): samples in rows, features in columns.
  blocks      <- lapply(omics_list, function(m) t(as.matrix(m)))
  group_sizes <- vapply(blocks, ncol, integer(1))
  cdata       <- do.call(cbind, blocks)

  mfa_fit <- FactoMineR::MFA(
    as.data.frame(cdata),
    group      = as.integer(group_sizes),
    type       = rep(if (scale) "s" else "c", length(blocks)),
    name.group = names(omics_list),
    ncp        = n_factors,
    graph      = FALSE
  )

  # Scores: individuals (samples) coordinates -> (k x samples)
  scores <- t(mfa_fit$ind$coord[, seq_len(n_factors), drop = FALSE])

  # Loadings: quantitative variable coordinates, split back into per-omics blocks
  W_full   <- mfa_fit$quanti.var$coord[, seq_len(n_factors), drop = FALSE]
  split_at <- cumsum(c(0, group_sizes))
  loadings <- lapply(seq_along(omics_list), function(m) {
    W_full[(split_at[m] + 1):split_at[m + 1], , drop = FALSE]
  })
  names(loadings) <- names(omics_list)

  list(scores = scores, loadings = loadings, model = mfa_fit)
}
