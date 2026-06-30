# R/methods/run_fabia.R
# FABIA adapted for multi-omics via concatenation (Equation 8 in paper)
# Hochreiter et al. (2010); Kasim et al. (2016)

#' Run FABIA on multi-omics data (concatenation strategy)
#'
#' Concatenates omics matrices row-wise and applies FABIA biclustering.
#' Recovers dataset-specific loading blocks post-hoc.
#'
#' @param omics_list  Named list of matrices, each (features × samples).
#'                    All matrices must share the same columns (samples).
#' @param n_factors   Integer. Number of latent factors (biclusters) to extract.
#' @param alpha       Numeric. FABIA regularisation parameter (sparsity strength). Default 0.1.
#' @param spz         Numeric. Sparsity of factor scores. Default 0.5.
#' @param spl         Numeric. Sparsity of factor loadings. Default 0.5.
#' @param max_iter    Integer. Maximum EM iterations. Default 500.
#' @param ...         Additional arguments passed to fabia::fabia()
#'
#' @return List with:
#'   \item{scores}{Matrix (k × samples) of factor scores (shared across omics)}
#'   \item{loadings}{Named list of loading matrices, one per omics layer (features_m × k)}
#'   \item{model}{Raw fabia model object}
run_fabia_multiomics <- function(omics_list,
                                  n_factors = 5,
                                  alpha     = 0.1,
                                  spz       = 0.5,
                                  spl       = 0.5,
                                  max_iter  = 500,
                                  ...) {
  if (!requireNamespace("fabia", quietly = TRUE)) {
    stop("Package 'fabia' is required. Install via BiocManager::install('fabia').")
  }

  # Validate: all matrices must share same number of columns
  n_cols <- sapply(omics_list, ncol)
  if (length(unique(n_cols)) != 1) {
    stop("All omics matrices must have the same number of samples (columns).")
  }

  # Concatenate row-wise: (sum(p_m) × n)
  X_concat <- do.call(rbind, omics_list)

  # Run FABIA
  message("Running FABIA with ", n_factors, " factors ...")
  model <- fabia::fabia(
    X       = X_concat,
    p       = n_factors,
    alpha   = alpha,
    cyc     = max_iter,
    spz     = spz,
    spl     = spl,
    ...
  )

  # Extract shared factor scores (k × n)
  scores <- fabia::extractBic(model)$bic

  # If extractBic fails, fall back to raw Z slot
  if (is.null(scores) || nrow(scores) == 0) {
    scores <- model@Z
  }

  # Extract loading matrix for full concatenated data ((sum p_m) × k)
  W_full <- model@L  # features × factors

  # Split loadings back into per-omics blocks
  n_features <- sapply(omics_list, nrow)
  split_idx  <- cumsum(c(0, n_features))

  loadings <- lapply(seq_along(omics_list), function(m) {
    rows <- (split_idx[m] + 1):split_idx[m + 1]
    W_full[rows, , drop = FALSE]
  })
  names(loadings) <- names(omics_list)

  list(
    scores   = scores,
    loadings = loadings,
    model    = model
  )
}

#' Explore FABIA sparsity parameter grid
#'
#' Reproduces Figure 8: proportion of non-zero loadings/scores across (alpha, spz, spl) grid.
#'
#' @param omics_list  Named list of omics matrices
#' @param n_factors   Number of factors
#' @param alphas      Vector of alpha values
#' @param spz_vals    Vector of spz values
#' @param spl_vals    Vector of spl values
#' @return Data frame with columns: alpha, spz, spl, prop_nonzero_scores,
#'         and one prop_nonzero_loadings_<omic> column per omics layer
run_fabia_sparsity_grid <- function(omics_list,
                                     n_factors = 5,
                                     alphas    = seq(0.01, 0.10, by = 0.01),
                                     spz_vals  = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                                     spl_vals  = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5)) {
  grid <- expand.grid(alpha = alphas, spz = spz_vals, spl = spl_vals,
                      stringsAsFactors = FALSE)

  results <- lapply(seq_len(nrow(grid)), function(i) {
    params <- grid[i, ]
    message(sprintf("Grid [%d/%d]: alpha=%.2f, spz=%.2f, spl=%.2f",
                    i, nrow(grid), params$alpha, params$spz, params$spl))

    fit <- tryCatch(
      run_fabia_multiomics(
        omics_list, n_factors = n_factors,
        alpha = params$alpha, spz = params$spz, spl = params$spl
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)

    # Proportion of non-zero scores (averaged across factors)
    prop_scores <- mean(fit$scores != 0)

    # Proportion of non-zero loadings per omics layer
    prop_loadings <- sapply(fit$loadings, function(L) mean(L != 0))

    row <- data.frame(params, prop_nonzero_scores = prop_scores,
                      stringsAsFactors = FALSE)
    for (nm in names(prop_loadings)) {
      row[[paste0("prop_nonzero_", nm)]] <- prop_loadings[[nm]]
    }
    row
  })

  do.call(rbind, Filter(Negate(is.null), results))
}
