# scripts/03_compute_metrics.R
# Compute cross-method agreement metrics on the fitted real-data models.
#
# For each dataset and each (aligned) factor:
#   - pairwise Jaccard of selected features across a sweep of cutoffs
#   - pairwise cosine similarity and Pearson correlation of raw loadings
#
# Input : results/<dataset>_fits.rds (from 02)
# Output: results/<dataset>_metrics.rds

source("R/utils/libs.R")
source("R/utils/metrics.R")
source("R/utils/factor_alignment.R")

dir.create("results", recursive = TRUE, showWarnings = FALSE)

#' Compute agreement metrics for one fitted dataset
#'
#' @param fits   Named list of method fits (from run_all_methods).
#' @param view   Which omics layer's loadings to compare (name or index).
#' @param factor Which factor index to compare across methods.
#' @param cutoffs Quantile cutoffs for the Jaccard sweep.
#' @return List with $pairwise (single-cutoff metrics) and $jaccard_curve.
compute_dataset_metrics <- function(fits, view = 1, factor = 1,
                                    cutoffs = seq(0.50, 0.95, by = 0.05)) {
  # Reference method for factor alignment (first available)
  ref <- names(fits)[1]
  ref_scores <- t(fits[[ref]]$scores)   # n x k

  # Extract the chosen factor's loadings for each method, sign-aligned to ref
  loading_vecs <- lapply(names(fits), function(m) {
    al  <- align_factors(ref_scores, t(fits[[m]]$scores))
    fac <- al$matching[factor]
    L   <- fits[[m]]$loadings[[view]]
    sgn <- sign(al$correlation[factor]); if (sgn == 0) sgn <- 1
    sgn * L[, fac]
  })
  names(loading_vecs) <- names(fits)

  list(
    pairwise = compute_pairwise_metrics(
      loading_vecs,
      lapply(loading_vecs, function(v) as.integer(abs(v) >= quantile(abs(v), 0.80)))
    ),
    jaccard_curve = jaccard_across_cutoffs(loading_vecs, cutoffs = cutoffs)
  )
}

run_metrics <- function(fits_path, out_path, view = 1, factor = 1) {
  if (!file.exists(fits_path)) {
    warning("Missing fits: ", fits_path, " (run scripts/02_run_methods.R first).")
    return(invisible(NULL))
  }
  fits <- readRDS(fits_path)
  if (length(fits) < 2) { warning("Need >= 2 methods to compare in ", fits_path); return(invisible(NULL)) }
  m <- compute_dataset_metrics(fits, view = view, factor = factor)
  saveRDS(m, out_path)
  message("Saved metrics -> ", out_path)
  invisible(m)
}

if (!interactive()) {
  run_metrics("results/cll_fits.rds",       "results/cll_metrics.rds",       view = 1, factor = 1)
  run_metrics("results/radiation_fits.rds", "results/radiation_metrics.rds", view = 1, factor = 1)
}
