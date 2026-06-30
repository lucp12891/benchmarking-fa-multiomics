# R/utils/preprocessing.R
# Data loading, harmonization, and normalization helpers shared by the
# real-data benchmarking scripts.

#' Harmonize a set of omics matrices to a common sample set
#'
#' Subsets and reorders every matrix to the intersection of sample names
#' (columns), so all views share identical, identically-ordered samples.
#'
#' @param omics_list Named list of matrices, each (features x samples) with column names.
#' @return Named list of matrices restricted to common samples, same column order.
harmonize_samples <- function(omics_list) {
  sample_sets <- lapply(omics_list, colnames)
  if (any(vapply(sample_sets, is.null, logical(1)))) {
    stop("Every omics matrix must have column (sample) names.")
  }
  common <- Reduce(intersect, sample_sets)
  if (length(common) == 0) stop("No samples shared across all omics layers.")
  lapply(omics_list, function(m) m[, common, drop = FALSE])
}

#' Feature-wise standardization (z-score) of a single omics matrix
#'
#' Scales each feature (row) to zero mean and unit variance. Features with zero
#' variance are set to 0. NA values are imputed to 0 after scaling.
#'
#' @param x Matrix (features x samples)
#' @return Standardized matrix of the same shape.
scale_features <- function(x) {
  x <- as.matrix(x)
  z <- t(scale(t(x)))          # scale operates column-wise; transpose to scale features
  z[is.na(z)] <- 0
  z
}

#' Filter to the most variable features
#'
#' @param x      Matrix (features x samples)
#' @param top_n  Keep this many highest-variance features. If NULL, keep all.
#' @return Filtered matrix (features x samples).
filter_top_variable <- function(x, top_n = NULL) {
  if (is.null(top_n) || top_n >= nrow(x)) return(x)
  v    <- apply(x, 1, var, na.rm = TRUE)
  keep <- order(v, decreasing = TRUE)[seq_len(top_n)]
  x[sort(keep), , drop = FALSE]
}

#' End-to-end preprocessing of a multi-omics list
#'
#' Harmonizes samples, optionally filters to top-variable features per view,
#' and feature-standardizes each view.
#'
#' @param omics_list Named list of (features x samples) matrices.
#' @param top_n      Named integer vector or single integer of top features per view (optional).
#' @param standardize Logical. Apply feature-wise z-scoring. Default TRUE.
#' @return Preprocessed named list of matrices.
preprocess_omics <- function(omics_list, top_n = NULL, standardize = TRUE) {
  omics_list <- harmonize_samples(omics_list)

  if (!is.null(top_n)) {
    tn <- if (length(top_n) == 1) setNames(rep(top_n, length(omics_list)), names(omics_list)) else top_n
    omics_list <- lapply(names(omics_list), function(nm) filter_top_variable(omics_list[[nm]], tn[[nm]]))
    names(omics_list) <- names(tn)
  }

  if (standardize) omics_list <- lapply(omics_list, scale_features)
  omics_list
}
