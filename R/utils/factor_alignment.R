# R/utils/factor_alignment.R
# Factor matching / alignment across methods and against ground truth.
#
# Latent factors are identifiable only up to permutation and sign. Before any
# pairwise or against-truth comparison, factors must be matched (which estimated
# factor corresponds to which reference factor) and sign-flipped to agree.

#' Align a set of estimated factors to a reference set
#'
#' Matches columns of `estimated` to columns of `reference` by maximizing the
#' sum of absolute correlations (greedy / Hungarian assignment), then flips the
#' sign of each matched factor so its correlation with the reference is positive.
#'
#' @param reference Matrix (n x k_ref): reference factors (e.g., ground-truth scores).
#' @param estimated Matrix (n x k_est): estimated factors to be aligned.
#' @return List with:
#'   \item{aligned}{Matrix (n x k_ref) of sign-corrected, reordered estimated factors}
#'   \item{matching}{Integer vector: estimated-column index matched to each reference column}
#'   \item{correlation}{Numeric vector of matched (signed) correlations}
align_factors <- function(reference, estimated) {
  reference <- as.matrix(reference)
  estimated <- as.matrix(estimated)

  k_ref <- ncol(reference)
  k_est <- ncol(estimated)

  # Correlation matrix between every reference and estimated factor
  C <- suppressWarnings(cor(reference, estimated, use = "pairwise.complete.obs"))
  C[is.na(C)] <- 0

  # Assignment: prefer Hungarian (clue::solve_LSAP) on |C|; fall back to greedy.
  matching <- rep(NA_integer_, k_ref)
  if (requireNamespace("clue", quietly = TRUE) && k_ref <= k_est) {
    assign <- clue::solve_LSAP(abs(C), maximum = TRUE)
    matching <- as.integer(assign)
  } else {
    available <- seq_len(k_est)
    for (j in seq_len(k_ref)) {
      scores <- ifelse(seq_len(k_est) %in% available, abs(C[j, ]), -Inf)
      pick   <- which.max(scores)
      matching[j]  <- pick
      available    <- setdiff(available, pick)
    }
  }

  aligned     <- matrix(NA_real_, nrow = nrow(estimated), ncol = k_ref)
  correlation <- numeric(k_ref)
  for (j in seq_len(k_ref)) {
    est_col        <- estimated[, matching[j]]
    sgn            <- sign(C[j, matching[j]]); if (sgn == 0) sgn <- 1
    aligned[, j]   <- sgn * est_col
    correlation[j] <- sgn * C[j, matching[j]]
  }
  colnames(aligned) <- colnames(reference)

  list(aligned = aligned, matching = matching, correlation = correlation)
}

#' Match factors between two methods (no ground truth)
#'
#' Convenience wrapper that aligns method B's factors onto method A's.
#'
#' @param scores_a Matrix (n x k) of method A factor scores
#' @param scores_b Matrix (n x k) of method B factor scores
#' @return Result of \code{align_factors(scores_a, scores_b)}
match_methods <- function(scores_a, scores_b) {
  align_factors(reference = scores_a, estimated = scores_b)
}
