# R/utils/metrics.R
# Agreement metrics: Jaccard index, cosine similarity, Pearson correlation
# Used across both real-data benchmarking and simulation studies

#' Compute Jaccard similarity between two binary vectors
#'
#' @param x Binary vector (0/1)
#' @param y Binary vector (0/1)
#' @return Scalar Jaccard index
jaccard_index <- function(x, y) {
  intersection <- sum(x == 1 & y == 1)
  union        <- sum(x == 1 | y == 1)
  if (union == 0) return(NA_real_)
  intersection / union
}

#' Compute cosine similarity between two numeric vectors
#'
#' @param x Numeric vector
#' @param y Numeric vector
#' @return Scalar cosine similarity
cosine_similarity <- function(x, y) {
  denom <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  if (denom == 0) return(NA_real_)
  sum(x * y) / denom
}

#' Compute all pairwise agreement metrics for a list of factor vectors
#'
#' @param factor_list Named list of numeric vectors (one per method)
#' @param binary_list Named list of binary vectors (binarized at cutoff τ)
#' @return Data frame with pairwise Jaccard, cosine, and Pearson values
compute_pairwise_metrics <- function(factor_list, binary_list) {
  methods <- names(factor_list)
  pairs   <- combn(methods, 2, simplify = FALSE)

  results <- lapply(pairs, function(p) {
    m1 <- p[1]; m2 <- p[2]
    data.frame(
      method1  = m1,
      method2  = m2,
      jaccard  = jaccard_index(binary_list[[m1]], binary_list[[m2]]),
      cosine   = cosine_similarity(factor_list[[m1]], factor_list[[m2]]),
      pearson  = cor(factor_list[[m1]], factor_list[[m2]], method = "pearson"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, results)
}

#' Compute Jaccard index across a range of quantile cutoffs
#'
#' @param factor_list Named list of numeric loading/score vectors (one per method)
#' @param cutoffs Numeric vector of quantile thresholds (e.g., seq(0.50, 0.95, by = 0.05))
#' @return Data frame: cutoff × method pair × Jaccard index
jaccard_across_cutoffs <- function(factor_list, cutoffs = seq(0.50, 0.95, by = 0.05)) {
  methods <- names(factor_list)
  pairs   <- combn(methods, 2, simplify = FALSE)

  results <- lapply(cutoffs, function(q) {
    binary_list <- lapply(factor_list, function(v) {
      threshold <- quantile(abs(v), probs = q)
      as.integer(abs(v) >= threshold)
    })
    metrics <- compute_pairwise_metrics(factor_list, binary_list)
    metrics$cutoff <- q
    metrics
  })
  do.call(rbind, results)
}

#' Classification metrics against ground truth (simulation use)
#'
#' @param predicted Binary vector of predicted signal indicators
#' @param truth     Binary vector of true signal indicators
#' @return Named numeric vector: sensitivity, specificity, precision, F1
classification_metrics <- function(predicted, truth) {
  TP <- sum(predicted == 1 & truth == 1)
  TN <- sum(predicted == 0 & truth == 0)
  FP <- sum(predicted == 1 & truth == 0)
  FN <- sum(predicted == 0 & truth == 1)

  sensitivity <- TP / (TP + FN + 1e-10)
  specificity <- TN / (TN + FP + 1e-10)
  precision   <- TP / (TP + FP + 1e-10)
  f1          <- 2 * precision * sensitivity / (precision + sensitivity + 1e-10)

  c(sensitivity = sensitivity, specificity = specificity,
    precision = precision, F1 = f1)
}
