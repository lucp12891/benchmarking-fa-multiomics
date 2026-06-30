# R/utils/thresholding.R
# Six cutoff strategies (τ1–τ6) for binarising continuous factor loadings/scores
# See Table 1 of Osang'ir et al. (2026)

#' τ1 — Normalized spread heuristic
#'
#' Normalizes vector to [0,1] and applies midpoint of the normalized range.
#' Scale-invariant; exploits post-normalization separation.
#'
#' @param u Numeric vector of loadings or scores
#' @return Binary signal vector (0/1)
threshold_tau1 <- function(u) {
  u_norm <- (u - min(u)) / (max(u) - min(u) + 1e-12)
  tau    <- 0.5 * (max(u_norm) - min(u_norm))
  as.integer(u_norm >= tau)
}

#' τ2 — Quantile-informed cutoff (oracle / simulation use)
#'
#' Uses known signal proportion π to define the cutoff.
#' Matches true sparsity — ceiling benchmark for simulation studies.
#'
#' @param u Numeric vector of loadings or scores
#' @param pi Known proportion of true signal features/samples (0 < pi < 1)
#' @return Binary signal vector (0/1)
threshold_tau2 <- function(u, pi) {
  tau <- quantile(abs(u), probs = 1 - pi)
  as.integer(abs(u) >= tau)
}

#' τ3 — Fixed 80th percentile
#'
#' Simple, scale-free rule retaining the top 20% of absolute values.
#' Default starting point recommended for exploratory analyses.
#'
#' @param u Numeric vector of loadings or scores
#' @return Binary signal vector (0/1)
threshold_tau3 <- function(u) {
  tau <- quantile(abs(u), probs = 0.80)
  as.integer(abs(u) >= tau)
}

#' τ4 — Robust Z/MAD rule
#'
#' Flags entries exceeding median + c * MAD. Robust to outliers.
#' c = 3 balances adaptivity and interpretability.
#'
#' @param u Numeric vector of loadings or scores
#' @param c Multiplier for MAD (default 3)
#' @return Binary signal vector (0/1)
threshold_tau4 <- function(u, c = 3) {
  center <- median(abs(u))
  mad_u  <- mad(abs(u), constant = 1)
  tau    <- center + c * mad_u
  as.integer(abs(u) >= tau)
}

#' τ5 — Empirical-null / FDR-based cutoff
#'
#' Fits a two-component mixture to |u|, identifies background null,
#' and calls signal at FDR ≤ alpha using the empirical null distribution.
#'
#' @param u     Numeric vector of loadings or scores
#' @param alpha FDR threshold (default 0.05)
#' @return Binary signal vector (0/1)
threshold_tau5 <- function(u, alpha = 0.05) {
  abs_u <- abs(u)

  # Fit two-component Gaussian mixture via EM (null + signal)
  # Fallback to simple percentile if mixture fails
  tryCatch({
    if (!requireNamespace("mclust", quietly = TRUE)) {
      warning("mclust not installed; falling back to tau3 for tau5.")
      return(threshold_tau3(u))
    }
    fit <- mclust::densityMclust(abs_u, G = 2, verbose = FALSE)
    # Identify null component (lower mean)
    null_comp <- which.min(fit$parameters$mean)
    # Posterior probability of belonging to signal component
    prob_signal <- fit$z[, -null_comp, drop = FALSE]
    if (ncol(prob_signal) > 1) prob_signal <- rowSums(prob_signal)
    # FDR = expected proportion of null among called positives
    # Call signal where local FDR <= alpha
    local_fdr <- 1 - prob_signal
    as.integer(local_fdr <= alpha)
  }, error = function(e) {
    warning("tau5 mixture fitting failed; falling back to tau3: ", conditionMessage(e))
    threshold_tau3(u)
  })
}

#' τ6 — Rolling-mean filter (locally adaptive)
#'
#' Flags entries exceeding their local neighbourhood mean.
#' Robust under structured noise at high sigma. Window size delta = 31.
#'
#' @param u     Numeric vector of loadings or scores
#' @param delta Rolling window size (default 31)
#' @return Binary signal vector (0/1)
threshold_tau6 <- function(u, delta = 31) {
  n     <- length(u)
  half  <- floor(delta / 2)
  R     <- numeric(n)
  for (i in seq_len(n)) {
    lo   <- max(1, i - half)
    hi   <- min(n, i + half)
    R[i] <- mean(u[lo:hi])
  }
  b       <- as.integer(u > R)
  flagged <- which(b == 1)
  if (length(flagged) == 0) return(b)
  tau <- mean(u[flagged])
  as.integer(u >= tau)
}

#' Apply all six threshold strategies to a vector
#'
#' @param u  Numeric vector of loadings or scores
#' @param pi Known signal proportion for tau2 (required for simulation use)
#' @return Data frame: one column per strategy (tau1 … tau6)
apply_all_thresholds <- function(u, pi = NULL) {
  res <- data.frame(
    tau1 = threshold_tau1(u),
    tau2 = if (!is.null(pi)) threshold_tau2(u, pi) else NA_integer_,
    tau3 = threshold_tau3(u),
    tau4 = threshold_tau4(u),
    tau5 = threshold_tau5(u),
    tau6 = threshold_tau6(u)
  )
  res
}
