# simulation/generate_simdata.R
# Multi-omics data generation for the simulation study.
#
# Preferred backend: the SUMO package (Osang'ir et al. 2025, CRAN), which the
# benchmark was built around. A self-contained fallback generator is provided so
# the pipeline runs even when SUMO is unavailable; it plants low-rank factor
# structure with controllable signal-sharing across two omics layers.

#' Generate one simulated multi-omics dataset
#'
#' @param scenario Character: "unique" (each factor drives one omic),
#'                 "shared" (factors drive both omics jointly), or
#'                 "mixed" (one shared + one unique factor).
#' @param sigma    Numeric: noise standard deviation added to each entry.
#' @param config   Simulation config list (from simulation_config.R).
#' @param use_sumo Logical: use SUMO if installed. Default TRUE.
#' @return List with:
#'   \item{data}{Named list: Omic1 (p1 x n), Omic2 (p2 x n)}
#'   \item{truth}{List: Z_true (k x n), W_true (list of p_m x k)}
generate_sim_data <- function(scenario, sigma, config, use_sumo = TRUE) {
  if (use_sumo && requireNamespace("SUMO", quietly = TRUE) &&
      exists("simulateMultiOmics", where = asNamespace("SUMO"))) {
    return(.generate_sim_data_sumo(scenario, sigma, config))
  }
  .generate_sim_data_fallback(scenario, sigma, config)
}

# ---- SUMO backend ---------------------------------------------------------
.generate_sim_data_sumo <- function(scenario, sigma, config) {
  sim <- SUMO::simulateMultiOmics(
    vector_features = as.integer(config$n_features),
    n_samples       = config$n_samples,
    n_factors       = config$n_factors,
    snr             = 1 / sigma,
    signal.samples  = config$signal_prop_samples,
    signal.features = config$signal_prop_features,
    factor_structure = scenario
  )
  # Harmonize SUMO's output to the common interface used downstream.
  list(
    data  = list(Omic1 = sim$omic1, Omic2 = sim$omic2),
    truth = list(Z_true = sim$scores_true, W_true = list(sim$loadings1_true, sim$loadings2_true))
  )
}

# ---- Self-contained fallback ---------------------------------------------
.generate_sim_data_fallback <- function(scenario, sigma, config) {
  p1 <- config$n_features[["omic1"]]
  p2 <- config$n_features[["omic2"]]
  n  <- config$n_samples
  k  <- config$n_factors

  n_sig_samp  <- max(1, round(config$signal_prop_samples  * n))
  n_sig_feat1 <- max(1, round(config$signal_prop_features * p1))
  n_sig_feat2 <- max(1, round(config$signal_prop_features * p2))

  # Latent scores: each factor active in a (possibly distinct) block of samples
  Z <- matrix(0, nrow = k, ncol = n)
  for (f in seq_len(k)) {
    start <- ((f - 1) * n_sig_samp) %% n + 1
    idx   <- ((start:(start + n_sig_samp - 1) - 1) %% n) + 1
    Z[f, idx] <- rnorm(length(idx), mean = 3, sd = 1)
  }

  # Loadings: which omic each factor drives depends on the scenario
  W1 <- matrix(0, nrow = p1, ncol = k)
  W2 <- matrix(0, nrow = p2, ncol = k)
  drive <- switch(scenario,
    unique = list(omic1 = 1, omic2 = 2),                 # factor1 -> omic1, factor2 -> omic2
    shared = list(omic1 = seq_len(k), omic2 = seq_len(k)),# all factors -> both omics
    mixed  = list(omic1 = c(1, 2),   omic2 = 2),          # factor1 shared dimension
    stop("Unknown scenario: ", scenario)
  )
  for (f in drive$omic1) {
    fi <- sample.int(p1, n_sig_feat1)
    W1[fi, f] <- rnorm(n_sig_feat1, mean = 2, sd = 0.5)
  }
  for (f in drive$omic2) {
    fi <- sample.int(p2, n_sig_feat2)
    W2[fi, f] <- rnorm(n_sig_feat2, mean = 2, sd = 0.5)
  }

  Omic1 <- W1 %*% Z + matrix(rnorm(p1 * n, sd = sigma), p1, n)
  Omic2 <- W2 %*% Z + matrix(rnorm(p2 * n, sd = sigma), p2, n)

  rownames(Omic1) <- paste0("omic1_feature_", seq_len(p1))
  rownames(Omic2) <- paste0("omic2_feature_", seq_len(p2))
  colnames(Omic1) <- colnames(Omic2) <- paste0("sample_", seq_len(n))

  list(
    data  = list(Omic1 = Omic1, Omic2 = Omic2),
    truth = list(Z_true = Z, W_true = list(W1, W2))
  )
}
