# simulation/simulation_config.R
# Central configuration for the simulation study:
# scenarios, noise levels, replicate count, factor count, and the set of
# cutoff (thresholding) functions evaluated.

source("R/utils/thresholding.R")

#' Build the cutoff function set evaluated in the simulation
#'
#' Each entry maps a continuous loading/score vector to a binary signal vector.
#' tau2 is the oracle (quantile-informed) cutoff and requires the known signal
#' proportion `pi`, supplied per factor at evaluation time.
#'
#' @return Named list of functions, each \code{function(u, pi = NULL) -> 0/1 vector}.
build_cutoff_fns <- function() {
  list(
    tau1 = function(u, pi = NULL) threshold_tau1(u),
    tau2 = function(u, pi = NULL) if (!is.null(pi)) threshold_tau2(u, pi) else threshold_tau3(u),
    tau3 = function(u, pi = NULL) threshold_tau3(u),
    tau4 = function(u, pi = NULL) threshold_tau4(u),
    tau5 = function(u, pi = NULL) threshold_tau5(u),
    tau6 = function(u, pi = NULL) threshold_tau6(u)
  )
}

#' Default simulation configuration
#'
#' @param n_replicates Replicates per (scenario, sigma) cell. Default 100.
#' @param n_factors    Number of latent factors simulated and fit. Default 2.
#' @return A configuration list consumed by run_full_simulation().
get_simulation_config <- function(n_replicates = 100, n_factors = 2) {
  list(
    # Signal-sharing scenarios across the two omics layers
    scenarios = c("unique", "shared", "mixed"),

    # 29 noise levels (signal-to-noise sweep)
    sigmas = round(seq(0.5, 8.5, length.out = 29), 3),

    n_replicates = n_replicates,
    n_factors    = n_factors,

    # Dimensions of the simulated omics matrices
    n_features = c(omic1 = 4000, omic2 = 3000),
    n_samples  = 50,

    # Proportion of true signal features / samples per factor
    signal_prop_features = 0.025,   # ~100 of 4000
    signal_prop_samples  = 0.18,    # ~9 of 50

    # Cutoff strategies evaluated
    cutoff_fns = build_cutoff_fns()
  )
}
