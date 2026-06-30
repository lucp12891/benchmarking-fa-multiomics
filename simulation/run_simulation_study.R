# simulation/run_simulation_study.R
# Full simulation pipeline: 3 scenarios × 29 noise levels × 100 replicates
# Uses SUMO (Osang'ir et al. 2025) for data generation

source("R/utils/metrics.R")
source("R/utils/thresholding.R")
source("R/utils/factor_alignment.R")
source("R/methods/run_mofa.R")
source("R/methods/run_mfa.R")
source("R/methods/run_gfa.R")
source("R/methods/run_fabia.R")
source("simulation/simulation_config.R")

#' Run one replicate of the simulation study
#'
#' @param scenario  Character: "unique", "shared", or "mixed"
#' @param sigma     Numeric: noise standard deviation
#' @param seed      Integer: random seed
#' @param config    List of simulation parameters (from simulation_config.R)
#' @return Data frame of performance metrics for all methods × cutoffs
run_one_replicate <- function(scenario, sigma, seed, config) {
  set.seed(seed)

  # --- 1. Generate data ---
  sim_data <- generate_sim_data(scenario = scenario, sigma = sigma, config = config)
  omics_list  <- sim_data$data       # named list: Omic1, Omic2
  ground_truth <- sim_data$truth     # list: Z_true (k × n), W_true (list of p_m × k)

  # --- 2. Run all four methods ---
  methods <- list(
    MOFA  = tryCatch(run_mofa_multiomics(omics_list,  n_factors = config$n_factors), error = function(e) NULL),
    MFA   = tryCatch(run_mfa_multiomics(omics_list,   n_factors = config$n_factors), error = function(e) NULL),
    GFA   = tryCatch(run_gfa_multiomics(omics_list,   n_factors = config$n_factors, n_init = 1), error = function(e) NULL),
    FABIA = tryCatch(run_fabia_multiomics(omics_list, n_factors = config$n_factors), error = function(e) NULL)
  )
  methods <- Filter(Negate(is.null), methods)
  if (length(methods) == 0) return(NULL)

  # --- 3. Align factors to ground truth ---
  gt_scores <- ground_truth$Z_true  # k × n

  results <- lapply(names(methods), function(mname) {
    fit <- methods[[mname]]

    # Align scores to ground truth
    aligned_scores <- tryCatch(
      align_factors(t(gt_scores), t(fit$scores))$aligned,
      error = function(e) t(fit$scores)
    )

    # Per-factor, per-cutoff metrics
    factor_results <- lapply(seq_len(config$n_factors), function(fac) {
      true_z  <- as.integer(abs(gt_scores[fac, ]) > 0)  # binary ground truth
      pred_z  <- aligned_scores[, fac]

      cutoff_results <- lapply(config$cutoff_fns, function(tau_fn) {
        pred_bin <- tau_fn(pred_z)
        jac <- jaccard_index(pred_bin, true_z)
        cls <- classification_metrics(pred_bin, true_z)
        pcc <- cor(pred_z, gt_scores[fac, ], use = "complete.obs")

        data.frame(
          method = mname, factor = fac,
          jaccard = jac, pearson = pcc,
          sensitivity = cls["sensitivity"],
          specificity = cls["specificity"],
          precision   = cls["precision"],
          F1          = cls["F1"],
          stringsAsFactors = FALSE
        )
      })
      do.call(rbind, cutoff_results)
    })
    do.call(rbind, factor_results)
  })

  result_df <- do.call(rbind, results)
  result_df$scenario <- scenario
  result_df$sigma    <- sigma
  result_df$seed     <- seed
  result_df
}

#' Run full simulation study across all scenarios and noise levels
#'
#' @param config  Simulation config list (from simulation_config.R)
#' @param outdir  Directory to save per-scenario results
run_full_simulation <- function(config, outdir = "results/simulation") {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  all_results <- list()
  total <- length(config$scenarios) * length(config$sigmas) * config$n_replicates
  message(sprintf("Starting simulation: %d total runs", total))

  for (scenario in config$scenarios) {
    for (sigma in config$sigmas) {
      cat(sprintf("\n[%s | sigma=%.1f] ", scenario, sigma))

      rep_results <- lapply(seq_len(config$n_replicates), function(rep) {
        cat(".")
        run_one_replicate(
          scenario = scenario,
          sigma    = sigma,
          seed     = rep * 100 + as.integer(sigma),
          config   = config
        )
      })
      rep_results <- Filter(Negate(is.null), rep_results)

      scenario_df <- do.call(rbind, rep_results)
      fname <- file.path(outdir, sprintf("sim_%s_sigma%.0f.rds", scenario, sigma))
      saveRDS(scenario_df, fname)
      all_results[[length(all_results) + 1]] <- scenario_df
    }
  }

  final <- do.call(rbind, all_results)
  saveRDS(final, file.path(outdir, "simulation_all_results.rds"))
  message("\nSimulation complete. Results saved to: ", outdir)
  invisible(final)
}

# --- Entry point ---
if (!interactive()) {
  config <- get_simulation_config()
  run_full_simulation(config)
}
