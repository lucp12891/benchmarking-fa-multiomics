# scripts/05_simulation_study.R
# End-to-end driver for the simulation study.
# Wraps simulation/run_simulation_study.R so the whole study runs from one call.
#
# WARNING: the full grid (3 scenarios x 29 noise levels x 100 replicates x 4 methods)
# is computationally heavy (hours). Use the `quick` flag for a smoke test.
#
# Output: results/simulation/*.rds (per-cell) and simulation_all_results.rds

source("R/utils/libs.R")
source("simulation/generate_simdata.R")
source("simulation/simulation_config.R")
source("simulation/run_simulation_study.R")

run_simulation_pipeline <- function(quick = FALSE) {
  config <- if (quick) {
    cfg <- get_simulation_config(n_replicates = 2, n_factors = 2)
    cfg$sigmas    <- c(1, 4, 8)
    cfg$scenarios <- c("unique", "shared")
    cfg
  } else {
    get_simulation_config(n_replicates = 100, n_factors = 2)
  }

  run_full_simulation(config, outdir = "results/simulation")
}

if (!interactive()) {
  quick <- isTRUE(as.logical(Sys.getenv("SIM_QUICK", "FALSE")))
  run_simulation_pipeline(quick = quick)
}
