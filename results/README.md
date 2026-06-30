# Results

Saved numerical outputs from the analysis scripts. These files are generated, not
committed (see `.gitignore`).

| File | Produced by | Contents |
|------|-------------|----------|
| `cll_fits.rds`              | `02_run_methods.R`     | MOFA/MFA/GFA/FABIA fits on the CLL data |
| `radiation_fits.rds`        | `02_run_methods.R`     | MOFA/MFA/GFA/FABIA fits on the radiation data |
| `cll_metrics.rds`           | `03_compute_metrics.R` | Pairwise agreement + Jaccard-vs-cutoff curves (CLL) |
| `radiation_metrics.rds`     | `03_compute_metrics.R` | Pairwise agreement + Jaccard-vs-cutoff curves (radiation) |
| `simulation/sim_<scenario>_sigma<σ>.rds` | `05_simulation_study.R` | Per-cell simulation metrics |
| `simulation/simulation_all_results.rds`  | `05_simulation_study.R` | Concatenated full-study metrics |

Each fit object is a named list (one entry per method):
`list(scores = k x n, loadings = named list of (features_m x k), model = raw fit)`.
