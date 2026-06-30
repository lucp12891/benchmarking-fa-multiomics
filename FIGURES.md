# Figure → Code Map

Every main-text figure in *Osang'ir et al. (2026)* and the script(s) that produce
it. File names are exactly as used during the analysis (unchanged). Folders follow
the manuscript sections.

> Figure 1 Panel A (schematic) is a drawn diagram, not script output.

## Methods / Simulation framework

| Figure | Description | Script(s) |
|--------|-------------|-----------|
| 1B–C | Concatenated-omics heatmaps for the three scenarios; signal-to-noise ratio vs noise variance | [`simulation_study/framework/simulate_nonoverlapping_signals_varying.R`](simulation_study/framework/simulate_nonoverlapping_signals_varying.R), [`trysim_image_plots.R`](simulation_study/framework/trysim_image_plots.R), [`snr_plot.R`](simulation_study/framework/snr_plot.R), [`snr_with_shift.R`](simulation_study/framework/snr_with_shift.R) |

## Application to CLL Data (`cll_study/`)

| Figure | Description | Script(s) |
|--------|-------------|-----------|
| 2 | Factor 1 scores: variance explained, sample scores (IGHV), Pearson + cosine | [`cll_study/cll_pipeline.R`](cll_study/cll_pipeline.R), [`cosine_cll.R`](cll_study/cosine_cll.R) |
| 3 | Factor 1 loadings: Pearson + cosine (CpG & drug-response) | [`cll_study/cosine_cll.R`](cll_study/cosine_cll.R) |
| 4 | Detected **samples** agreement: Jaccard across cutoffs, UpSet, Venn | [`cll_study/jaccard_cll.R`](cll_study/jaccard_cll.R), [`venn_diagram.R`](cll_study/venn_diagram.R) |
| 5 | Detected **features** agreement: Jaccard (methylation & drugs) | [`cll_study/jaccard_cll.R`](cll_study/jaccard_cll.R) |
| 6 | FABIA sparsity control over (α, spz, spl) | [`cll_study/monitor_fabia_parameters_cll.R`](cll_study/monitor_fabia_parameters_cll.R), [`alpha_spl_spz_cll.R`](cll_study/alpha_spl_spz_cll.R), [`cll_pipeline_fabia_cutoff.R`](cll_study/cll_pipeline_fabia_cutoff.R) |
| 7 | Top-10 loadings (lollipop) + cross-method concordance after cutoff | [`cll_study/lollipop plots of loadings.R`](cll_study/lollipop%20plots%20of%20loadings.R), [`cosine_cll.R`](cll_study/cosine_cll.R), [`cll_pipeline_fabia_cutoff.R`](cll_study/cll_pipeline_fabia_cutoff.R) |

Supporting: [`app_to_cll_data.R`](cll_study/app_to_cll_data.R) (data prep), [`analysis_cll.R`](cll_study/analysis_cll.R).

## Application to Radiation Data (`radiation_study/`)

| Figure | Description | Script(s) |
|--------|-------------|-----------|
| 8 | Factor 1 scores: variance explained, sample scores (folic-acid group), Pearson + cosine | [`radiation_study/app_to_rad_data.R`](radiation_study/app_to_rad_data.R), [`cosine_rad.R`](radiation_study/cosine_rad.R) |
| 9 | Factor 1 loadings: Pearson + cosine (genes & proteins) | [`radiation_study/cosine_rad.R`](radiation_study/cosine_rad.R) |

Supporting: [`rad_pipeline.R`](radiation_study/rad_pipeline.R), [`jaccard_rad.R`](radiation_study/jaccard_rad.R), [`venn_diagram_rad.R`](radiation_study/venn_diagram_rad.R), [`alpha_spl_spz_rad.R`](radiation_study/alpha_spl_spz_rad.R) (FABIA cutoff), [`lollipop plots of loadings_rad.R`](radiation_study/lollipop%20plots%20of%20loadings_rad.R), [`analysis_rad.R`](radiation_study/analysis_rad.R).

## Simulation Study (`simulation_study/`)

| Figure | Description | Script(s) |
|--------|-------------|-----------|
| 10 | Single simulation instance (Setting C): Jaccard, F1, precision, sensitivity, specificity, GT Pearson | [`simulation_study/single_instance/Paper_II_Analysis_One_Instance_App.R`](simulation_study/single_instance/Paper_II_Analysis_One_Instance_App.R), [`one_data_example_appendix.R`](simulation_study/single_instance/one_data_example_appendix.R), [`latent.R`](simulation_study/single_instance/latent.R), [`parallel_single_unique.R`](simulation_study/single_instance/parallel_single_unique.R) |
| 11 | Jaccard for Omic-1 loadings across cutoff strategies × noise levels | [`simulation_study/full_grid/utilitiesR/plot_jaccard_index_sim.R`](simulation_study/full_grid/utilitiesR/plot_jaccard_index_sim.R), [`varphi_multisim.R`](simulation_study/full_grid/utilitiesR/varphi_multisim.R), [`retrieve_JI_PM_data.R`](simulation_study/full_grid/utilitiesR/retrieve_JI_PM_data.R) + the HPC runners below |

### Simulation engine

- **Shared function library** (`functions/`): [`libs.R`](functions/libs.R), [`factor_methods.R`](functions/factor_methods.R) (FABIA/MOFA/MFA/GFA), [`similarity_metrics.R`](functions/similarity_metrics.R), [`cutoffs.R`](functions/cutoffs.R) (the 12 `varphi` cutoff rules), [`helpers.R`](functions/helpers.R), [`multiple_factor.R`](functions/multiple_factor.R).
- **HPC grid runs** (`simulation_study/full_grid/`): per-scenario array jobs — `singleR/parallel_single_sim_*.R` (omic-specific), `sharedR/parallel_shared_sim_*.R` (shared), `mixedR/` + `mixedR/hpc/parallel_multi_sim_*.R` and `simulate_MOmics_main*.R` (mixed). These generate the per-noise-level results feeding Figure 11 and Supplementary Figures S28–S60.
- **Result retrieval & plotting** (`simulation_study/full_grid/utilitiesR/`): `retrieve_JI_PM_data.R`, `retrieve_JI_PM_shared_data.R`, `plot_jaccard_index_sim.R`, `plot_performanceMetrics.R`, `varphi_multisim.R`, `trysim_image_plots.R`.

## Notes

- Where a script appeared in several dated working folders, the **most recent**
  version was taken (CLL: *Mar 02, 2026*; single instance: *Jan 27, 2026*;
  full grid: *Run_02_09_2025*). Earlier snapshots and scratch `dumpsite/` copies
  were not carried over.
- Tables 2 and 3 (Jaccard at q = 0.60) are produced alongside Figures 4–5 / 8–9
  by the same `jaccard_*` scripts.
