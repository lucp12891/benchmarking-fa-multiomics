# Benchmarking Factor-Based Models for Multi-Omics Integration

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-%3E%3D4.2.0-blue)](https://www.r-project.org/)
[![SUMO on CRAN](https://img.shields.io/cran/v/SUMO)](https://CRAN.R-project.org/package=SUMO)

> **Osang'ir et al. (2026)** — *Benchmarking Factor-Based Models for Integrative Analysis of Large-Scale Biomolecular Data: Performance, Robustness, and Interpretability*

This repository contains the analysis code behind the paper, which benchmarks four
factor analysis (FA)-based methods — **MOFA**, **MFA**, **GFA**, and **FABIA** —
across two real-world multi-omics datasets (CLL, radiation) and a systematic
simulation study.

The code is organized by **manuscript section**. For an exact figure-by-figure
map, see **[FIGURES.md](FIGURES.md)**.

---

## Repository Structure

```
benchmarking-fa-multiomics/
├── README.md
├── FIGURES.md                  # Figure → script map (main text)
├── LICENSE  ·  CITATION.cff
│
├── functions/                  # Shared simulation engine
│   ├── libs.R                  # Package loading
│   ├── factor_methods.R        # FABIA / MOFA / MFA / GFA fitting + extractors
│   ├── similarity_metrics.R    # Jaccard, cosine, classification metrics
│   ├── cutoffs.R               # The varphi cutoff rules (τ)
│   ├── helpers.R
│   └── multiple_factor.R
│
├── cll_study/                  # Application to CLL data — Figures 2–7, Table 2
│   ├── app_to_cll_data.R       # Data preparation
│   ├── cll_pipeline.R          # Fit all four methods
│   ├── cll_pipeline_fabia_cutoff.R
│   ├── cosine_cll.R            # Figs 2, 3, 7  (Pearson/cosine score & loading agreement)
│   ├── jaccard_cll.R           # Figs 4, 5     (Jaccard across cutoffs)
│   ├── venn_diagram.R          # Fig 4         (Venn / overlap)
│   ├── monitor_fabia_parameters_cll.R , alpha_spl_spz_cll.R   # Fig 6 (FABIA sparsity)
│   ├── lollipop plots of loadings.R                            # Fig 7 (top-10 loadings)
│   └── analysis_cll.R , monitor_fabia_parameters.R
│
├── radiation_study/            # Application to radiation data — Figures 8–9, Table 3
│   ├── app_to_rad_data.R       # Data prep + fitting
│   ├── rad_pipeline.R
│   ├── cosine_rad.R            # Figs 8, 9
│   ├── jaccard_rad.R , venn_diagram_rad.R
│   ├── alpha_spl_spz_rad.R     # FABIA sparsity (radiation)
│   ├── lollipop plots of loadings_rad.R
│   └── analysis_rad.R
│
└── simulation_study/           # Simulation study — Figures 1, 10, 11 + Supplementary
    ├── framework/              # Fig 1   (scenario heatmaps, SNR)
    ├── single_instance/        # Fig 10  (one Setting-C instance, all metrics)
    └── full_grid/              # Fig 11 + Suppl S28–S60 (HPC grid over noise levels)
        ├── singleR/  sharedR/  mixedR/   # per-scenario array jobs
        ├── mixedR/hpc/                   # simulate_MOmics_main*, parallel_multi_sim_*
        └── utilitiesR/                   # result retrieval + plotting
```

---

## Datasets

### CLL Dataset
- **Source**: Available via the [MOFAdata](https://bioconductor.org/packages/MOFAdata) Bioconductor package
- **Modalities**: DNA methylation + drug response, n = 200 samples
- Loaded in `cll_study/app_to_cll_data.R`

### Radiation Exposure Dataset
- **Transcriptomics**: NCBI SRA BioProject [PRJNA1431456](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1431456)
- **Proteomics**: PRIDE repository [PXD075411](https://www.ebi.ac.uk/pride/archive/projects/PXD075411)
- **Modalities**: RNA-seq (18,329 genes) + proteomics (3,011 proteins), n = 16 samples
- Loaded in `radiation_study/app_to_rad_data.R`

---

## Installation

```r
install.packages(c("SUMO", "FactoMineR", "fabia", "GFA", "ggplot2",
                   "ggupset", "VennDiagram", "pheatmap", "dplyr", "tidyr",
                   "data.table", "pROC", "mclust", "zoo"))

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("MOFA2", "MOFAdata"))
```

R version ≥ 4.2.0 is recommended. Simulation data are generated with
[**SUMO**](https://CRAN.R-project.org/package=SUMO) (Osang'ir et al., 2025).

---

## Methods Overview

| Method | Type | Inference | Sparsity | Reference |
|--------|------|-----------|----------|-----------|
| MOFA | Probabilistic FA | Variational Bayes | ARD + element-wise | Argelaguet et al. (2018) |
| MFA | Classical PCA extension | SVD | None | Abdi et al. (2007) |
| GFA | Bayesian group FA | Approx. Bayes | ARD (dataset-wise) | Klami et al. (2014) |
| FABIA | Sparse biclustering FA | EM algorithm | Laplace (loadings + scores) | Hochreiter et al. (2010) |

---

## Citation

See [`CITATION.cff`](CITATION.cff).

```bibtex
@article{osangir2026benchmarking,
  title   = {Benchmarking Factor-Based Models for Integrative Analysis of
             Large-Scale Biomolecular Data: Performance, Robustness, and Interpretability},
  author  = {Osang'ir, Bernard Isekah and Ahmed, Mohamed Mysara and Mastroleo, Felice
             and Leys, Natalie and Claesen, Jurgen and Gupta, Surya and Shkedy, Ziv},
  journal = {BMC Bioinformatics},
  year    = {2026}
}
```

## License

MIT License. See `LICENSE`.

## Contact

Correspondence: **Surya Gupta** (surya.gupta@sckcen.be) · **Ziv Shkedy** (ziv.shkedy@uhasselt.be)
