# Benchmarking Factor-Based Models for Multi-Omics Integration

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-%3E%3D4.2.0-blue)](https://www.r-project.org/)
[![SUMO on CRAN](https://img.shields.io/cran/v/SUMO)](https://CRAN.R-project.org/package=SUMO)

> **Osang'ir et al. (2026)** — *Benchmarking Factor-Based Models for Integrative Analysis of Large-Scale Biomolecular Data: Performance, Robustness, and Interpretability*

This repository contains all code to reproduce the analyses, figures, and simulation results from the paper. It benchmarks four factor analysis (FA)-based methods — **MOFA**, **MFA**, **GFA**, and **FABIA** — across two real-world multi-omics datasets and a systematic simulation study.

---

## Repository Structure

```
benchmarking-fa-multiomics/
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
│
├── data/
│   ├── README.md               # Instructions to download CLL and radiation data
│   ├── CLL/                    # CLL dataset (downloaded via MOFAdata)
│   └── radiation/              # Radiation dataset (NCBI + PRIDE links)
│
├── R/
│   ├── utils/
│   │   ├── metrics.R           # Jaccard, cosine similarity, Pearson correlation
│   │   ├── thresholding.R      # Six cutoff strategies (τ1–τ6)
│   │   ├── factor_alignment.R  # Factor matching across methods
│   │   └── preprocessing.R     # Data loading and normalization
│   ├── methods/
│   │   ├── run_mofa.R          # MOFA wrapper
│   │   ├── run_mfa.R           # MFA wrapper
│   │   ├── run_gfa.R           # GFA wrapper
│   │   └── run_fabia.R         # FABIA wrapper (multi-omics adaptation)
│   └── plots/
│       ├── plot_scores.R       # Factor score plots (Figures 3, 10)
│       ├── plot_loadings.R     # Factor loading heatmaps (Figures 4, 11)
│       ├── plot_jaccard.R      # Jaccard index curves (Figures 5–7)
│       ├── plot_simulation.R   # Simulation performance plots (Figures 12–13)
│       └── plot_fabia_sparsity.R # FABIA sparsity grid (Figure 8)
│
├── simulation/
│   ├── simulation_config.R     # Parameters: noise levels, scenarios, replicates
│   ├── generate_simdata.R      # SUMO-based data generation
│   └── run_simulation_study.R  # Full simulation pipeline (100 replicates × σ)
│
├── scripts/
│   ├── 01_prepare_data.R       # Load and preprocess CLL and radiation data
│   ├── 02_run_methods.R        # Run all four FA methods on both datasets
│   ├── 03_compute_metrics.R    # Compute all agreement metrics
│   ├── 04_generate_plots.R     # Reproduce all main text figures
│   └── 05_simulation_study.R   # End-to-end simulation analysis
│
├── results/
│   └── README.md               # Description of saved output files
│
└── figures/
    └── README.md               # Description of all generated figures
```

---

## Datasets

### CLL Dataset
- **Source**: Available via the [MOFAdata](https://bioconductor.org/packages/MOFAdata) Bioconductor package
- **Modalities**: DNA methylation (4,248 CpG sites) + drug response (310 compounds), n = 200 samples
- **Download**: Run `scripts/01_prepare_data.R` — automatically fetches via `MOFAdata`

### Radiation Exposure Dataset
- **Transcriptomics**: NCBI SRA BioProject [PRJNA1431456](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1431456)
- **Proteomics**: PRIDE repository [PXD075411](https://www.ebi.ac.uk/pride/archive/projects/PXD075411)
- **Modalities**: RNA-seq (18,329 genes) + proteomics (3,011 proteins), n = 16 samples
- See `data/radiation/README.md` for full download and preprocessing instructions

---

## Installation

```r
# Install required packages
install.packages(c("SUMO", "FactoMineR", "FABIA", "GFA", "ggplot2",
                   "ggupset", "VennDiagram", "pheatmap", "dplyr", "tidyr"))

# Bioconductor packages
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("MOFA2", "MOFAdata"))
```

R version ≥ 4.2.0 is recommended. A full `renv.lock` is provided for exact reproducibility.

```r
# Restore exact package environment
install.packages("renv")
renv::restore()
```

---

## Reproducing the Paper

Run scripts in order:

```r
source("scripts/01_prepare_data.R")       # ~2 min
source("scripts/02_run_methods.R")        # ~5 min
source("scripts/03_compute_metrics.R")   # ~2 min
source("scripts/04_generate_plots.R")    # ~3 min
source("scripts/05_simulation_study.R")  # ~2–4 hours (100 replicates × 29 noise levels)
```

All figures are saved to `figures/`. All numerical results are saved to `results/`.

---

## Methods Overview

| Method | Type | Inference | Sparsity | Reference |
|--------|------|-----------|----------|-----------|
| MOFA | Probabilistic FA | Variational Bayes | ARD (dataset-wise) + element-wise | Argelaguet et al. (2018) |
| MFA | Classical PCA extension | SVD | None | Abdi et al. (2007) |
| GFA | Bayesian group FA | Approx. Bayes | ARD (dataset-wise) | Klami et al. (2014) |
| FABIA | Sparse biclustering FA | EM algorithm | Laplace (loadings + scores) | Hochreiter et al. (2010) |

---

## Thresholding Strategies

Six complementary cutoff strategies (τ1–τ6) are implemented in `R/utils/thresholding.R`:

| ID | Name | Rationale |
|----|------|-----------|
| τ1 | Normalized spread | Scale-invariant heuristic |
| τ2 | Quantile-informed | Oracle benchmark matching known sparsity |
| τ3 | Fixed 80th percentile | Simple scale-free rule |
| τ4 | Z/MAD rule | Robust deviation-from-center |
| τ5 | Empirical null / FDR | Controls false discoveries |
| τ6 | Rolling-mean filter | Locally adaptive; robust under structured noise |

---

## Citation

If you use this code, please cite:

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

---

## License

MIT License. See `LICENSE` for details.

## Contact

Correspondence: **Surya Gupta** (surya.gupta@sckcen.be) · **Ziv Shkedy** (ziv.shkedy@uhasselt.be)
