# Data

Raw data are **not** committed to this repository. This folder documents where
to obtain each dataset and what the preprocessing scripts expect. Running
`scripts/01_prepare_data.R` populates the `.rds` files consumed downstream.

## CLL dataset (`CLL/`)

DNA methylation + drug response, n = 200 chronic lymphocytic leukaemia samples
(Dietrich et al., 2018). Distributed with the **MOFAdata** Bioconductor package —
no manual download needed.

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("MOFAdata")
```

`scripts/01_prepare_data.R` loads `CLL_data`, selects the Methylation and Drugs
views, preprocesses them, and writes `data/CLL/cll_omics.rds`.

## Radiation dataset (`radiation/`)

RNA-seq + proteomics, n = 16 samples. See `radiation/README.md` for the NCBI/PRIDE
accessions and the expected raw-file layout.

## Produced files

| File | Produced by | Contents |
|------|-------------|----------|
| `CLL/cll_omics.rds`             | `01_prepare_data.R` | list(Methylation, Drugs) — features x samples |
| `radiation/radiation_omics.rds` | `01_prepare_data.R` | list(RNAseq, Proteomics) — features x samples |
