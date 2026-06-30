# CLL Dataset

The chronic lymphocytic leukaemia (CLL) multi-omics dataset (Dietrich et al.,
2018) is **not stored as a file** — it ships with the **MOFAdata** Bioconductor
package and is loaded directly in R. This is exactly how `cll_study/app_to_cll_data.R`
obtains it.

- **Modalities used:** DNA methylation + drug response (n = 200 samples)
- **Package:** [MOFAdata](https://bioconductor.org/packages/MOFAdata)
- **Sample metadata:** EBI MOFA vignette FTP (see below)

## Load it

```r
# install.packages("BiocManager"); BiocManager::install("MOFAdata")
library(MOFAdata)
library(data.table)

data("CLL_data")                 # list: Drugs, Methylation, mRNA, Mutations
CLL_data2   <- CLL_data[c(1, 2)] # keep Drugs + Methylation (as used in the paper)
drugs       <- CLL_data2[[1]]    # 310 x 200
methylation <- CLL_data2[[2]]    # 4,248 x 200

# Sample annotation (IGHV status etc.) used in Figure 2
CLL_metadata <- fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/sample_metadata.txt")
```

Full clinical metadata is also available via the Bioconductor package
`BloodCancerMultiOmics2017` (`data("patmeta")`).
