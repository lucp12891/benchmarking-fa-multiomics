# scripts/01_prepare_data.R
# Load and preprocess the two real-world multi-omics datasets:
#   1. CLL  (DNA methylation + drug response)  -- via the MOFAdata package
#   2. Radiation (RNA-seq + proteomics)        -- from local files in data/radiation/
#
# Output: data/CLL/cll_omics.rds and data/radiation/radiation_omics.rds, each a
# named list of (features x samples) matrices on a common sample set.

source("R/utils/libs.R")
source("R/utils/preprocessing.R")

dir.create("data/CLL",       recursive = TRUE, showWarnings = FALSE)
dir.create("data/radiation", recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. CLL dataset (Dietrich et al. 2018), distributed with MOFAdata
# ---------------------------------------------------------------------------
prepare_cll <- function(top_n = c(mRNA = 5000, Methylation = 5000)) {
  if (!requireNamespace("MOFAdata", quietly = TRUE)) {
    stop("Install MOFAdata: BiocManager::install('MOFAdata').")
  }
  data("CLL_data", package = "MOFAdata", envir = environment())

  # CLL_data is a list of (features x samples) matrices: Drugs, Methylation, mRNA, Mutations
  omics <- list(
    Methylation = as.matrix(CLL_data$Methylation),
    Drugs       = as.matrix(CLL_data$Drugs)
  )
  omics <- preprocess_omics(omics, top_n = NULL, standardize = TRUE)
  saveRDS(omics, "data/CLL/cll_omics.rds")
  message("CLL data prepared: ",
          paste(names(omics), sapply(omics, nrow), "features", collapse = "; "))
  invisible(omics)
}

# ---------------------------------------------------------------------------
# 2. Radiation dataset (RNA-seq: PRJNA1431456; proteomics: PXD075411)
#    Expected raw inputs (see data/radiation/README.md):
#      data/radiation/rnaseq_counts.csv     (genes    x samples)
#      data/radiation/proteomics.csv        (proteins x samples)
# ---------------------------------------------------------------------------
prepare_radiation <- function(rna_path = "data/radiation/rnaseq_counts.csv",
                              pro_path = "data/radiation/proteomics.csv") {
  if (!file.exists(rna_path) || !file.exists(pro_path)) {
    warning("Radiation raw files not found. See data/radiation/README.md for downloads.")
    return(invisible(NULL))
  }
  rna <- as.matrix(read.csv(rna_path, row.names = 1, check.names = FALSE))
  pro <- as.matrix(read.csv(pro_path, row.names = 1, check.names = FALSE))

  # log-transform RNA-seq counts; keep proteomics on its reported scale
  rna <- log2(rna + 1)

  omics <- preprocess_omics(
    list(RNAseq = rna, Proteomics = pro),
    top_n       = c(RNAseq = 5000, Proteomics = NULL),
    standardize = TRUE
  )
  saveRDS(omics, "data/radiation/radiation_omics.rds")
  message("Radiation data prepared: ",
          paste(names(omics), sapply(omics, nrow), "features", collapse = "; "))
  invisible(omics)
}

if (!interactive()) {
  prepare_cll()
  prepare_radiation()
}
