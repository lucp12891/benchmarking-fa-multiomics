# Radiation Exposure Dataset

Multi-omics dataset of n = 16 samples profiled by transcriptomics and proteomics.

## Accessions

| Modality | Repository | Accession |
|----------|------------|-----------|
| Transcriptomics (RNA-seq) | NCBI SRA BioProject | [PRJNA1431456](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1431456) |
| Proteomics                | PRIDE               | [PXD075411](https://www.ebi.ac.uk/pride/archive/projects/PXD075411) |

## Expected raw files

After download and quantification, place the following matrices here (samples as
columns, features as rows; first column = feature id):

```
data/radiation/rnaseq_counts.csv   # genes    x samples (raw or expected counts)
data/radiation/proteomics.csv      # proteins x samples (normalized intensities)
```

Sample column names **must match** across the two files so the views can be
harmonized to a common sample set.

## Preprocessing

`scripts/01_prepare_data.R::prepare_radiation()` log2-transforms the RNA-seq
counts, keeps the top-variable genes, feature-standardizes both views, and writes
`radiation_omics.rds`.
