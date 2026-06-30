# Radiation Exposure Dataset

Transcriptomics + proteomics from n = 16 mice (folic-acid supplementation study).
The normalized matrices used in the analysis are stored here directly.

## Accessions (raw data)

| Modality | Repository | Accession |
|----------|------------|-----------|
| Transcriptomics (RNA-seq) | NCBI SRA BioProject | [PRJNA1431456](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1431456) |
| Proteomics                | PRIDE               | [PXD075411](https://www.ebi.ac.uk/pride/archive/projects/PXD075411) |

## Files in this folder

| File | Description |
|------|-------------|
| `D17_mRNA_Cortex_normalized_tmm.csv`        | RNA-seq, cortex — TMM-normalized counts (genes × samples) |
| `D17_mRNA_Hippocampus_normalized_tmm.csv`   | RNA-seq, hippocampus — TMM-normalized counts |
| `D17_Protein_Cortex_normalized_quant.csv`      | Proteomics, cortex — normalized quant (proteins × samples) |
| `D17_Protein_Hippocampus_normalized_quant.csv` | Proteomics, hippocampus — normalized quant |
| `D17_cortex_metadata.csv`        | Sample annotation, cortex |
| `D17_hippocampus_metadata.csv`   | Sample annotation, hippocampus |
| `rad_metadata.csv`               | Combined sample metadata (treatment groups) |
| `submission/Model.organism.animal.1.0_*.xlsx` | Formatted data-submission tables (Cortex, Hippocampus) |

These feed `radiation_study/app_to_rad_data.R` (Figures 8–9, Table 3).
