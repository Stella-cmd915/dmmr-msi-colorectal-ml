# dMMR/MSI Colorectal Cancer — Multi-omics ML Pipeline

An R pipeline that downloads TCGA/GEO colorectal cancer data, preprocesses multi-omics layers and identifies dMMR/MSI-related molecular subtypes for prognosis and treatment-sensitivity analysis.

## Tech stack
- R · Bioconductor (TCGAbiolinks, maftools, survival, GSVA, …)

## Features
- Automated package setup
- TCGA/GEO download and preprocessing (expression, mutations, methylation)
- Consensus clustering of molecular subtypes
- Survival analysis, immune infiltration and pathway (ssGSEA) characterisation

## Getting started
```r
source("00_setup_packages.R")
Sys.setenv(PROJECT_DIR = "path/to/workdir")
source("scripts/01_download_data.R")  # then 02 → 04
```

---
Part of the **SIGMA LAB** project portfolio.
