# =============================================================
# 00_setup_packages.R
# Εγκατάσταση όλων των πακέτων για το project dMMR/MSI - CRC
# Τρέξε το ΜΙΑ φορά. Μπορεί να πάρει 20-40 λεπτά.
# =============================================================

# --- 0. Ρυθμίσεις ---
options(timeout = 100000)              # αποφυγή timeout σε μεγάλα downloads
Sys.setenv(R_INSTALL_STAGED = FALSE)

# --- 1. BiocManager (διαχειριστής πακέτων Bioconductor) ---
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# --- 2. Πακέτα CRAN ---
cran_pkgs <- c(
  "tidyverse", "dplyr", "data.table", "ggplot2", "survival",
  "survminer", "glmnet", "randomForestSRC", "timeROC", "rms",
  "pheatmap", "remotes", "devtools"
)
new_cran <- cran_pkgs[!cran_pkgs %in% installed.packages()[, "Package"]]
if (length(new_cran)) install.packages(new_cran)

# --- 3. Πακέτα Bioconductor ---
bioc_pkgs <- c(
  "TCGAbiolinks", "GEOquery", "limma", "sva", "maftools",
  "GSVA", "clusterProfiler", "org.Hs.eg.db", "ComplexHeatmap",
  "SummarizedExperiment", "edgeR"
)
new_bioc <- bioc_pkgs[!bioc_pkgs %in% installed.packages()[, "Package"]]
if (length(new_bioc)) BiocManager::install(new_bioc, update = FALSE, ask = FALSE)

# --- 4. Πακέτα από GitHub (ο πυρήνας της ανάλυσης) ---
# MOVICS: multi-omics consensus clustering (Κεφ. 4-5)
if (!requireNamespace("MOVICS", quietly = TRUE))
  remotes::install_github("xlucpu/MOVICS")

# Mime: benchmark ML προγνωστικών μοντέλων (Κεφ. 6)
if (!requireNamespace("Mime1", quietly = TRUE))
  remotes::install_github("l-magnificence/Mime")

# oncoPredict: πρόβλεψη ευαισθησίας φαρμάκων (Κεφ. 8)
if (!requireNamespace("oncoPredict", quietly = TRUE))
  BiocManager::install("oncoPredict", update = FALSE, ask = FALSE)

# estimate: stromal/immune scores (Κεφ. 5)
if (!requireNamespace("estimate", quietly = TRUE))
  install.packages("estimate", repos = "http://r-forge.r-project.org", dependencies = TRUE)

# --- 5. Έλεγχος ότι όλα φορτώνουν ---
check <- c("TCGAbiolinks", "GEOquery", "sva", "limma", "MOVICS",
           "survival", "glmnet", "maftools", "GSVA", "timeROC", "rms")
for (p in check) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("%-16s : %s\n", p, ifelse(ok, "OK", "ΛΕΙΠΕΙ - ξανατρέξε")))
}

cat("\nΈτοιμο. Αν κάποιο λέει 'ΛΕΙΠΕΙ', ξανατρέξε τη σχετική γραμμή.\n")
