# =============================================================
# 01_download_data.R  —  Κατέβασμα όλων των δεδομένων (Κεφ. 3)
# Τρέχει ΜΙΑ φορά. Χρειάζεται ~3 GB χώρο και καλή σύνδεση.
#
# Πηγές (όλες δημόσιες):
#  * TCGA COAD+READ PanCancer Atlas (mRNA RSEM, μεθυλίωση 450K, CNA,
#    κλινικά/επιβίωση, MSI scores, TMB)  -> cBioPortal DataHub
#  * Σωματικές μεταλλάξεις TCGA MC3 (COAD, READ) -> TCGAmutations (GitHub)
#  * GEO κοόρτες (GSE39582, GSE17536, GSE17537, GSE14333, GSE13294)
#    ήδη επιμελημένες (gene symbols, κλινικά) -> Bioconductor curatedCRCData
#  * Αναφορές: LM22 (CIBERSORT), IPS γονίδια, GDSC 2016 (φάρμακα)
# =============================================================
source("scripts/00_config.R")

dir.create("data/tcga", recursive = TRUE, showWarnings = FALSE)
dir.create("data/geo",  recursive = TRUE, showWarnings = FALSE)
dir.create("data/mut",  recursive = TRUE, showWarnings = FALSE)
dir.create("data/ref",  recursive = TRUE, showWarnings = FALSE)

dl <- function(url, dest) {
  if (!file.exists(dest)) {
    message("Κατέβασμα: ", basename(dest))
    download.file(url, dest, mode = "wb", quiet = FALSE)
  }
}

# ---- 1. TCGA (cBioPortal DataHub, Git LFS) ----
B <- "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/coadread_tcga_pan_can_atlas_2018"
for (f in c("data_clinical_patient.txt", "data_clinical_sample.txt",
            "data_mrna_seq_v2_rsem.txt", "data_log2_cna.txt",
            "data_methylation_hm450.txt"))          # το τελευταίο ~1.4 GB
  dl(file.path(B, f), file.path("data/tcga", f))
dl(file.path(B, "normals/data_mrna_seq_v2_rsem_normal_samples.txt"),
   "data/tcga/normals_rsem.txt")

# ---- 2. Μεταλλάξεις TCGA MC3 ----
M <- "https://github.com/PoisonAlien/TCGAmutations/raw/master/inst/extdata/MC3"
dl(file.path(M, "COAD.RDs"), "data/mut/COAD.RDs")
dl(file.path(M, "READ.RDs"), "data/mut/READ.RDs")

# ---- 3. GEO κοόρτες (curatedCRCData) ----
if (!requireNamespace("curatedCRCData", quietly = TRUE))
  BiocManager::install("curatedCRCData", update = FALSE, ask = FALSE)
library(curatedCRCData)
for (g in c("GSE39582", "GSE17536", "GSE17537", "GSE14333", "GSE13294")) {
  obj <- paste0(g, "_eset")
  data(list = obj, package = "curatedCRCData", envir = environment())
  assign(obj, get(obj)); save(list = obj, file = sprintf("data/geo/%s_eset.rda", g))
}

# ---- 4. Αρχεία αναφοράς ----
dl("https://raw.githubusercontent.com/SiYangming/CIBERSORT-DATA/main/LM22.txt", "data/ref/LM22.txt")
dl("https://raw.githubusercontent.com/icbi-lab/Immunophenogram/master/IPS_genes.txt", "data/ref/IPS_genes.txt")
# GDSC 2016 (έρχεται μαζί με το MOVICS — δεν χρειάζεται download)

cat("\n=== ΟΛΟΚΛΗΡΩΘΗΚΕ. Δες τους φακέλους data/tcga, data/geo, data/mut, data/ref ===\n")
