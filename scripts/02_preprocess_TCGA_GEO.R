# =============================================================
# 02_preprocess_TCGA_GEO.R  —  Κεφάλαιο 3
# Καθαρισμός/οργάνωση TCGA (mRNA, μεθυλίωση, CNA, μεταλλάξεις,
# κλινικά), GEO κοόρτες, META-cohort με ComBat, gene set dMMR/MSI
# =============================================================
source("scripts/00_config.R")
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(Biobase); library(sva)
  library(ggplot2); library(maftools)
})

# -------------------------------------------------------------
# 3.2.1 TCGA mRNA (RSEM) -> log2(x+1), μόνο πρωτογενείς όγκοι
# -------------------------------------------------------------
rsem <- fread("data/tcga/data_mrna_seq_v2_rsem.txt", data.table = FALSE)
rsem <- rsem[rsem$Hugo_Symbol != "" & !is.na(rsem$Hugo_Symbol), ]
m <- as.matrix(rsem[, -(1:2)]); rownames(m) <- NULL
# διπλά σύμβολα γονιδίων: κρατάμε τη γραμμή με τη μεγαλύτερη μέση έκφραση
keep <- !duplicated(rsem$Hugo_Symbol[order(-rowMeans(m, na.rm = TRUE))])
ord  <- order(-rowMeans(m, na.rm = TRUE))
m <- m[ord, ][keep, ]; rownames(m) <- rsem$Hugo_Symbol[ord][keep]
m[is.na(m)] <- 0
m <- m[, grepl("-01$", colnames(m))]                    # primary tumour
colnames(m) <- substr(colnames(m), 1, 12)               # patient barcode
expr_tcga <- log2(m + 1)
# φιλτράρισμα χαμηλής έκφρασης: έκφραση >1 (log2) σε >=10% δειγμάτων
expr_tcga <- expr_tcga[rowMeans(expr_tcga > 1) >= 0.10, ]
cat("TCGA mRNA:", dim(expr_tcga), "\n")

# φυσιολογικοί ιστοί (για σύγκριση όγκος vs φυσιολογικό)
nrm <- fread("data/tcga/normals_rsem.txt", data.table = FALSE)
nrm <- nrm[nrm$Hugo_Symbol %in% gene_set & !duplicated(nrm$Hugo_Symbol), ]
normal_expr <- log2(as.matrix(nrm[, -(1:2)]) + 1); rownames(normal_expr) <- nrm$Hugo_Symbol

# -------------------------------------------------------------
# 3.2.6 Κλινικά & επιβίωση
# -------------------------------------------------------------
pat <- fread("data/tcga/data_clinical_patient.txt", skip = 4, data.table = FALSE)
smp <- fread("data/tcga/data_clinical_sample.txt",  skip = 4, data.table = FALSE)
smp <- smp[grepl("-01$", smp$SAMPLE_ID), ]
num <- function(x) suppressWarnings(as.numeric(x))
st  <- function(x) as.integer(substr(x, 1, 1))

clin_tcga <- pat %>%
  transmute(
    ID       = PATIENT_ID,
    project  = CANCER_TYPE_ACRONYM,
    tcga_subtype = sub("^(COAD|READ)_", "", SUBTYPE),
    age      = num(AGE),
    gender   = SEX,
    stage_raw = AJCC_PATHOLOGIC_TUMOR_STAGE,
    stage    = case_when(grepl("STAGE IV", stage_raw) ~ "IV",
                         grepl("STAGE III", stage_raw) ~ "III",
                         grepl("STAGE II", stage_raw) ~ "II",
                         grepl("STAGE I", stage_raw) ~ "I", TRUE ~ NA_character_),
    T = sub("^(T[0-4]).*", "\\1", PATH_T_STAGE),
    N = sub("^(N[0-3]).*", "\\1", PATH_N_STAGE),
    M = sub("^(M[01]).*", "\\1", PATH_M_STAGE),
    OS.time  = num(OS_MONTHS) * 30.44,  OS  = st(OS_STATUS),
    DSS.time = num(DSS_MONTHS) * 30.44, DSS = st(DSS_STATUS),
    PFI.time = num(PFS_MONTHS) * 30.44, PFI = st(PFS_STATUS)
  ) %>%
  left_join(smp %>% transmute(ID = PATIENT_ID, site = TUMOR_TISSUE_SITE,
                              MSI_MANTIS = num(MSI_SCORE_MANTIS),
                              MSIsensor = num(MSI_SENSOR_SCORE),
                              TMB = num(TMB_NONSYNONYMOUS),
                              aneuploidy = num(ANEUPLOIDY_SCORE)), by = "ID") %>%
  mutate(
    # MSI-H κατά MANTIS (κατώφλι 0.4, Bonneville et al. 2017)
    MSI_status = case_when(MSI_MANTIS >= 0.4 ~ "MSI-H", MSI_MANTIS < 0.4 ~ "MSS/MSI-L", TRUE ~ NA_character_),
    T = ifelse(T %in% c("T1","T2","T3","T4"), T, NA),
    N = ifelse(N %in% c("N0","N1","N2"), N, NA),
    M = ifelse(M %in% c("M0","M1"), M, NA)
  )

# Δείγματα με έκφραση + επιβίωση (OS >= 30 ημέρες, όπως στη μελέτη-αναφορά)
clin_tcga <- clin_tcga %>% filter(ID %in% colnames(expr_tcga), !is.na(OS.time), !is.na(OS), OS.time >= 30)
expr_tcga <- expr_tcga[, clin_tcga$ID]
rownames(clin_tcga) <- clin_tcga$ID
cat("TCGA ασθενείς με mRNA + OS:", nrow(clin_tcga), "\n")
print(table(clin_tcga$project)); print(table(clin_tcga$MSI_status, useNA = "ifany"))

# -------------------------------------------------------------
# 3.2.4 Μεθυλίωση DNA: promoter probes (TSS1500/TSS200/5'UTR/1stExon)
#       των γονιδίων dMMR/MSI -> μέσο β ανά γονίδιο
# -------------------------------------------------------------
meth <- fread("data/tcga/data_methylation_hm450.txt", data.table = FALSE)
meth <- meth[meth$NAME %in% gene_set & grepl("TSS1500|TSS200|5'UTR|1stExon", meth$DESCRIPTION), ]
mm <- as.matrix(meth[, -(1:4)])
mm <- mm[, grepl("-01$", colnames(mm))]; colnames(mm) <- substr(colnames(mm), 1, 12)
meth_probe <- mm; rownames(meth_probe) <- meth$ENTITY_STABLE_ID
meth_probe_gene <- setNames(meth$NAME, meth$ENTITY_STABLE_ID)
meth_gene <- t(sapply(split(seq_len(nrow(mm)), meth$NAME), function(i) colMeans(mm[i, , drop = FALSE], na.rm = TRUE)))
meth_gene <- meth_gene[, colSums(is.na(meth_gene)) == 0, drop = FALSE]
cat("Μεθυλίωση (γονίδια x δείγματα):", dim(meth_gene), "\n")
rm(meth, mm); gc()

# -------------------------------------------------------------
# CNA (log2 ratio ανά γονίδιο)
# -------------------------------------------------------------
cna <- fread("data/tcga/data_log2_cna.txt", data.table = FALSE)
cna <- cna[cna$Hugo_Symbol %in% gene_set & !duplicated(cna$Hugo_Symbol), ]
cna_gene <- as.matrix(cna[, -(1:2)]); rownames(cna_gene) <- cna$Hugo_Symbol
cna_gene <- cna_gene[, grepl("-01$", colnames(cna_gene))]; colnames(cna_gene) <- substr(colnames(cna_gene), 1, 12)

# -------------------------------------------------------------
# 3.2.5 Σωματικές μεταλλάξεις (TCGA MC3)
# -------------------------------------------------------------
maf_coad <- readRDS("data/mut/COAD.RDs"); maf_read <- readRDS("data/mut/READ.RDs")
maf_all  <- merge_mafs(list(maf_coad, maf_read), verbose = FALSE)
nonsyn <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins",
            "In_Frame_Del","In_Frame_Ins","Splice_Site","Translation_Start_Site","Nonstop_Mutation")
mdat <- as.data.frame(maf_all@data)
mdat$ID <- substr(mdat$Tumor_Sample_Barcode, 1, 12)
mut_ids <- unique(substr(as.character(maf_all@clinical.data$Tumor_Sample_Barcode), 1, 12))
mut_bin <- sapply(gene_set, function(g) as.integer(mut_ids %in% mdat$ID[mdat$Hugo_Symbol == g & mdat$Variant_Classification %in% nonsyn]))
mut_bin <- t(mut_bin); colnames(mut_bin) <- mut_ids
cat("Δείγματα με μεταλλάξεις:", length(mut_ids), "\n")
print(sort(rowMeans(mut_bin), decreasing = TRUE))

# -------------------------------------------------------------
# 3.3 GEO κοόρτες (curatedCRCData)
# -------------------------------------------------------------
load_eset <- function(g) { e <- new.env(); load(sprintf("data/geo/%s_eset.rda", g), envir = e); get(ls(e)[1], e) }
geo_ids <- c("GSE39582", "GSE17536", "GSE17537")
# Σημ.: το GSE14333 ΔΕΝ μπαίνει στο META-cohort, επειδή περιέχει δείγματα
# του Moffitt Cancer Center που επικαλύπτονται με το GSE17536.
geo <- lapply(setNames(geo_ids, geo_ids), function(g) {
  x <- load_eset(g); p <- pData(x)
  d <- data.frame(
    ID = paste0(g, "_", colnames(x)), dataset = g,
    age = as.numeric(p$age_at_initial_pathologic_diagnosis),
    gender = ifelse(p$gender == "m", "Male", ifelse(p$gender == "f", "Female", NA)),
    stage = c(A = "I", B = "II", C = "III", D = "IV")[as.character(p$Dstage)],
    MSI_status = ifelse(p$msi == "MSI", "MSI-H", ifelse(p$msi == "MSS", "MSS/MSI-L", NA)),
    OS.time = as.numeric(p$days_to_death),
    OS = ifelse(p$vital_status == "deceased", 1, ifelse(p$vital_status == "living", 0, NA)),
    RFS.time = as.numeric(p$days_to_recurrence_or_death),
    RFS = ifelse(p$dfs_status == "deceased_or_recurrence", 1, ifelse(p$dfs_status == "living_norecurrence", 0, NA))
  )
  ex <- exprs(x); colnames(ex) <- d$ID
  ex <- ex[rowSums(is.na(ex)) == 0, ]
  list(expr = ex, clin = d)
})
for (g in geo_ids) cat(g, ": n =", nrow(geo[[g]]$clin), " OS events =", sum(geo[[g]]$clin$OS, na.rm = TRUE),
                       " RFS events =", sum(geo[[g]]$clin$RFS, na.rm = TRUE), "\n")

# -------------------------------------------------------------
# 3.4 META-cohort: κοινά γονίδια + ComBat (batch = dataset)
# -------------------------------------------------------------
common <- Reduce(intersect, lapply(geo, function(z) rownames(z$expr)))
meta_raw  <- do.call(cbind, lapply(geo, function(z) z$expr[common, ]))
meta_clin <- do.call(rbind, lapply(geo, function(z) z$clin)); rownames(meta_clin) <- meta_clin$ID
batch <- meta_clin$dataset
meta_expr <- ComBat(dat = meta_raw, batch = batch, par.prior = TRUE)
cat("META-cohort:", dim(meta_expr), "\n")

# PCA πριν / μετά το ComBat (Εικ. 3.x)
pca_plot <- function(mat, title) {
  v <- apply(mat, 1, var); top <- names(sort(v, decreasing = TRUE))[1:2000]
  pc <- prcomp(t(mat[top, ]), scale. = TRUE)
  ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  ggplot(data.frame(pc$x[, 1:2], dataset = batch), aes(PC1, PC2, colour = dataset)) +
    geom_point(size = 1.4, alpha = .75) + theme_bw(base_size = 12) +
    labs(title = title, x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)")) +
    scale_colour_manual(values = c("#2E6FA7", "#D1495B", "#EDAE49"))
}
p <- cowplot::plot_grid(pca_plot(meta_raw, "Before ComBat"), pca_plot(meta_expr, "After ComBat"), ncol = 2)
ggsave2(p, "Fig3_1_META_PCA_ComBat", w = 12, h = 5)

# -------------------------------------------------------------
# 3.5 dMMR/MSI gene set: παρουσία & έκφραση όγκος vs φυσιολογικό
# -------------------------------------------------------------
gs_tab <- data.frame(gene = gene_set,
                     in_TCGA = gene_set %in% rownames(expr_tcga),
                     in_META = gene_set %in% rownames(meta_expr),
                     in_methylation = gene_set %in% rownames(meth_gene),
                     mut_freq_TCGA = round(100 * rowMeans(mut_bin[gene_set, ]), 1))
write_tab(gs_tab, "Table3_4_geneset_coverage")
print(gs_tab)

gg <- intersect(gene_set, rownames(normal_expr))
long <- rbind(
  data.frame(gene = rep(gg, ncol(expr_tcga)), value = as.vector(expr_tcga[gg, ]), type = "Tumour"),
  data.frame(gene = rep(gg, ncol(normal_expr)), value = as.vector(normal_expr[gg, ]), type = "Normal"))
p <- ggplot(long, aes(gene, value, fill = type)) + geom_boxplot(outlier.size = .3) +
  ggpubr::stat_compare_means(aes(group = type), label = "p.signif", method = "wilcox.test") +
  scale_fill_manual(values = c(Normal = "#66A182", Tumour = "#D1495B")) +
  theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "log2(RSEM + 1)", fill = NULL, title = "dMMR/MSI genes: tumour vs normal tissue (TCGA)")
ggsave2(p, "Fig3_2_geneset_tumour_vs_normal", w = 10, h = 5)

# Πίνακας 3.x: περιγραφικά κοορτών
desc <- rbind(
  data.frame(cohort = "TCGA-COAD/READ", n = nrow(clin_tcga), OS_events = sum(clin_tcga$OS),
             MSI_H = sum(clin_tcga$MSI_status == "MSI-H", na.rm = TRUE),
             median_age = median(clin_tcga$age, na.rm = TRUE)),
  do.call(rbind, lapply(geo_ids, function(g) { d <- geo[[g]]$clin
    data.frame(cohort = g, n = nrow(d), OS_events = sum(d$OS, na.rm = TRUE),
               MSI_H = sum(d$MSI_status == "MSI-H", na.rm = TRUE), median_age = median(d$age, na.rm = TRUE)) })))
desc$RFS_events <- c(NA, sapply(geo_ids, function(g) sum(geo[[g]]$clin$RFS, na.rm = TRUE)))
write_tab(desc, "Table3_1_cohorts")
print(desc)

# ---- Αποθήκευση ----
saveRDS(list(expr = expr_tcga, clin = clin_tcga, meth_gene = meth_gene, meth_probe = meth_probe,
             meth_probe_gene = meth_probe_gene, cna = cna_gene, mut = mut_bin, normal = normal_expr),
        "data/processed/tcga.rds")
saveRDS(maf_all, "data/processed/tcga_maf.rds")
saveRDS(list(expr = meta_expr, clin = meta_clin, raw = meta_raw), "data/processed/meta.rds")
cat("\n=== Κεφ. 3 ολοκληρώθηκε ===\n")
