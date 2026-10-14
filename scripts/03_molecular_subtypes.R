# =============================================================
# 03_molecular_subtypes.R  —  Κεφάλαιο 4
# Multi-omics consensus clustering (MOVICS) με βάση το dMMR/MSI
# gene set, επιβίωση υποτύπων, εξωτερική επικύρωση (NTP, PAM)
# =============================================================
source("scripts/00_config.R")
suppressPackageStartupMessages({
  library(MOVICS); library(survival); library(survminer); library(dplyr)
  library(ggplot2); library(CMScaller)
})
tcga <- readRDS("data/processed/tcga.rds")
meta <- readRDS("data/processed/meta.rds")
FP <- "results/figures"   # φάκελος εικόνων για το MOVICS

# -------------------------------------------------------------
# 4.2.1 Multi-omics είσοδοι & ευθυγράμμιση δειγμάτων
# -------------------------------------------------------------
ids <- Reduce(intersect, list(tcga$clin$ID, colnames(tcga$meth_gene), colnames(tcga$mut), colnames(tcga$cna)))
mp  <- tcga$meth_probe[, ids]; mp <- mp[rowSums(is.na(mp)) == 0, ]
mut <- tcga$mut[, ids]; mut <- mut[rowSums(mut) >= 3, ]          # μεταλλάξεις σε >=3 όγκους
mo_data <- list(mRNA = tcga$expr[gene_set, ids],
                Methylation = mp,
                CNA = tcga$cna[gene_set, ids],
                Mutation = mut)
cat("Δείγματα με πλήρη multi-omics:", length(ids), "\n"); print(sapply(mo_data, nrow))
surv_mo <- tcga$clin[ids, c("OS.time", "OS")]; colnames(surv_mo) <- c("futime", "fustat")

# -------------------------------------------------------------
# 4.2.2 Προγνωστική αξιολόγηση χαρακτηριστικών (μονοπαραγοντική Cox)
#   Επειδή το gene set είναι προκαθορισμένο (15 γονίδια), ΔΕΝ πετάμε
#   χαρακτηριστικά· αναφέρουμε όμως την προγνωστική τους αξία.
# -------------------------------------------------------------
cox_one <- function(x, s) {
  if (length(unique(x)) < 2) return(c(HR = NA, lower = NA, upper = NA, p = NA))
  f <- summary(coxph(Surv(s$futime, s$fustat) ~ x))
  c(HR = f$conf.int[1], lower = f$conf.int[3], upper = f$conf.int[4], p = f$coefficients[5])
}
feat_cox <- do.call(rbind, lapply(names(mo_data), function(L) {
  X <- mo_data[[L]]
  r <- t(apply(X, 1, function(x) cox_one(if (L == "Mutation") x else as.numeric(scale(x)), surv_mo)))
  data.frame(layer = L, feature = rownames(X), round(r, 4))
}))
feat_cox$gene <- ifelse(feat_cox$layer == "Methylation", tcga$meth_probe_gene[feat_cox$feature], feat_cox$feature)
write_tab(feat_cox, "Table4_1_feature_univariate_cox")
cat("Χαρακτηριστικά με p<0.05:", sum(feat_cox$p < 0.05, na.rm = TRUE), "\n")

# -------------------------------------------------------------
# 4.2.3 Βέλτιστος αριθμός clusters (CPI + Gap statistic)
# -------------------------------------------------------------
optk <- getClustNum(data = mo_data, is.binary = c(FALSE, FALSE, FALSE, TRUE),
                    try.N.clust = 2:8, fig.path = FP, fig.name = "Fig4_1_optimal_k")
K <- 2   # CPI/Gap -> 2· συμφωνεί με τη βιολογία dMMR (MSI-H vs MSS)

# -------------------------------------------------------------
# 4.2.4 Ensemble clustering με 10 αλγορίθμους + consensus
# -------------------------------------------------------------
algs <- c("SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering",
          "IntNMF", "CIMLR", "MoCluster", "iClusterBayes")
moic_file <- "data/processed/moic_list.rds"   # αν υπάρχει ήδη, δεν ξανατρέχει (~6 λεπτά)
moic_list <- if (file.exists(moic_file)) readRDS(moic_file) else list()
for (a in setdiff(algs, names(moic_list))) {
  cat("\n>>", a, "\n")
  r <- tryCatch(getMOIC(data = mo_data, methodslist = a, N.clust = K,
                        type = c("gaussian", "gaussian", "gaussian", "binomial")),
                error = function(e) { message("  παράλειψη ", a, ": ", conditionMessage(e)); NULL })
  if (!is.null(r)) moic_list[[a]] <- r
}
cat("Αλγόριθμοι που έτρεξαν:", names(moic_list), "\n")
saveRDS(moic_list, "data/processed/moic_list.rds")

cmoic <- getConsensusMOIC(moic.res.list = moic_list, fig.name = "Fig4_2_consensus_heatmap",
                          fig.path = FP, distance = "euclidean", linkage = "average")
getSilhouette(sil = cmoic$sil, fig.path = FP, fig.name = "Fig4_3_silhouette", height = 5.5, width = 5)

# Ονομασία: CS1 = υποτύπος με τη χαμηλότερη έκφραση MLH1 (dMMR-like)
cl <- cmoic$clust.res
mlh1_by <- tapply(tcga$expr["MLH1", cl$samID], cl$clust, median)
lab <- setNames(paste0("CS", rank(mlh1_by)), names(mlh1_by))
cl$Subtype <- lab[as.character(cl$clust)]
cl$clust <- as.integer(sub("CS", "", cl$Subtype))
cmoic$clust.res <- cl
print(table(cl$Subtype))

# Multi-omics heatmap
plotdata <- getStdiz(data = mo_data, halfwidth = c(2, 2, 2, NA), centerFlag = c(TRUE, TRUE, TRUE, FALSE),
                     scaleFlag = c(TRUE, TRUE, TRUE, FALSE))
ann <- tcga$clin[ids, c("MSI_status", "stage")]
ann$MSI_status[is.na(ann$MSI_status)] <- "NA"; ann$stage[is.na(ann$stage)] <- "NA"
ann_col <- list(MSI_status = c("MSI-H" = "#7B3294", "MSS/MSI-L" = "#C2A5CF", "NA" = "grey85"),
                stage = c(I = "#FEE5D9", II = "#FCAE91", III = "#FB6A4A", IV = "#CB181D", "NA" = "grey85"))
col_list <- list(c("#00FF00", "#008000", "#000000", "#800000", "#FF0000"),
                 c("#6699CC", "white", "#FF3C38"),
                 c("#2166AC", "white", "#B2182B"),
                 c("grey90", "black"))
getMoHeatmap(data = plotdata, row.title = names(mo_data), is.binary = c(FALSE, FALSE, FALSE, TRUE),
             legend.name = c("mRNA z", "beta z", "CNA z", "Mut"), color = col_list,
             clust.res = cl, clust.col = unname(COL_CS[1:K]), annCol = ann, annColors = ann_col,
             width = 10, height = 8, fig.path = FP, fig.name = "Fig4_4_multiomics_heatmap")

# -------------------------------------------------------------
# 4.2.5 Διαφορές επιβίωσης (OS, PFI) μεταξύ υποτύπων
# -------------------------------------------------------------
km_plot <- function(df, time, event, group, title, pal, name, xlab = "Time (years)") {
  df$tt <- df[[time]] / 365.25; df$ee <- df[[event]]; df$gg <- factor(df[[group]])
  df <- df[is.finite(df$tt) & !is.na(df$ee) & !is.na(df$gg), ]
  fit <- survfit(Surv(tt, ee) ~ gg, data = df)
  names(fit$strata) <- sub("gg=", "", names(fit$strata))
  p <- ggsurvplot(fit, data = df, pval = TRUE, risk.table = TRUE, palette = unname(pal[levels(df$gg)]),
                  xlab = xlab, ylab = "Survival probability", title = title,
                  legend.labs = levels(df$gg), legend.title = "", risk.table.height = .25,
                  ggtheme = theme_bw(base_size = 12))
  pdf(fig(paste0(name, ".pdf")), width = 6.5, height = 6.5, onefile = FALSE); print(p); dev.off()
  png(fig(paste0(name, ".png")), width = 6.5, height = 6.5, units = "in", res = 300); print(p); dev.off()
  sd <- survdiff(Surv(tt, ee) ~ gg, data = df)
  hr <- if (nlevels(df$gg) == 2) summary(coxph(Surv(tt, ee) ~ gg, data = df))$conf.int[1, c(1, 3, 4)] else c(NA, NA, NA)
  data.frame(analysis = name, n = nrow(df), events = sum(df$ee), logrank_p = 1 - pchisq(sd$chisq, length(sd$n) - 1),
             HR = hr[1], HR_low = hr[2], HR_up = hr[3])
}
tc <- tcga$clin[ids, ]; tc$Subtype <- cl$Subtype[match(ids, cl$samID)]
km_res <- rbind(
  km_plot(tc, "OS.time", "OS", "Subtype", "TCGA — Overall survival (OS)", COL_CS, "Fig4_5a_KM_TCGA_OS"),
  km_plot(tc, "PFI.time", "PFI", "Subtype", "TCGA — Progression-free interval (PFI)", COL_CS, "Fig4_5b_KM_TCGA_PFI"),
  km_plot(tc, "DSS.time", "DSS", "Subtype", "TCGA — Disease-specific survival (DSS)", COL_CS, "Fig4_5c_KM_TCGA_DSS"))

# Κλινικά χαρακτηριστικά ανά υποτύπο (+ CMS ταξινόμηση)
cms <- CMScaller(tcga$expr[, ids] , RNAseq = FALSE, rowNames = "symbol", doPlot = FALSE, seed = 1)
tc$CMS <- as.character(cms$prediction); tc$CMS[is.na(tc$CMS)] <- "NA"
tc$MSI_MANTIS_cont <- tc$MSI_MANTIS
clin_tab <- tc %>% group_by(Subtype) %>% summarise(
  n = n(), age_median = median(age, na.rm = TRUE), female_pct = round(100 * mean(gender == "Female"), 1),
  stageIII_IV_pct = round(100 * mean(stage %in% c("III", "IV")), 1),
  MSI_H_pct = round(100 * mean(MSI_status == "MSI-H", na.rm = TRUE), 1),
  right_colon_pct = NA, TMB_median = median(TMB, na.rm = TRUE),
  MANTIS_median = median(MSI_MANTIS, na.rm = TRUE), rectum_pct = round(100 * mean(project == "READ"), 1),
  CMS1_pct = round(100 * mean(CMS == "CMS1"), 1), CMS2_pct = round(100 * mean(CMS == "CMS2"), 1),
  CMS3_pct = round(100 * mean(CMS == "CMS3"), 1), CMS4_pct = round(100 * mean(CMS == "CMS4"), 1))
clin_tab$right_colon_pct <- NULL
pv <- c(MSI = fisher.test(table(tc$Subtype, tc$MSI_status))$p.value,
        stage = fisher.test(table(tc$Subtype, tc$stage %in% c("III", "IV")))$p.value,
        CMS = fisher.test(table(tc$Subtype, tc$CMS), simulate.p.value = TRUE, B = 1e4)$p.value,
        TCGA_subtype = fisher.test(table(tc$Subtype, tc$tcga_subtype), simulate.p.value = TRUE, B = 1e4)$p.value,
        TMB = wilcox.test(TMB ~ Subtype, data = tc)$p.value,
        age = wilcox.test(age ~ Subtype, data = tc)$p.value)
write_tab(clin_tab, "Table4_2_subtype_clinical"); write_tab(data.frame(test = names(pv), p = signif(pv, 3)), "Table4_2b_subtype_clinical_pvalues")
print(clin_tab); print(signif(pv, 3))

# Ραβδόγραμμα MSI / CMS / TCGA subtype ανά CS
bar_df <- rbind(
  data.frame(Subtype = tc$Subtype, var = "MSI status", level = tc$MSI_status),
  data.frame(Subtype = tc$Subtype, var = "CMS", level = tc$CMS),
  data.frame(Subtype = tc$Subtype, var = "TCGA subtype", level = tc$tcga_subtype))
bar_df <- bar_df[!is.na(bar_df$level) & bar_df$level != "", ]
p <- ggplot(bar_df, aes(Subtype, fill = level)) + geom_bar(position = "fill") + facet_wrap(~var) +
  theme_bw(base_size = 12) + labs(y = "Proportion", x = NULL, fill = NULL) +
  scale_fill_manual(values = c("#7B3294","#C2A5CF","#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","grey80","#1F78B4"))
ggsave2(p, "Fig4_6_subtype_MSI_CMS", w = 10, h = 4.5)

# -------------------------------------------------------------
# Επέκταση ετικετών σε ΟΛΟΥΣ τους TCGA ασθενείς (n=560) με NTP
# -------------------------------------------------------------
# 4.2.6 Εξωτερική επικύρωση: templates από mRNA (όλο το μεταγράφωμα)
runDEA(dea.method = "limma", expr = tcga$expr[, ids], moic.res = cmoic, prefix = "TCGA",
       overwt = TRUE, res.path = "data/processed")
mk <- runMarker(moic.res = cmoic, dea.method = "limma", prefix = "TCGA",
                dat.path = "data/processed", res.path = "results/tables",
                p.cutoff = 0.05, p.adj.cutoff = 0.05, dirct = "up", n.marker = 100, doplot = FALSE)
# heatmap των marker γονιδίων (Εικ. 4.7)
mg <- mk$templates$probe
hm <- t(scale(t(tcga$expr[mg, cl$samID[order(cl$Subtype)]]))); hm[hm > 2] <- 2; hm[hm < -2] <- -2
ha <- ComplexHeatmap::HeatmapAnnotation(Subtype = sort(cl$Subtype),
        MSI = tcga$clin[cl$samID[order(cl$Subtype)], "MSI_status"],
        col = list(Subtype = COL_CS[1:K], MSI = c("MSI-H" = "#7B3294", "MSS/MSI-L" = "#C2A5CF")), na_col = "grey85")
save_pdf_png(function() ComplexHeatmap::draw(ComplexHeatmap::Heatmap(hm, name = "z", top_annotation = ha,
  cluster_columns = FALSE, show_column_names = FALSE, row_names_gp = grid::gpar(fontsize = 4),
  row_split = mk$templates$class, col = circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")))),
  "Fig4_7_marker_heatmap", w = 8, h = 9)
tmpl <- mk$templates
# όλοι οι TCGA
ntp_tcga <- runNTP(expr = tcga$expr, templates = tmpl, scaleFlag = TRUE, centerFlag = TRUE, doPlot = FALSE)
all_sub <- data.frame(ID = rownames(ntp_tcga$clust.res), Subtype = paste0("CS", ntp_tcga$clust.res$clust))
all_sub$Subtype[match(cl$samID, all_sub$ID)] <- cl$Subtype   # οι 350 κρατούν την ετικέτα του consensus
# META-cohort: NTP
ntp_meta <- runNTP(expr = meta$expr, templates = tmpl, scaleFlag = TRUE, centerFlag = TRUE,
                   doPlot = TRUE, fig.path = FP, fig.name = "Fig4_8_NTP_META")
mc <- meta$clin; mc$Subtype <- paste0("CS", ntp_meta$clust.res[mc$ID, "clust"])
# PAM (Kappa συμφωνίας)
pam_meta <- runPAM(train.expr = tcga$expr[, ids], moic.res = cmoic, test.expr = meta$expr)
mc$Subtype_PAM <- paste0("CS", pam_meta$clust.res[mc$ID, "clust"])
kap <- runKappa(subt1 = ntp_meta$clust.res$clust, subt2 = pam_meta$clust.res$clust,
                subt1.lab = "NTP", subt2.lab = "PAM", fig.path = FP, fig.name = "Fig4_9_kappa_NTP_PAM")

km_res <- rbind(km_res,
  km_plot(mc, "RFS.time", "RFS", "Subtype", "META-cohort — Relapse-free survival (RFS)", COL_CS, "Fig4_10a_KM_META_RFS"),
  km_plot(mc[!is.na(mc$OS), ], "OS.time", "OS", "Subtype", "META (GSE17536+GSE17537) — OS", COL_CS, "Fig4_10b_KM_META_OS"),
  km_plot(mc[mc$dataset == "GSE39582", ], "RFS.time", "RFS", "Subtype", "GSE39582 — RFS", COL_CS, "Fig4_10c_KM_GSE39582_RFS"))
write_tab(km_res, "Table4_3_subtype_survival")
print(km_res)
# MSI συμφωνία στο GSE39582
msi_tab <- table(mc$Subtype, mc$MSI_status); print(msi_tab)
write_tab(as.data.frame.matrix(msi_tab) %>% tibble::rownames_to_column("Subtype"), "Table4_4_GSE39582_subtype_vs_MSI")
cat("Fisher p (GSE39582 CS vs MSI):", fisher.test(msi_tab)$p.value, "\n")

# -------------------------------------------------------------
# 4.2.7 Απεικόνιση χαμηλής διάστασης (PCA, t-SNE, UMAP)
# -------------------------------------------------------------
low_dim <- function(X, grp, title) {
  X <- t(scale(t(X))); X[is.na(X)] <- 0
  pc <- prcomp(t(X)); set.seed(1)
  ts <- Rtsne::Rtsne(t(X), perplexity = 30, check_duplicates = FALSE)$Y
  um <- uwot::umap(t(X), n_neighbors = 15, min_dist = .3)
  d <- rbind(data.frame(x = pc$x[, 1], y = pc$x[, 2], m = "PCA"), data.frame(x = ts[, 1], y = ts[, 2], m = "t-SNE"),
             data.frame(x = um[, 1], y = um[, 2], m = "UMAP")); d$Subtype <- rep(grp, 3)
  ggplot(d, aes(x, y, colour = Subtype)) + geom_point(size = 1.2, alpha = .8) + facet_wrap(~m, scales = "free") +
    scale_colour_manual(values = COL_CS) + theme_bw(base_size = 12) + labs(x = NULL, y = NULL, title = title)
}
p1 <- low_dim(tcga$expr[gene_set, ids], cl$Subtype[match(ids, cl$samID)], "TCGA (dMMR/MSI genes)")
g2 <- intersect(gene_set, rownames(meta$expr))
p2 <- low_dim(meta$expr[g2, mc$ID], mc$Subtype, "META-cohort (dMMR/MSI genes)")
ggsave2(cowplot::plot_grid(p1, p2, ncol = 1), "Fig4_11_lowdim", w = 11, h = 8)

# ---- Αποθήκευση ----
saveRDS(list(cmoic = cmoic, ids = ids, tcga_all = all_sub, tcga_mo = tc, meta = mc, templates = tmpl, mo_data = mo_data),
        "data/processed/subtypes.rds")
cat("\n=== Κεφ. 4 ολοκληρώθηκε ===\n")
