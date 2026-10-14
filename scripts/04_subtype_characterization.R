# =============================================================
# 04_subtype_characterization.R  —  Κεφάλαιο 5
# Μοριακός χαρακτηρισμός υποτύπων CS1 / CS2:
#  5.2 διαφορική έκφραση · 5.3 μεταγραφικοί παράγοντες & chromatin
#  remodeling · 5.4 immune checkpoints & μικροπεριβάλλον (ESTIMATE,
#  CIBERSORT) · 5.5 ssGSEA μονοπατιών KEGG
# =============================================================
source("scripts/00_config.R")
suppressPackageStartupMessages({
  library(limma); library(ggplot2); library(dplyr); library(tidyr); library(ggpubr)
  library(ComplexHeatmap); library(circlize); library(msigdbr); library(org.Hs.eg.db)
})
tcga <- readRDS("data/processed/tcga.rds")
sub  <- readRDS("data/processed/subtypes.rds")
cl   <- sub$cmoic$clust.res
ids  <- cl$samID; grp <- factor(cl$Subtype, levels = c("CS1", "CS2"))
X    <- tcga$expr[, ids]

# ---- βοηθητικό: boxplots ανά υποτύπο με Wilcoxon ----
box_by <- function(mat, group, title, ylab, pal = COL_CS, ncol = NULL, fname, w = 10, h = 5) {
  d <- data.frame(t(mat), group = group, check.names = FALSE) %>%
    pivot_longer(-group, names_to = "feature", values_to = "value")
  d$feature <- factor(d$feature, levels = rownames(mat))
  p <- ggplot(d, aes(feature, value, fill = group)) +
    geom_boxplot(outlier.size = .3, width = .7) +
    stat_compare_means(aes(group = group), label = "p.signif", method = "wilcox.test", hide.ns = FALSE) +
    scale_fill_manual(values = pal) + theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = NULL, y = ylab, fill = NULL, title = title)
  ggsave2(p, fname, w = w, h = h)
  # πίνακας p-values
  data.frame(feature = rownames(mat),
             median_CS1 = apply(mat[, group == "CS1", drop = FALSE], 1, median),
             median_CS2 = apply(mat[, group == "CS2", drop = FALSE], 1, median),
             p = apply(mat, 1, function(v) wilcox.test(v ~ group)$p.value)) %>%
    mutate(FDR = p.adjust(p, "BH"))
}

# -------------------------------------------------------------
# 5.2 Διαφορική έκφραση CS1 vs CS2 (limma)
# -------------------------------------------------------------
design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
fit <- eBayes(contrasts.fit(lmFit(X, design), makeContrasts(CS1 - CS2, levels = design)))
deg <- topTable(fit, number = Inf) %>% tibble::rownames_to_column("gene")
deg$sig <- with(deg, ifelse(adj.P.Val < 0.05 & logFC > 1, "Up in CS1",
                     ifelse(adj.P.Val < 0.05 & logFC < -1, "Up in CS2", "NS")))
write_tab(deg, "Table5_1_DEG_CS1_vs_CS2")
cat("DEGs (|log2FC|>1, FDR<0.05):"); print(table(deg$sig))

lab <- deg %>% filter(sig != "NS") %>% arrange(adj.P.Val) %>% head(20)
p <- ggplot(deg, aes(logFC, -log10(adj.P.Val), colour = sig)) + geom_point(size = .8, alpha = .7) +
  scale_colour_manual(values = c("Up in CS1" = COL_CS[["CS1"]], "Up in CS2" = COL_CS[["CS2"]], NS = "grey75")) +
  geom_vline(xintercept = c(-1, 1), lty = 2) + geom_hline(yintercept = -log10(.05), lty = 2) +
  ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 3, colour = "black", max.overlaps = 30) +
  theme_bw(base_size = 12) + labs(x = "log2 fold change (CS1 vs CS2)", y = "-log10(FDR)", colour = NULL,
                                  title = "Differentially expressed genes between subtypes")
ggsave2(p, "Fig5_1_volcano", w = 7.5, h = 6)

top <- deg %>% filter(sig != "NS") %>% group_by(sig) %>% slice_min(adj.P.Val, n = 25) %>% pull(gene)
o <- order(grp); hm <- t(scale(t(X[top, o]))); hm[hm > 2] <- 2; hm[hm < -2] <- -2
ha <- HeatmapAnnotation(Subtype = grp[o], MSI = tcga$clin[ids[o], "MSI_status"],
                        col = list(Subtype = COL_CS[1:2], MSI = c("MSI-H" = "#7B3294", "MSS/MSI-L" = "#C2A5CF")), na_col = "grey85")
save_pdf_png(function() draw(Heatmap(hm, name = "z", top_annotation = ha, cluster_columns = FALSE,
  show_column_names = FALSE, row_names_gp = grid::gpar(fontsize = 6),
  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")))), "Fig5_2_DEG_heatmap", w = 8, h = 8)

# -------------------------------------------------------------
# 5.3 Μεταγραφικοί παράγοντες & chromatin remodeling (GO)
# -------------------------------------------------------------
go_genes <- function(go) unique(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = go,
                         keytype = "GOALL", columns = "SYMBOL"))$SYMBOL)
tf_genes  <- go_genes("GO:0003700")    # DNA-binding transcription factor activity
chr_genes <- go_genes("GO:0006338")    # chromatin remodeling
tf_deg  <- deg %>% filter(gene %in% tf_genes,  adj.P.Val < 0.05) %>% arrange(desc(abs(logFC)))
chr_deg <- deg %>% filter(gene %in% chr_genes, adj.P.Val < 0.05) %>% arrange(desc(abs(logFC)))
write_tab(tf_deg,  "Table5_2_TF_DEG"); write_tab(chr_deg, "Table5_3_chromatin_DEG")
cat("TF DEGs (FDR<0.05):", nrow(tf_deg), "| chromatin remodeling DEGs:", nrow(chr_deg), "\n")
pl <- rbind(head(tf_deg, 20) %>% mutate(class = "Transcription factors"),
            head(chr_deg, 20) %>% mutate(class = "Chromatin remodeling"))
p <- ggplot(pl, aes(reorder(gene, logFC), logFC, fill = logFC > 0)) + geom_col() + coord_flip() +
  facet_wrap(~class, scales = "free_y") + scale_fill_manual(values = c(`TRUE` = COL_CS[["CS1"]], `FALSE` = COL_CS[["CS2"]]), guide = "none") +
  theme_bw(base_size = 11) + labs(x = NULL, y = "log2 fold change (CS1 vs CS2)",
                                  title = "Top differentially expressed TFs and chromatin regulators")
ggsave2(p, "Fig5_3_TF_chromatin", w = 10, h = 6)

# -------------------------------------------------------------
# 5.4 Immune checkpoints & μικροπεριβάλλον όγκου
# -------------------------------------------------------------
cp <- intersect(CHECKPOINTS, rownames(X))
cp_tab <- box_by(X[cp, ], grp, "Immune checkpoint expression", "log2(RSEM+1)", fname = "Fig5_4_checkpoints", w = 11)
write_tab(cp_tab, "Table5_4_checkpoints")

# ESTIMATE (stromal / immune / estimate score, tumour purity)
est_in <- tidyestimate::filter_common_genes(as.data.frame(X), id = "hgnc_symbol", tidy = FALSE, tell_missing = FALSE, find_alias = TRUE)
est <- tidyestimate::estimate_score(est_in, is_affymetrix = FALSE)
est_m <- t(as.matrix(est[, c("stromal", "immune", "estimate")])); colnames(est_m) <- est$sample
est_m <- est_m[, ids]; rownames(est_m) <- c("StromalScore", "ImmuneScore", "ESTIMATEScore")
est_tab <- box_by(est_m, grp, "ESTIMATE scores", "Score", fname = "Fig5_5_ESTIMATE", w = 6, h = 5)
write_tab(est_tab, "Table5_5_ESTIMATE")

# CIBERSORT (LM22, 22 τύποι ανοσοκυττάρων)
source("data/ref/CIBERSORT.R")
suppressPackageStartupMessages({ library(tibble); library(e1071) })
lm22 <- as.matrix(read.delim("data/ref/LM22.txt", row.names = 1, check.names = FALSE))
ciber <- CIBERSORT(sig_matrix = lm22, mixture_file = as.data.frame(2^X - 1), perm = 0, QN = FALSE)
ciber <- ciber[ids, colnames(lm22)]
cib_tab <- box_by(t(ciber), grp, "CIBERSORT immune cell fractions (LM22)", "Fraction", fname = "Fig5_6_CIBERSORT", w = 13, h = 5.5)
write_tab(cib_tab, "Table5_6_CIBERSORT")

# -------------------------------------------------------------
# 5.5 ssGSEA μονοπατιών KEGG
# -------------------------------------------------------------
kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
kegg_sets <- split(kegg$gene_symbol, kegg$gs_name)
kegg_sc <- run_ssgsea(X, kegg_sets)
rownames(kegg_sc) <- sub("^KEGG_", "", rownames(kegg_sc))
kfit <- eBayes(contrasts.fit(lmFit(kegg_sc, design), makeContrasts(CS1 - CS2, levels = design)))
kres <- topTable(kfit, number = Inf) %>% tibble::rownames_to_column("pathway")
write_tab(kres, "Table5_7_KEGG_ssGSEA")
kt <- kres %>% filter(adj.P.Val < 0.05) %>% arrange(t)
kt <- rbind(head(kt, 15), tail(kt, 15)) %>% distinct()
hm <- t(scale(t(kegg_sc[kt$pathway, o])))
hm[hm > 2] <- 2; hm[hm < -2] <- -2
save_pdf_png(function() draw(Heatmap(hm, name = "ssGSEA z", top_annotation = HeatmapAnnotation(Subtype = grp[o], col = list(Subtype = COL_CS[1:2])),
  cluster_columns = FALSE, show_column_names = FALSE, row_names_gp = grid::gpar(fontsize = 7),
  row_names_max_width = unit(9, "cm"), col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")))),
  "Fig5_7_KEGG_ssGSEA_heatmap", w = 10, h = 8)

saveRDS(list(deg = deg, estimate = est_m, cibersort = ciber, kegg = kegg_sc), "data/processed/characterization.rds")
cat("\n=== Κεφ. 5 ολοκληρώθηκε ===\n")
