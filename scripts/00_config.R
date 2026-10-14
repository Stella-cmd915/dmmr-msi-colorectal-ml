# =============================================================
# 00_config.R  —  Κοινές ρυθμίσεις για όλα τα scripts
# Το φορτώνει αυτόματα κάθε script με source("scripts/00_config.R")
# =============================================================

# ---- Φάκελος εργασίας (ΑΛΛΑΞΕ ΤΟ στον δικό σου υπολογιστή) ----
# π.χ. PROJECT_DIR <- "C:/Users/you/dMMR_MSI_project"
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = "~/dMMR_MSI_project")
setwd(PROJECT_DIR)

dir.create("data/processed",   recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables",   recursive = TRUE, showWarnings = FALSE)

set.seed(20260927)
# Το MOVICS αποθηκεύει εικόνες με dev.copy2pdf(), που θέλει "οθόνη".
# Στο RStudio αυτό ισχύει ήδη· σε Rscript/διακομιστή χρησιμοποιούμε X11.
if (!interactive() && nzchar(Sys.getenv("DISPLAY")))
  options(device = function(...) grDevices::x11(...))
options(stringsAsFactors = FALSE, timeout = 100000)

# ---- Χρώματα υποτύπων / ομάδων κινδύνου ----
COL_CS   <- c(CS1 = "#2E6FA7", CS2 = "#D1495B", CS3 = "#EDAE49", CS4 = "#66A182")
COL_RISK <- c(Low = "#2E6FA7", High = "#D1495B")

# ---- Βοηθητικές συναρτήσεις ----
fig  <- function(name) file.path("results/figures", name)
tab  <- function(name) file.path("results/tables",  name)
save_pdf_png <- function(plot_fun, name, w = 7, h = 6) {
  pdf(fig(paste0(name, ".pdf")), width = w, height = h); plot_fun(); dev.off()
  png(fig(paste0(name, ".png")), width = w, height = h, units = "in", res = 300); plot_fun(); dev.off()
  invisible(TRUE)
}
ggsave2 <- function(p, name, w = 7, h = 6) {
  ggplot2::ggsave(fig(paste0(name, ".pdf")), p, width = w, height = h)
  ggplot2::ggsave(fig(paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}
write_tab <- function(df, name) {
  write.csv(df, tab(paste0(name, ".csv")), row.names = FALSE)
}
# C-index (Harrell) — μεγαλύτερο σκορ = μεγαλύτερος κίνδυνος
cindex <- function(time, status, score) {
  ok <- is.finite(time) & is.finite(status) & is.finite(score)
  survival::concordance(survival::Surv(time[ok], status[ok]) ~ score[ok], reverse = TRUE)$concordance
}

gene_set <- trimws(readLines("data/dMMR_MSI_geneset.txt"))
gene_set <- gene_set[gene_set != ""]

# ssGSEA που δουλεύει και με παλιές (<1.50) και με νέες (>=1.52) εκδόσεις GSVA
run_ssgsea <- function(X, sets) {
  sets <- lapply(sets, function(s) intersect(s, rownames(X)))
  sets <- sets[lengths(sets) >= 3]
  if ("ssgseaParam" %in% getNamespaceExports("GSVA"))
    GSVA::gsva(GSVA::ssgseaParam(as.matrix(X), sets), verbose = FALSE)
  else GSVA::gsva(as.matrix(X), sets, method = "ssgsea", verbose = FALSE)
}
# Λίστα immune checkpoints (Κεφ. 5 και 8)
CHECKPOINTS <- c("CD274","PDCD1","PDCD1LG2","CTLA4","LAG3","HAVCR2","TIGIT","IDO1","CD80","CD86",
                 "ICOS","TNFRSF9","TNFRSF4","CD276","VSIR","BTLA","CD27","CD40","CD70","SIGLEC15")
