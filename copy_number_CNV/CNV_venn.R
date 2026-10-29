# CNV vs transcript (DEG) overlap Venn diagrams: DXB11, DG44, S vs K1.
# Recovered from the lab's R console history and cleaned (file.choose() removed, paths via config.R).
# Inputs (place in data/): CNV_DESeq2_Integrated_Analysis_Updated_latest.xlsx,
#                          Cell_line_specific_merged_K1.xlsx
source("config.R")
library(readxl); library(dplyr); library(VennDiagram); library(grid)

CNV_FILE <- file.path(DATA_DIR, "CNV_DESeq2_Integrated_Analysis_Updated_latest.xlsx")
DEG_FILE <- file.path(DATA_DIR, "Cell_line_specific_merged_K1.xlsx")

cnv_s     <- read_excel(CNV_FILE, sheet = "CNV_DESeq2_CHO_S")
cnv_dxb11 <- read_excel(CNV_FILE, sheet = "CNV_DESeq2_CHO_DXB11")
cnv_dg44  <- read_excel(CNV_FILE, sheet = "CNV_DESeq2_CHO_DG44")
deg_s_up     <- read_excel(DEG_FILE, sheet = "Upregulated_K1 vs S")
deg_s_dn     <- read_excel(DEG_FILE, sheet = "Downregulated_K1 vs S")
deg_dxb_up   <- read_excel(DEG_FILE, sheet = "Upregulated_K1 vs DXB11")
deg_dxb_dn   <- read_excel(DEG_FILE, sheet = "Downregulated_K1 vs DXB11")
deg_dg44_up  <- read_excel(DEG_FILE, sheet = "Upregulated_K1 vs DG44")
deg_dg44_dn  <- read_excel(DEG_FILE, sheet = "Downregulated_K1 vs DG44")

make_venn <- function(deg_up, deg_dn, cnv_df, cnv_col, comparison, out_file) {
  transcript_genes <- unique(c(deg_up$Gene, deg_dn$Gene))
  cnv_genes <- cnv_df %>% filter(.data[[cnv_col]] != "Normal") %>% pull(gene) %>% unique()
  t_only  <- length(setdiff(transcript_genes, cnv_genes))
  c_only  <- length(setdiff(cnv_genes, transcript_genes))
  overlap <- length(intersect(transcript_genes, cnv_genes))
  cat(comparison, "-> Transcript only:", t_only, "| Overlap:", overlap, "| CNV only:", c_only, "\n")
  pdf(out_file, width = 6, height = 5)
  draw.pairwise.venn(area1 = t_only + overlap, area2 = c_only + overlap, cross.area = overlap,
    category = c("", ""), fill = c("#B3CDE3", "#FBB4AE"), alpha = 0.5,
    col = c("#2166AC", "#B2182B"), lwd = 1.5, cex = 1.6, fontfamily = "sans",
    scaled = TRUE, euler.d = TRUE, print.mode = "raw")
  grid.text(comparison, x = 0.15, y = 0.45, gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "sans"))
  grid.text("Transcript only", x = 0.42, y = 0.92, gp = gpar(fontsize = 11, fontface = "bold", col = "#B2182B", fontfamily = "sans"))
  grid.text("CNV only", x = 0.72, y = 0.08, gp = gpar(fontsize = 11, fontface = "bold", col = "#B2182B", fontfamily = "sans"))
  grid.text("Transcript and\nCNV", x = 0.55, y = 0.20, gp = gpar(fontsize = 9, fontface = "bold", fontfamily = "sans"))
  dev.off()
  message("Saved -> ", out_file)
}
make_venn(deg_dxb_up,  deg_dxb_dn,  cnv_dxb11, "CHO_DXB11_CNV", "DXB11 vs K1", file.path(OUT_DIR, "Venn_DXB11_vs_K1.pdf"))
make_venn(deg_dg44_up, deg_dg44_dn, cnv_dg44,  "CHO_DG44_CNV",  "DG44 vs K1",  file.path(OUT_DIR, "Venn_DG44_vs_K1.pdf"))
make_venn(deg_s_up,    deg_s_dn,    cnv_s,     "CHO_S_CNV",     "S vs K1",     file.path(OUT_DIR, "Venn_S_vs_K1.pdf"))
