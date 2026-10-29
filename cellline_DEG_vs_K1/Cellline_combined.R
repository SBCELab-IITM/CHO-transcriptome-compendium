# Paths come from config.R in the repository root; run this script from the repository root.
source("config.R")
# DESeq2: K1 (merged with K1-GS) vs S, DXB11, DG44
#
# K1 and K1(GS) are merged into a single "K1" reference group
# based on their minimal transcriptomic divergence.
#
# Samples removed: SRR29672926, SRR29672927
#   (CHO-DP12 in serum DMEM, misannotated as K1 NP, Batch 24)
#
# Output:
#   - Summary sheet
#   - Sheet 1-3  : All genes (K1 vs S, K1 vs DXB11, K1 vs DG44)
#   - Sheet 4-6  : Upregulated
#   - Sheet 7-9  : Downregulated
#
# Cutoffs: padj < 0.01, |log2FC| > 2

library(DESeq2)
library(dplyr)
library(openxlsx)

# 1. LOAD DATA

counts_raw <- read.csv(file.path(DATA_DIR, "Count_raw_all.csv"),
                       row.names = 1, header = TRUE, check.names = FALSE)

metadata <- read.delim(file.path(DATA_DIR, "Final_Metadata_06092025-2.txt"),
                       header = TRUE, stringsAsFactors = FALSE)

cat("Counts:", nrow(counts_raw), "genes x", ncol(counts_raw), "samples\n")
cat("Metadata:", nrow(metadata), "rows\n")

# 2. REMOVE BAD SAMPLES

bad_samples <- c("SRR29672926", "SRR29672927")
counts_raw  <- counts_raw[, !colnames(counts_raw) %in% bad_samples]
metadata    <- metadata[!metadata$Run %in% bad_samples, ]
cat("After removal:", ncol(counts_raw), "samples\n")

# 3. SUBSET TO CHO CELL LINES & MERGE K1 + K1(GS)

target_cells <- c("K1", "K1 (GS)", "S", "DXB11", "DG44")
metadata     <- metadata[metadata$Cell %in% target_cells, ]
counts_raw   <- counts_raw[, colnames(counts_raw) %in% metadata$Run]

# *** MERGE: Relabel K1(GS) as K1 ***
metadata$Cell[metadata$Cell == "K1 (GS)"] <- "K1"

# Align order
metadata           <- metadata[match(colnames(counts_raw), metadata$Run), ]
rownames(metadata) <- metadata$Run

cat("\nSample distribution (after merging K1 + K1-GS):\n")
print(table(metadata$Cell))

# 4. LOW-COUNT FILTER

keep        <- rowSums(counts_raw) >= 10
counts_filt <- counts_raw[keep, ]
cat("\nGenes after low-count filter:", nrow(counts_filt), "\n")

# CUTOFFS

padj_cutoff <- 0.01
lfc_cutoff  <- 2.0

# 5. DESeq2 COMBINED MODEL (4 cell lines, K1 as reference)

cat("\n=== COMBINED MODEL RUN ===\n")

metadata$Cell <- factor(metadata$Cell, levels = c("K1", "S", "DXB11", "DG44"))

dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = metadata,
  design    = ~ Cell
)
dds <- DESeq(dds)

# Define comparisons (K1 is reference - no K1 vs K1-GS anymore)
comparisons <- list(
  list(label = "K1 vs S",      test = "S",      ref = "K1"),
  list(label = "K1 vs DXB11",  test = "DXB11",  ref = "K1"),
  list(label = "K1 vs DG44",   test = "DG44",   ref = "K1")
)

all_results <- list()

for (comp in comparisons) {
  cat("\n", comp$label, "...\n")
  
  res    <- results(dds, contrast = c("Cell", comp$test, comp$ref))
  res_df <- as.data.frame(res)
  res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
  
  up   <- res_df %>%
    filter(!is.na(padj), log2FoldChange >  lfc_cutoff, padj < padj_cutoff)
  down <- res_df %>%
    filter(!is.na(padj), log2FoldChange < -lfc_cutoff, padj < padj_cutoff)
  
  all_results[[comp$label]] <- list(
    all = res_df, up = up, down = down,
    test = comp$test, ref = comp$ref
  )
  
  cat("  Upregulated  :", nrow(up), "\n")
  cat("  Downregulated:", nrow(down), "\n")
}

# 6. BUILD EXCEL

cat("\n=== Writing Excel ===\n")

wb <- createWorkbook()

add_data_sheet <- function(wb, sheet_name, df) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, cbind(Gene = rownames(df), df), rowNames = FALSE)
}

# summary
addWorksheet(wb, "Summary")

writeData(wb, "Summary",
          data.frame(V1 = "K1 and K1(GS) merged as K1 reference | Comparisons: K1 vs S, DXB11, DG44"),
          startRow = 1, startCol = 1, colNames = FALSE)

sum_df <- data.frame(
  "K1 vs"  = c("Upregulated", "Downregulated"),
  "S"      = c(nrow(all_results[["K1 vs S"]]$up),
               nrow(all_results[["K1 vs S"]]$down)),
  "DXB11"  = c(nrow(all_results[["K1 vs DXB11"]]$up),
               nrow(all_results[["K1 vs DXB11"]]$down)),
  "DG44"   = c(nrow(all_results[["K1 vs DG44"]]$up),
               nrow(all_results[["K1 vs DG44"]]$down)),
  check.names = FALSE
)
writeData(wb, "Summary", sum_df, startRow = 2, startCol = 1)

contents <- data.frame(
  Sheet   = paste("Sheet", 1:9),
  Content = c(
    "K1 vs S_All",        "K1 vs DXB11_All",        "K1 vs DG44_All",
    "Upregulated_K1 vs S","Upregulated_K1 vs DXB11","Upregulated_K1 vs DG44",
    "Downregulated_K1 vs S","Downregulated_K1 vs DXB11","Downregulated_K1 vs DG44"
  )
)
writeData(wb, "Summary", contents, startRow = 6, startCol = 1)

writeData(wb, "Summary",
          data.frame(Cutoffs = c("Pvalue", "LogFC"),
                     Value   = c(padj_cutoff, paste("< or >", lfc_cutoff))),
          startRow = 6, startCol = 4)

writeData(wb, "Summary",
          data.frame(Note = c(
            "K1(GS) samples relabeled as K1 and merged into reference group",
            "Removed: SRR29672926, SRR29672927 (CHO-DP12 misannotated as K1 NP, Batch 24)"
          )),
          startRow = 18, startCol = 1, colNames = FALSE)

# all
add_data_sheet(wb, "K1 vs S_All",     all_results[["K1 vs S"]]$all)
add_data_sheet(wb, "K1 vs DXB11_All", all_results[["K1 vs DXB11"]]$all)
add_data_sheet(wb, "K1 vs DG44_All",  all_results[["K1 vs DG44"]]$all)

# upregulated
add_data_sheet(wb, "Upregulated_K1 vs S",     all_results[["K1 vs S"]]$up)
add_data_sheet(wb, "Upregulated_K1 vs DXB11", all_results[["K1 vs DXB11"]]$up)
add_data_sheet(wb, "Upregulated_K1 vs DG44",  all_results[["K1 vs DG44"]]$up)

# downregulated
add_data_sheet(wb, "Downregulated_K1 vs S",     all_results[["K1 vs S"]]$down)
add_data_sheet(wb, "Downregulated_K1 vs DXB11", all_results[["K1 vs DXB11"]]$down)
add_data_sheet(wb, "Downregulated_K1 vs DG44",  all_results[["K1 vs DG44"]]$down)
# save
saveWorkbook(wb, file.path(OUT_DIR, "Cell_line_specific_merged_K1.xlsx"), overwrite = TRUE)
cat("Saved: Cell_line_specific_merged_K1.xlsx\n")

# 7. QC PLOTS

vsd <- vst(dds, blind = FALSE)

pdf(file.path(OUT_DIR, "PCA_cell_line_merged_K1.pdf"), width = 8, height = 6)
print(
  plotPCA(vsd, intgroup = "Cell") +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("PCA - K1 (merged with K1-GS) vs S, DXB11, DG44")
)
dev.off()

pdf(file.path(OUT_DIR, "MA_plots_cell_line_merged_K1.pdf"), width = 12, height = 5)
par(mfrow = c(1, 3))
for (comp in comparisons) {
  plotMA(results(dds, contrast = c("Cell", comp$test, comp$ref)),
         main = comp$label, ylim = c(-8, 8))
}
dev.off()

cat("QC plots saved\n")

# FINAL SUMMARY

cat("\n========== FINAL SUMMARY ==========\n")
cat("Design : K1 + K1(GS) merged as reference\n")
cat("Removed: SRR29672926, SRR29672927\n")
cat("Samples:", ncol(counts_filt), "\n")
cat("Genes  :", nrow(counts_filt), "\n")
cat("Cutoffs: padj <", padj_cutoff, "| |log2FC| >", lfc_cutoff, "\n\n")

for (nm in names(all_results)) {
  cat(sprintf("  %-20s  UP: %4d   DOWN: %4d\n", nm,
              nrow(all_results[[nm]]$up),
              nrow(all_results[[nm]]$down)))
}