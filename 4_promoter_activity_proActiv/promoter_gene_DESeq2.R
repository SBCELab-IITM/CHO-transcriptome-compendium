# promoter and gene level DESeq2, each cell line vs K1
# input: processed_promoter_data_final.csv (proActiv raw counts), Final_Metadata_06092025-2.txt
# output goes to results/promoter_deseq2/
# gene counts = sum of promoter counts per gene, K1 and K1(GS) treated as K1
# keeping features with >=10 counts in >=20 samples, design ~ Cell_line
# significant = padj < 0.05 and |log2FC| > 1
# (Study is not in the model because it overlaps with cell line)
source("config.R")
library(DESeq2)

MIN_COUNT <- 10; MIN_SAMPLES <- 20
outdir <- file.path(OUT_DIR, "promoter_deseq2"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# count matrices
pa_file <- file.path(DATA_DIR, "processed_promoter_data_final.csv")
hdr <- names(read.csv(pa_file, nrows = 1, check.names = FALSE))
raw_cols  <- grep("\\.raw_counts$", hdr, value = TRUE)
anno_cols <- c("promoterId", "geneId", "seqnames", "start", "strand", "internalPromoter", "promoterPosition")
pa <- read.csv(pa_file, check.names = FALSE)[, c(anno_cols, raw_cols)]
anno   <- pa[, anno_cols]
counts <- as.matrix(pa[, raw_cols]); colnames(counts) <- sub("\\..*$", "", raw_cols)
rm(pa)
cat(nrow(counts), "promoters,", ncol(counts), "samples\n")

meta <- read.delim(file.path(DATA_DIR, "Final_Metadata_06092025-2.txt"), stringsAsFactors = FALSE)
rownames(meta) <- meta$Run
samples <- colnames(counts)[colnames(counts) %in% rownames(meta)]
stopifnot(length(samples) == ncol(counts))              
counts <- counts[, samples]
coldata <- data.frame(Cell_line = sub("^K1 \\(GS\\)$", "K1", meta[samples, "Cell"]),
                      Study = meta[samples, "Study"], Media = meta[samples, "Media"], row.names = samples)
coldata$Cell_line <- factor(coldata$Cell_line, levels = c("K1", "S", "DG44", "DXB11"))
stopifnot(!anyNA(coldata$Cell_line)); print(table(coldata$Cell_line))

rownames(counts) <- anno$promoterId
gene <- rowsum(counts, group = anno$geneId, reorder = TRUE)                
keep_filter <- function(m, label) {
  k <- rowSums(m >= MIN_COUNT) >= MIN_SAMPLES
  cat(label, ":", nrow(m), "->", sum(k), "after filtering\n"); m[k, , drop = FALSE]
}
gene     <- keep_filter(gene, "gene")
promoter <- keep_filter(counts, "promoter")
write.csv(gene,     file.path(outdir, "counts_gene.csv"))
write.csv(promoter, file.path(outdir, "counts_promoter.csv"))
write.csv(coldata,  file.path(outdir, "sample_sheet.csv"), row.names = TRUE)
write.csv(anno,     file.path(outdir, "promoter_annotation.csv"), row.names = FALSE)

# deseq2
run_level <- function(level, cts) {
  mode(cts) <- "integer"
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Cell_line)
  dds <- DESeq(dds)
  rows <- list()
  for (a in c("S", "DG44", "DXB11")) {
    res <- results(dds, contrast = c("Cell_line", a, "K1"), alpha = 0.05)
    df  <- as.data.frame(res); df <- cbind(feature = rownames(df), df)
    if (level == "promoter") {
      m <- anno[match(df$feature, anno$promoterId), c("geneId", "promoterPosition", "internalPromoter")]
      df <- cbind(df[, "feature", drop = FALSE], m, df[, -1])
    }
    write.csv(df[order(df$padj), ], file.path(outdir, sprintf("%s_%s_vs_K1.csv", level, a)), row.names = FALSE)
    sig <- !is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) > 1
    rows[[a]] <- data.frame(contrast = paste(a, "vs K1"), level = level, tested = nrow(df), DE = sum(sig),
                            up_in_test = sum(sig & df$log2FoldChange > 0), up_in_K1 = sum(sig & df$log2FoldChange < 0))
  }
  do.call(rbind, rows)
}
summary_all <- rbind(run_level("gene", gene), run_level("promoter", promoter))
write.csv(summary_all, file.path(outdir, "DE_summary_vs_K1.csv"), row.names = FALSE)
print(summary_all)

# genes with more than one promoter
# only non-internal promoters; contrasting = gene has both an up and a down promoter
multi_genes <- names(which(table(anno$geneId) > 1))
alt <- do.call(rbind, lapply(c("S", "DG44", "DXB11"), function(a) {
  p <- read.csv(file.path(outdir, sprintf("promoter_%s_vs_K1.csv", a)))
  p <- p[p$geneId %in% multi_genes & !as.logical(p$internalPromoter) & !is.na(p$padj) & p$padj < 0.05 & abs(p$log2FoldChange) > 1, ]
  up <- unique(p$geneId[p$log2FoldChange > 0]); dn <- unique(p$geneId[p$log2FoldChange < 0])
  data.frame(contrast = paste(a, "vs K1"), sig_promoters = nrow(p), genes_up_only = length(setdiff(up, dn)),
             genes_down_only = length(setdiff(dn, up)), genes_contrasting = length(intersect(up, dn)))
}))
write.csv(alt, file.path(outdir, "alternative_promoter_counts.csv"), row.names = FALSE)
print(alt)
