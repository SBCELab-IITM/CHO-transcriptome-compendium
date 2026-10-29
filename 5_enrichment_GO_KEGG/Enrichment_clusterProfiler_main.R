# Paths come from config.R in the repository root; run this script from the repository root.
source("config.R")
# --- Install and load required packages ---
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pkgs <- c("clusterProfiler", "AnnotationHub", "dplyr", "enrichplot")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
  library(pkg, character.only = TRUE)
}

# --- Load input gene list ---
input_file <- file.path(DATA_DIR, "6.10.26 S PRM FC NO GENE ACT MULTI - Sheet1.csv")
gene_data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
gene_list <- unique(gene_data[[1]])
cat("Number of input genes:", length(gene_list), "\n")

# --- Retrieve Cricetulus griseus annotation ---
ah <- AnnotationHub::AnnotationHub()

cho_db <- ah[["AH120126"]]
cat("Annotation data retrieved:", metadata(cho_db)$Title, "\n")

# --- GO Enrichment Analysis (BP, MF, CC) ---
ontologies <- c("BP", "MF", "CC")
go_results <- list()

for (ont in ontologies) {
  cat("\nRunning GO enrichment for:", ont, "...\n")
  
  ego <- enrichGO(
    gene          = gene_list,
    OrgDb         = cho_db,
    keyType       = "SYMBOL",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  ego_simplified <- simplify(
    ego,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang"
  )
  
  go_results[[ont]] <- ego_simplified
  cat(ont, "GO terms enriched:", nrow(as.data.frame(ego_simplified)), "\n")
}

# --- Convert gene names to Entrez IDs before KEGG ---
cat("\nConverting gene symbols to Entrez IDs...\n")
gene_conversion <- clusterProfiler::bitr(
  gene_list,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = cho_db
)

cat("Gene name to Entrez ID mapping completed.\n")
cat("Mapped genes:", nrow(gene_conversion), "out of", length(gene_list), "\n")

# Extract unique Entrez IDs for KEGG
gene_entrez <- unique(gene_conversion$ENTREZID)

# --- KEGG Enrichment Analysis ---
cat("\nRunning KEGG enrichment...\n")
ekegg <- enrichKEGG(
  gene          = gene_entrez,
  organism      = "cge",  # Cricetulus griseus KEGG code
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

# --- Add gene SYMBOL mapping to KEGG output ---
if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
  ekegg_df <- as.data.frame(ekegg)
  
  ekegg_df$Gene_Symbols <- sapply(ekegg_df$geneID, function(x) {
    ids <- unlist(strsplit(x, "/"))
    matched <- gene_conversion$SYMBOL[match(ids, gene_conversion$ENTREZID)]
    paste(unique(na.omit(matched)), collapse = ", ")
  })
  
  cat("KEGG enrichment successful with", nrow(ekegg_df), "pathways.\n")
} else {
  ekegg_df <- data.frame()
  cat("No significant KEGG pathways found.\n")
}

# --- Save results ---
output_dir <- file.path(OUT_DIR, "6.10.26 S PRM FC NO GENE ACT")
dir.create(output_dir, showWarnings = FALSE)

# Save GO enrichment results
write.csv(as.data.frame(go_results$BP),
          file.path(output_dir, "CHO_GO_BP_enrichment.csv"),
          row.names = FALSE)

write.csv(as.data.frame(go_results$MF),
          file.path(output_dir, "CHO_GO_MF_enrichment.csv"),
          row.names = FALSE)

write.csv(as.data.frame(go_results$CC),
          file.path(output_dir, "CHO_GO_CC_enrichment.csv"),
          row.names = FALSE)

# Save KEGG enrichment with gene symbols
write.csv(ekegg_df,
          file.path(output_dir, "CHO_KEGG_enrichment_with_genes.csv"),
          row.names = FALSE)

cat("\n Enrichment results saved in:", output_dir, "\n")

# --- Visualization ---
cat("Generating dot plots...\n")

if (nrow(as.data.frame(go_results$BP)) > 0) {
  pdf(file.path(output_dir, "GO_BP_dotplot.pdf"))
  print(dotplot(go_results$BP, showCategory = 20, font.size = 10, title = "GO: Biological Process"))
  dev.off()
}

if (nrow(as.data.frame(go_results$MF)) > 0) {
  pdf(file.path(output_dir, "GO_MF_dotplot.pdf"))
  print(dotplot(go_results$MF, showCategory = 20, font.size = 10, title = "GO: Molecular Function"))
  dev.off()
}

if (nrow(as.data.frame(go_results$CC)) > 0) {
  pdf(file.path(output_dir, "GO_CC_dotplot.pdf"))
  print(dotplot(go_results$CC, showCategory = 20, font.size = 10, title = "GO: Cellular Component"))
  dev.off()
}

if (nrow(ekegg_df) > 0) {
  pdf(file.path(output_dir, "KEGG_dotplot.pdf"))
  print(dotplot(ekegg, showCategory = 20, font.size = 10, title = "KEGG Pathways"))
  dev.off()
}

cat("\n All analyses and plots completed successfully.\n")