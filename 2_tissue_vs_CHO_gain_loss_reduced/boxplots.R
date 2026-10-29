# Paths come from config.R in the repository root; run this script from the repository root.
source("config.R")
# Box plots: tissue vs CHO for the "gained" (4 genes) and "lost"
# (representative serpin / viral genes) sets.
#
# Inputs (data/): expression_selected_genes.csv, sample_metadata.csv, gene_sets.csv
# Outputs go to results/figures and results/stats
#
# Design (matches the subset analysis in CHO_G_L_R_subsets):
#   * Tissue box  = the 16 tissue runs.
#   * CHO box     = the 64 CHO runs (the 4 random subsets x 16 runs, pooled).
#   * Test        = Mann-Whitney U (wilcox.test), done PER SUBSET
#                   (16 CHO runs vs the same 16 tissue runs), Holm-adjusted
#                   across all gene x subset comparisons within a figure.
#   * Stars       = the LEAST significant of the 4 per-subset adjusted
#                   p-values for that gene (conservative).
#   * Values      = log2(normalized counts + 1)
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
})

# ---- settings you may want to change -----------------------------
col_tissue <- "#CC79A7"
col_cho    <- "#009E73"
dpi_png    <- 600
out_fig    <- file.path(OUT_DIR, "figures")
out_stats  <- file.path(OUT_DIR, "stats")

dir.create(out_fig,   showWarnings = FALSE, recursive = TRUE)
dir.create(out_stats, showWarnings = FALSE, recursive = TRUE)

expr <- read_csv(file.path(DATA_DIR, "expression_selected_genes.csv"), show_col_types = FALSE)
meta <- read_csv(file.path(DATA_DIR, "sample_metadata.csv"),           show_col_types = FALSE)
sets <- read_csv(file.path(DATA_DIR, "gene_sets.csv"),                 show_col_types = FALSE)

long <- expr %>%
  pivot_longer(-matrix_row, names_to = "run", values_to = "norm_counts") %>%
  inner_join(meta, by = "run") %>%
  inner_join(sets, by = "matrix_row") %>%
  mutate(log2expr = log2(norm_counts + 1),
         group    = factor(group, levels = c("Tissue", "CHO")))

star <- function(p) {
  ifelse(p < 1e-4, "****",
  ifelse(p < 1e-3, "***",
  ifelse(p < 1e-2, "**",
  ifelse(p < 0.05, "*", "ns"))))
}

all_stats <- list()

for (fig in unique(sets$figure)) {
  d <- long %>% filter(figure == fig)
  gene_levels <- sets %>% filter(figure == fig) %>% arrange(order) %>% pull(gene)
  d$gene <- factor(d$gene, levels = gene_levels)

  # ---- per-subset tests: 16 CHO vs the 16 tissues -------------------
  res <- d %>% filter(group == "CHO") %>% distinct(gene, subset) %>%
    arrange(gene, subset) %>% as.data.frame()
  res$p <- NA_real_
  res$CHO_median_log2 <- NA_real_
  res$tissue_median_log2 <- NA_real_
  for (i in seq_len(nrow(res))) {
    g <- res$gene[i]; s <- res$subset[i]
    x <- d$log2expr[d$gene == g & d$group == "CHO"    & d$subset == s]
    y <- d$log2expr[d$gene == g & d$group == "Tissue"]
    res$p[i] <- suppressWarnings(wilcox.test(x, y, exact = FALSE, correct = TRUE)$p.value)
    res$CHO_median_log2[i]    <- median(x)
    res$tissue_median_log2[i] <- median(y)
  }
  res$Holm_adj_p <- p.adjust(res$p, method = "holm")
  res$figure <- fig
  res$n_CHO <- 16; res$n_tissue <- 16
  all_stats[[fig]] <- res

  worst <- res %>% group_by(gene) %>%
    summarise(worst_adj_p = max(Holm_adj_p), .groups = "drop") %>%
    mutate(stars = star(worst_adj_p))

  # ---- bracket / star positions ------------------------------------
  top <- d %>% group_by(gene) %>% summarise(ymax = max(log2expr), .groups = "drop") %>%
    left_join(worst, by = "gene") %>%
    mutate(x = as.numeric(gene), y0 = ymax + 0.25, y1 = ymax + 0.5)
  ylim_top <- max(top$ymax) * 1.12 + 1.5

  p <- ggplot(d, aes(x = gene, y = log2expr, fill = group)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.65,
                 outlier.shape = NA, alpha = 0.6, linewidth = 0.4)+
    geom_point(aes(group = group),
               position = position_jitterdodge(jitter.width = 0.18, dodge.width = 0.8,
                                               seed = 7),
               size = 0.6, alpha = 0.45, colour = "black", stroke = 0) +
    geom_segment(data = top, inherit.aes = FALSE,
                 aes(x = x - 0.2, xend = x - 0.2, y = y0, yend = y1), linewidth = 0.3) +
    geom_segment(data = top, inherit.aes = FALSE,
                 aes(x = x + 0.2, xend = x + 0.2, y = y0, yend = y1), linewidth = 0.3) +
    geom_segment(data = top, inherit.aes = FALSE,
                 aes(x = x - 0.2, xend = x + 0.2, y = y1, yend = y1), linewidth = 0.3) +
    geom_text(data = top, inherit.aes = FALSE,
              aes(x = x, y = y1 + 0.05, label = stars),
              vjust = 0, size = 3.6, fontface = "bold") +
    scale_fill_manual(values = c(Tissue = col_tissue, CHO = col_cho),
                      labels = c(Tissue = "Tissues (n=16)",
                                 CHO = "CHO, 4 subsets pooled (n=64)")) +
    scale_y_continuous(limits = c(-0.5, ylim_top), expand = c(0, 0)) +
    labs(x = NULL, y = expression(bold(log[2]*"(normalized counts + 1)")),
         title = unique(d$title), fill = NULL) +
    theme_classic(base_size = 10) +
    theme(axis.title   = element_text(face = "bold"),
          axis.text    = element_text(face = "bold", colour = "black"),
          axis.text.x  = element_text(angle = 35, hjust = 1),
          plot.title   = element_text(face = "bold", hjust = 0.5, size = 10),
          legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
          legend.direction = "horizontal", legend.background = element_blank(),
          legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 7.5))

  w <- 1.25 * length(gene_levels) + 1.3; h <- 4.4
  png_dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggsave(file.path(out_fig, paste0("final_", fig, ".png")), p, width = w, height = h,
         dpi = dpi_png, device = png_dev, bg = "white")
  ggsave(file.path(out_fig, paste0("final_", fig, ".pdf")), p, width = w, height = h,
         device = "pdf")
}

stats <- bind_rows(all_stats) %>%
  mutate(stars = star(Holm_adj_p)) %>%
  select(figure, gene, subset, n_CHO, n_tissue, CHO_median_log2, tissue_median_log2,
         MWU_p = p, Holm_adj_p, stars)
write_csv(stats, file.path(out_stats, "per_subset_stats.csv"))

summary_stats <- stats %>% group_by(figure, gene) %>%
  summarise(worst_subset_Holm_adj_p = max(Holm_adj_p), .groups = "drop") %>%
  mutate(stars = star(worst_subset_Holm_adj_p))
write_csv(summary_stats, file.path(out_stats, "stars_used_in_figures.csv"))
print(as.data.frame(summary_stats))
