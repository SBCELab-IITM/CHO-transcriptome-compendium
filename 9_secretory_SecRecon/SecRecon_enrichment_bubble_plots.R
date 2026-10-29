# Paths come from config.R in the repository root; run this script from the repository root.
source("config.R")
library(tidyverse)
library(ggplot2)
library(stringr)
library(readr)
library(readxl)

# =====================================================
# MASTER THEME - Matches GO Enrichment Loss reference
# =====================================================
custom_palette <- c("#E06663", "#C36E79", "#A3758E", "#7A7AA4", "#327EBA")

paper_theme <- theme_minimal(base_family = "sans") +
  theme(
    plot.title       = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    axis.text.y      = element_text(size = 9, lineheight = 1.3, color = "black", hjust = 1, margin = margin(r = 4)),
    axis.text.x      = element_text(size = 9, color = "black"),
    axis.title.x     = element_text(size = 11, face = "bold", margin = margin(t = 10)),
    axis.title.y     = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linetype = "solid", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "gray80", linetype = "dashed", linewidth = 0.3),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.position  = "right",
    legend.box       = "vertical",
    legend.text      = element_text(size = 9),
    legend.title     = element_text(size = 10, face = "bold"),
    legend.spacing.y = unit(0.5, "cm"),
    plot.margin      = margin(l = 10, r = 15, t = 15, b = 10)
  )

# =====================================================
# STANDARDIZED SAVE - fixed width, dynamic height
# =====================================================
PLOT_WIDTH <- 12
PLOT_DPI   <- 300
ROW_HEIGHT <- 0.55
MIN_HEIGHT <- 8

save_plot <- function(p, out_pdf, n_rows) {
  h <- max(MIN_HEIGHT, n_rows * ROW_HEIGHT + 3.5)
  ggsave(out_pdf, p, width = PLOT_WIDTH, height = h, dpi = PLOT_DPI)
  message("Saved -> ", out_pdf, "  (", round(h, 1), " in tall)")
  print(p)
}

# =====================================================
# HELPER FUNCTIONS
# =====================================================
parse_gene_ratio <- function(x) {
  if (grepl("/", x)) {
    p <- strsplit(x, "/")[[1]]
    as.numeric(p[1]) / as.numeric(p[2])
  } else {
    as.numeric(x)
  }
}

match_compulsory <- function(description, compulsory_terms) {
  any(str_detect(tolower(description), tolower(compulsory_terms)))
}

force_top15_total <- function(df, compulsory_terms, total_n = 15) {
  compulsory_hits <- df %>%
    filter(map_lgl(Description, match_compulsory, compulsory_terms = compulsory_terms))
  remaining <- df %>%
    filter(!Description %in% compulsory_hits$Description) %>%
    arrange(p.adjust)
  cc_count  <- nrow(compulsory_hits %>% filter(GO_type == "CC"))
  needed_cc <- max(0, 2 - cc_count)
  cc_pathways <- remaining %>% filter(GO_type == "CC") %>% slice_head(n = needed_cc)
  remaining_other <- remaining %>% filter(!Description %in% cc_pathways$Description)
  bind_rows(
    compulsory_hits,
    cc_pathways,
    remaining_other %>% slice_head(n = max(0, total_n - nrow(compulsory_hits) - nrow(cc_pathways)))
  ) %>%
    distinct(Description, .keep_all = TRUE) %>%
    arrange(p.adjust)
}

# =====================================================
# INPUT FILES - UPDATE THESE PATHS
# =====================================================
LOSS_FILE      <- file.path(DATA_DIR, "Final_loss.csv")
RED_FILE       <- file.path(DATA_DIR, "Final_Enrich_Red.csv")
KEGG_LOSS_FILE <- file.path(DATA_DIR, "CHO_KEGG_enrichment_with_genes_loss.csv")
KEGG_RED_FILE  <- file.path(DATA_DIR, "CHO_KEGG_enrichment_with_genes_red.csv")
SYSTEM_FILE    <- file.path(DATA_DIR, "Loss_Reduced_System_Subsystem_Analysis.xlsx")
SUBSYS_FILE    <- file.path(DATA_DIR, "Red_Subsystem_Intersection.csv")

OUT_DIR <- OUT_DIR

# =====================================================
# COMPULSORY PATHWAYS
# =====================================================
loss_compulsory <- c(
  "Sialic acid binding", "response to lipid", "cellular response to lipid",
  "oxygen carrier activity", "sodium independent organic anion", "smooth muscle contraction"
)
reduced_compulsory <- c(
  "aldehyde dehydrogenase", "2'-5'-oligoadenylate synthetase",
  "cyclic nucleotide metabolic", "regulation of viral life cycle",
  "cyclic purine nucleotide metabolic"
)

# =============================================================================
# PLOT 1 & 2: GO ENRICHMENT (LOSS & REDUCED)
# =============================================================================
plot_go <- function(df, title, out_pdf) {
  
  df <- df %>%
    arrange(p.adjust) %>%
    mutate(
      desc_wrapped = str_wrap(Description, width = 38),
      label = if ("ID" %in% colnames(.)) paste0(ID, "\n", desc_wrapped) else desc_wrapped,
      label   = factor(label, levels = rev(unique(label))),
      GO_type = factor(GO_type, levels = c("BP", "CC", "MF")),
      bg_group = as.integer((row_number() - 1)) %% 2
    )
  
  p <- ggplot(df, aes(x = GeneRatio, y = label, fill = p.adjust, shape = GO_type, size = Count)) +
    geom_tile(aes(fill = as.factor(bg_group)), height = 0.92, alpha = 0.07, color = NA, show.legend = FALSE) +
    scale_fill_manual(values = c("0" = "white", "1" = "#f0f0f0"), guide = "none") +
    geom_point(aes(fill = p.adjust), color = "black", stroke = 0.5) +
    scale_shape_manual(values = c(BP = 21, CC = 22, MF = 24), name = "Ontology") +
    scale_fill_gradientn(
      colours = custom_palette,
      limits  = c(min(df$p.adjust), max(df$p.adjust)),
      trans   = "reverse",
      guide   = guide_colorbar(reverse = TRUE, title = "p.adj", title.position = "top", barwidth = 1.5, barheight = 8),
      name    = "p.adjust"
    ) +
    scale_size_continuous(
      name = "Count", breaks = c(5, 10, 15, 20, 25), range = c(2, 7),
      guide = guide_legend(title = "Count", title.position = "top")
    ) +
    labs(title = title, x = "GeneRatio", y = NULL) +
    paper_theme
  
  save_plot(p, out_pdf, nrow(df))
}

# --- Run GO Loss ---
df_loss <- read_csv(LOSS_FILE, show_col_types = FALSE) %>%
  mutate(GeneRatio = map_dbl(GeneRatio, parse_gene_ratio),
         Count = as.numeric(Count), p.adjust = as.numeric(p.adjust), GO_type = Term) %>%
  filter(!is.na(GeneRatio), !is.na(p.adjust), Count > 0)

df_loss_final <- force_top15_total(df_loss, loss_compulsory)
plot_go(df_loss_final, "GO Enrichment \u2014 Loss", paste0(OUT_DIR, "GO_dotplot_Loss.pdf"))

# --- Run GO Reduced ---
df_red <- read_csv(RED_FILE, show_col_types = FALSE) %>%
  mutate(GeneRatio = map_dbl(GeneRatio, parse_gene_ratio),
         Count = as.numeric(Count), p.adjust = as.numeric(p.adjust), GO_type = Term) %>%
  filter(!is.na(GeneRatio), !is.na(p.adjust), Count > 0)

df_red_final <- force_top15_total(df_red, reduced_compulsory)
plot_go(df_red_final, "GO Enrichment \u2014 Reduced", paste0(OUT_DIR, "GO_dotplot_Reduced.pdf"))

# =============================================================================
# PLOT 3 & 4: KEGG PATHWAY (LOSS & REDUCED)
# =============================================================================
plot_kegg <- function(df, title, out_pdf) {
  
  df <- df %>%
    arrange(p.adjust) %>%
    mutate(
      Description = factor(Description, levels = rev(unique(Description))),
      bg_group = as.integer((row_number() - 1)) %% 2
    )
  
  p <- ggplot(df, aes(x = GeneRatio, y = Description, fill = p.adjust, size = Count)) +
    geom_point(color = "black", stroke = 0.6, shape = 21) +
    scale_fill_gradientn(
      colours = rev(custom_palette),
      limits  = c(0, max(df$p.adjust) * 1.1),
      guide   = guide_colorbar(reverse = FALSE, title = "p.adj", title.position = "top", barwidth = 1.5, barheight = 8),
      name    = "p.adjust"
    ) +
    scale_size_continuous(
      name = "Count", breaks = seq(0, max(df$Count), by = 5), range = c(2, 8),
      guide = guide_legend(title = "Count", title.position = "top", override.aes = list(fill = "black", color = "black"))
    ) +
    labs(title = title, x = "GeneRatio", y = NULL) +
    paper_theme
  
  save_plot(p, out_pdf, nrow(df))
}

# --- Run KEGG Loss ---
df_kegg_loss <- read_csv(KEGG_LOSS_FILE, show_col_types = FALSE) %>%
  mutate(GeneRatio = map_dbl(GeneRatio, parse_gene_ratio),
         Count = as.numeric(Count), p.adjust = as.numeric(p.adjust)) %>%
  filter(!is.na(GeneRatio), !is.na(p.adjust), Count > 0) %>%
  arrange(p.adjust) %>% slice_head(n = 15)

plot_kegg(df_kegg_loss, "KEGG Pathway Enrichment (Top 15) \u2013 Loss", paste0(OUT_DIR, "KEGG_dotplot_Loss.pdf"))

# --- Run KEGG Reduced ---
df_kegg_red <- read_csv(KEGG_RED_FILE, show_col_types = FALSE) %>%
  mutate(GeneRatio = map_dbl(GeneRatio, parse_gene_ratio),
         Count = as.numeric(Count), p.adjust = as.numeric(p.adjust)) %>%
  filter(!is.na(GeneRatio), !is.na(p.adjust), Count > 0) %>%
  arrange(p.adjust) %>% slice_head(n = 15)

plot_kegg(df_kegg_red, "KEGG Pathway Enrichment (Top 15) \u2013 Reduced", paste0(OUT_DIR, "KEGG_dotplot_Reduced.pdf"))

# =============================================================================
# PLOT 5: SYSTEM ENRICHMENT (REDUCED)
# =============================================================================
df_system_raw <- read_excel(SYSTEM_FILE, sheet = "Reduced_System_Intersection")
total_genes_sys <- n_distinct(df_system_raw$`CHO GENE SYMBOL`)

df_system_top <- df_system_raw %>%
  group_by(System) %>%
  summarise(Count = n_distinct(`CHO GENE SYMBOL`), p.adjust = min(`p-val`, na.rm = TRUE), .groups = "drop") %>%
  mutate(GeneRatio = Count / total_genes_sys) %>%
  arrange(p.adjust) %>%
  slice_head(n = 15) %>%
  mutate(System = factor(System, levels = rev(System[order(GeneRatio)])))

p_sys <- ggplot(df_system_top, aes(x = GeneRatio, y = System, fill = p.adjust, size = Count)) +
  geom_point(shape = 21, color = "black", stroke = 0.6) +
  scale_fill_gradientn(
    colours = custom_palette,
    limits  = c(min(df_system_top$p.adjust), max(df_system_top$p.adjust)),
    guide   = guide_colorbar(title = "p.adj", title.position = "top", barwidth = 1.5, barheight = 8),
    name    = "p.adjust"
  ) +
  scale_size_continuous(
    name = "Gene count", range = c(2, 8),
    guide = guide_legend(title = "Gene count", title.position = "top", override.aes = list(fill = "black", color = "black"))
  ) +
  labs(title = "System Enrichment (Reduced Genes)", x = "Gene Ratio", y = NULL) +
  paper_theme

save_plot(p_sys, paste0(OUT_DIR, "System_Enrichment_Reduced.pdf"), nrow(df_system_top))

# =============================================================================
# PLOT 6: SUBSYSTEM ENRICHMENT (REDUCED)
# =============================================================================
red_subsystem_raw <- read.csv(SUBSYS_FILE, stringsAsFactors = FALSE)
gene_col <- colnames(red_subsystem_raw)[grep("CHO|GENE|SYMBOL", colnames(red_subsystem_raw), ignore.case = TRUE)][1]
p_col    <- colnames(red_subsystem_raw)[grep("p.val|p-val|pval", colnames(red_subsystem_raw), ignore.case = TRUE)][1]
total_genes_sub <- n_distinct(red_subsystem_raw[[gene_col]])

red_subsystem_top <- red_subsystem_raw %>%
  group_by(Subsystem) %>%
  summarise(Count = n_distinct(.data[[gene_col]]), p.adjust = min(.data[[p_col]], na.rm = TRUE), .groups = "drop") %>%
  mutate(GeneRatio = Count / total_genes_sub) %>%
  arrange(p.adjust) %>%
  slice_head(n = 15) %>%
  mutate(Subsystem = factor(Subsystem, levels = rev(Subsystem[order(GeneRatio)])))

p_sub <- ggplot(red_subsystem_top, aes(x = GeneRatio, y = Subsystem, fill = p.adjust, size = Count)) +
  geom_point(shape = 21, color = "black", stroke = 0.6, alpha = 0.8) +
  scale_fill_gradientn(
    colours = custom_palette,
    limits  = c(min(red_subsystem_top$p.adjust), max(red_subsystem_top$p.adjust)),
    guide   = guide_colorbar(title = "p.adj", title.position = "top", barwidth = 1.5, barheight = 8),
    name    = "p.adjust"
  ) +
  scale_size_continuous(
    name = "Gene count", range = c(3, 8),
    guide = guide_legend(title = "Gene count", title.position = "top", override.aes = list(fill = "black", color = "black"))
  ) +
  labs(title = "Subsystem Enrichment (Reduced Genes)", x = "Gene Ratio", y = NULL) +
  paper_theme

save_plot(p_sub, paste0(OUT_DIR, "Subsystem_Enrichment_Reduced.pdf"), nrow(red_subsystem_top))

cat("\n\u2705 All 6 bubble plots generated with standardized theme\n")
