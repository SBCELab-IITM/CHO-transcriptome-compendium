# Paths come from config.R in the repository root; run this script from the repository root.
source("config.R")
library(tidyverse)
library(ggplot2)
library(stringr)
library(readr)

# INPUT FILES
LOSS_FILE = file.path(DATA_DIR, "Final_loss.csv")
RED_FILE  = file.path(DATA_DIR, "Final_Enrich_Red.csv")

OUT_PDF_LOSS = file.path(OUT_DIR, "GO_dotplot_Loss.pdf")
OUT_PDF_RED  = file.path(OUT_DIR, "GO_dotplot_Reduced.pdf")

# COMPULSORY PATHWAYS (PARTIAL MATCH ALLOWED)
loss_compulsory <- c(
  "Sialic acid binding",
  "response to lipid",
  "cellular response to lipid",
  "oxygen carrier activity",
  "sodium independent organic anion",
  "smooth muscle contraction"
)

reduced_compulsory <- c(
  "aldehyde dehydrogenase",
  "2'-5'-oligoadenylate synthetase",
  "cyclic nucleotide metabolic",
  "regulation of viral life cycle",
  "cyclic purine nucleotide metabolic"
)

# CUSTOM COLOR PALETTE
custom_palette <- c(
  "#E06663",
  "#C36E79",
  "#A3758E",
  "#7A7AA4",
  "#327EBA"
)

# HELPER FUNCTIONS

parse_gene_ratio <- function(x) {
  if (grepl("/", x)) {
    p <- strsplit(x, "/")[[1]]
    as.numeric(p[1]) / as.numeric(p[2])
  } else {
    as.numeric(x)
  }
}

infer_go_type <- function(desc) {
  desc <- tolower(desc)
  case_when(
    str_detect(desc, "activity") ~ "MF",
    str_detect(desc, "process|regulation|response") ~ "BP",
    TRUE ~ "CC"
  )
}

match_compulsory <- function(description, compulsory_terms) {
  any(str_detect(
    tolower(description),
    tolower(compulsory_terms)
  ))
}

force_top15_total <- function(df, compulsory_terms, total_n = 15) {
  
  compulsory_hits <- df %>%
    filter(map_lgl(
      Description,
      match_compulsory,
      compulsory_terms = compulsory_terms
    ))
  
  remaining <- df %>%
    filter(!Description %in% compulsory_hits$Description) %>%
    arrange(p.adjust)
  
  # Ensure at least 1-2 CC (Cellular Component) pathways
  cc_count <- nrow(compulsory_hits %>% filter(GO_type == "CC"))
  needed_cc <- max(0, 2 - cc_count)
  
  cc_pathways <- remaining %>%
    filter(GO_type == "CC") %>%
    slice_head(n = needed_cc)
  
  remaining_other <- remaining %>%
    filter(!Description %in% cc_pathways$Description)
  
  # Combine all and fill remaining slots
  result <- bind_rows(
    compulsory_hits,
    cc_pathways,
    remaining_other %>% slice_head(n = max(0, total_n - nrow(compulsory_hits) - nrow(cc_pathways)))
  ) %>%
    distinct(Description, .keep_all = TRUE) %>%
    arrange(p.adjust)
  
  return(result)
}

# IMPROVED PLOTTING FUNCTION WITH ROW SEPARATION
plot_go <- function(df, title, out_pdf) {
  
  df <- df %>%
    arrange(p.adjust) %>%
    mutate(
      Description = factor(Description, levels = rev(unique(Description))),
      GO_type = factor(GO_type, levels = c("BP", "CC", "MF")),
      row_number = row_number()
    )
  
  # Create alternating background for row separation
  df <- df %>%
    mutate(
      bg_group = as.integer((row_number() - 1) / 1) %% 2
    )
  
  p <- ggplot(
    df,
    aes(
      x = GeneRatio,
      y = Description,
      fill = p.adjust,
      shape = GO_type,
      size = Count
    )
  ) +
    # Alternating row background for clarity
    geom_tile(
      aes(fill = as.factor(bg_group)),
      height = 0.95,
      alpha = 0.08,
      color = NA,
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = c("0" = "white", "1" = "gray"),
      guide = "none"
    ) +
    
    # Actual data points on top
    geom_point(aes(fill = p.adjust), color = "black", stroke = 0.5) +
    
    scale_shape_manual(
      values = c(BP = 21, CC = 22, MF = 24),
      name = "Ontology"
    ) +
    
    # Custom color palette for p.adjust
    scale_fill_gradientn(
      colours = custom_palette,
      limits = c(min(df$p.adjust), max(df$p.adjust)),
      trans = "reverse",
      guide = guide_colorbar(
        reverse = TRUE,
        title = "p.adj",
        title.position = "top",
        barwidth = 1.5,
        barheight = 8
      ),
      name = "p.adjust"
    ) +
    
    scale_size_continuous(
      name = "Count",
      breaks = c(5, 10, 15, 20, 25),
      range = c(2, 7),
      guide = guide_legend(
        title = "Count",
        title.position = "top"
      )
    ) +
    
    labs(
      title = title,
      x = "GeneRatio",
      y = NULL
    ) +
    
    theme_minimal() +
    theme(
      # Title
      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 15)
      ),
      
      # Axes
      axis.text.y = element_text(
        size = 9,
        lineheight = 1.3,
        color = "black"
      ),
      axis.text.x = element_text(size = 9, color = "black"),
      axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_blank(),
      
      # Grid
      panel.grid.major.y = element_line(
        color = "gray90",
        linetype = "solid",
        size = 0.3
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(
        linetype = "dashed",
        color = "gray80",
        size = 0.3
      ),
      
      # Panel border
      panel.border = element_rect(
        color = "black",
        fill = NA,
        size = 0.6
      ),
      
      # Legend
      legend.position = "right",
      legend.box = "vertical",
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10, face = "bold"),
      legend.spacing.y = unit(0.5, "cm"),
      
      # Plot margins
      plot.margin = margin(l = 10, r = 15, t = 15, b = 10)
    )
  
  # Dynamic height based on number of pathways
  plot_height <- max(10, nrow(df) * 0.5)
  
  ggsave(
    out_pdf,
    p,
    width = 12,
    height = plot_height,
    dpi = 300
  )
  
  print(p)
}

# LOSS PIPELINE
df_loss <- read_csv(LOSS_FILE, show_col_types = FALSE) %>%
  mutate(
    GeneRatio = map_dbl(GeneRatio, parse_gene_ratio),
    Count = as.numeric(Count),
    p.adjust = as.numeric(p.adjust),
    GO_type = Term
  ) %>%
  filter(!is.na(GeneRatio), !is.na(p.adjust), Count > 0)

df_loss_final <- force_top15_total(df_loss, loss_compulsory)

cat("Loss pathways selected:\n")
print(df_loss_final %>% select(Description, GO_type, p.adjust, Count))

plot_go(
  df_loss_final,
  "GO Enrichment",
  OUT_PDF_LOSS
)

# REDUCED PIPELINE
df_red <- read_csv(RED_FILE, show_col_types = FALSE) %>%
  mutate(
    GeneRatio = map_dbl(GeneRatio, parse_gene_ratio),
    Count = as.numeric(Count),
    p.adjust = as.numeric(p.adjust),
    GO_type = Term
  ) %>%
  filter(!is.na(GeneRatio), !is.na(p.adjust), Count > 0)

df_red_final <- force_top15_total(df_red, reduced_compulsory)

cat("\nReduced pathways selected:\n")
print(df_red_final %>% select(Description, GO_type, p.adjust, Count))

plot_go(
  df_red_final,
  "GO Enrichment",
  OUT_PDF_RED
)

cat("\nLoss and Reduced GO plots generated successfully\n")

# KEGG PATHWAY

library(tidyverse)
library(ggplot2)
library(stringr)
library(readr)

# INPUT FILES
KEGG_LOSS_FILE <- file.path(DATA_DIR, "CHO_KEGG_enrichment_with_genes_loss.csv")
KEGG_RED_FILE  <- file.path(DATA_DIR, "CHO_KEGG_enrichment_with_genes_red.csv")

OUT_PDF_KEGG_LOSS <- file.path(OUT_DIR, "KEGG_dotplot_Loss_1.pdf")
OUT_PDF_KEGG_RED  <- file.path(OUT_DIR, "KEGG_dotplot_Reduced_1.pdf")

# CUSTOM COLOR PALETTE
custom_palette <- c(
  "#E06663",
  "#C36E79",
  "#A3758E",
  "#7A7AA4",
  "#327EBA"
)

# HELPER FUNCTIONS

parse_gene_ratio <- function(x) {
  if (grepl("/", x)) {
    p <- strsplit(x, "/")[[1]]
    as.numeric(p[1]) / as.numeric(p[2])
  } else {
    as.numeric(x)
  }
}

# SELECT TOP N PATHWAYS
select_top_pathways <- function(df, n = 15) {
  df %>%
    arrange(p.adjust) %>%
    slice_head(n = n) %>%
    arrange(p.adjust)
}

# KEGG PLOTTING FUNCTION
plot_kegg <- function(df, title, out_pdf) {
  
  df <- df %>%
    mutate(
      Description = factor(Description, levels = rev(unique(Description))),
      row_number = row_number()
    )
  
  # Create alternating background for row separation
  df <- df %>%
    mutate(
      bg_group = as.integer((row_number() - 1) / 1) %% 2
    )
  
  p <- ggplot(
    df,
    aes(
      x = GeneRatio,
      y = Description,
      fill = p.adjust,
      size = Count
    )
  ) +
    # Actual data points - circles with p.adjust color FILLED
    geom_point(color = "black", stroke = 0.6, shape = 21) +
    
    # Color scale for p.adjust - RED = Highly Significant
    scale_fill_gradientn(
      colours = rev(custom_palette),
      limits = c(0, max(df$p.adjust) * 1.1),
      guide = guide_colorbar(
        reverse = FALSE,
        title = "p.adj",
        title.position = "top",
        barwidth = 1.5,
        barheight = 8
      ),
      name = "p.adjust"
    ) +
    
    scale_size_continuous(
      name = "Count",
      breaks = seq(0, max(df$Count), by = 5),
      range = c(2, 8),
      guide = guide_legend(
        title = "Count",
        title.position = "top",
        override.aes = list(fill = "black", color = "black")
      )
    ) +
    
    labs(
      title = title,
      x = "GeneRatio",
      y = NULL
    ) +
    
    theme_minimal() +
    theme(
      # Title
      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 15)
      ),
      
      # Axes
      axis.text.y = element_text(
        size = 9,
        lineheight = 1.3,
        color = "black"
      ),
      axis.text.x = element_text(size = 9, color = "black"),
      axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_blank(),
      
      # Grid
      panel.grid.major.y = element_line(
        color = "gray90",
        linetype = "solid",
        size = 0.3
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(
        linetype = "dashed",
        color = "gray80",
        size = 0.3
      ),
      
      # Panel border
      panel.border = element_rect(
        color = "black",
        fill = NA,
        size = 0.6
      ),
      
      # Legend
      legend.position = "right",
      legend.box = "vertical",
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10, face = "bold"),
      legend.spacing.y = unit(0.5, "cm"),
      
      # Plot margins
      plot.margin = margin(l = 10, r = 15, t = 15, b = 10)
    )
  
  # Dynamic height based on number of pathways
  plot_height <- max(8, nrow(df) * 0.5)
  
  ggsave(
    out_pdf,
    p,
    width = 12,
    height = plot_height,
    dpi = 300
  )
  
  print(p)
}

# KEGG LOSS PIPELINE
cat("Processing KEGG Loss file...\n")
df_kegg_loss <- read_csv(KEGG_LOSS_FILE, show_col_types = FALSE) %>%
  mutate(
    GeneRatio = map_dbl(GeneRatio, parse_gene_ratio),
    Count = as.numeric(Count),
    p.adjust = as.numeric(p.adjust)
  ) %>%
  filter(!is.na(GeneRatio), !is.na(p.adjust), Count > 0)

df_kegg_loss_final <- select_top_pathways(df_kegg_loss, n = 15)

cat("Loss pathways selected:\n")
print(df_kegg_loss_final %>% select(Description, p.adjust, Count, GeneRatio))

plot_kegg(
  df_kegg_loss_final,
  "KEGG Pathway Enrichment (Top 15) - Loss",
  OUT_PDF_KEGG_LOSS
)

# KEGG REDUCED PIPELINE
cat("\nProcessing KEGG Reduced file...\n")
df_kegg_red <- read_csv(KEGG_RED_FILE, show_col_types = FALSE) %>%
  mutate(
    GeneRatio = map_dbl(GeneRatio, parse_gene_ratio),
    Count = as.numeric(Count),
    p.adjust = as.numeric(p.adjust)
  ) %>%
  filter(!is.na(GeneRatio), !is.na(p.adjust), Count > 0)

df_kegg_red_final <- select_top_pathways(df_kegg_red, n = 15)

cat("Reduced pathways selected:\n")
print(df_kegg_red_final %>% select(Description, p.adjust, Count, GeneRatio))

plot_kegg(
  df_kegg_red_final,
  "KEGG Pathway Enrichment ",
  OUT_PDF_KEGG_RED
)

cat("\nKEGG plots generated successfully!\n")

#version 3
# LIBRARIES
library(tidyverse)
library(ggplot2)
library(readxl)

# INPUT / OUTPUT
SYSTEM_FILE <- file.path(DATA_DIR, "Loss_Reduced_System_Subsystem_Analysis.xlsx")
OUT_PDF_SYSTEM <- file.path(OUT_DIR, "System_Enrichment_Reduced.pdf")

# COLOR PALETTE
custom_palette <- c(
  "#E06663",
  "#C36E79",
  "#A3758E",
  "#7A7AA4",
  "#327EBA"
)

# SYSTEM ENRICHMENT DATA PREP
df_system_raw <- read_excel(
  SYSTEM_FILE,
  sheet = "Reduced_System_Intersection"
)

# Total reduced genes
total_genes <- n_distinct(df_system_raw$`CHO GENE SYMBOL`)

# Collapse gene-level -> system-level
df_system_enrich <- df_system_raw %>%
  group_by(System) %>%
  summarise(
    Count = n_distinct(`CHO GENE SYMBOL`),
    p.adjust = min(`p-val`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    GeneRatio = Count / total_genes,
    System = factor(System, levels = rev(System[order(GeneRatio)]))
  ) %>%
  arrange(p.adjust)

# Select top systems (optional)
df_system_top <- df_system_enrich %>%
  slice_head(n = 15)

# SYSTEM ENRICHMENT DOTPLOT
p <- ggplot(
  df_system_top,
  aes(
    x = GeneRatio,
    y = System,
    fill = p.adjust,
    size = Count
  )
) +
  geom_point(shape = 21, color = "black", stroke = 0.6) +
  
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    name = "p.adj"
  ) +
  
  scale_size_continuous(
    name = "Gene count",
    range = c(2, 8),
    guide = guide_legend(
      override.aes = list(
        fill = "black",
        color = "black"
      )
    )
  ) +
  
  labs(
    title = "System Enrichment (Reduced Genes)",
    x = "Gene Ratio",
    y = NULL
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold"),
    panel.grid.major.x = element_line(linetype = "dashed", color = "gray80"),
    panel.border = element_rect(color = "black", fill = NA)
  )

ggsave(
  OUT_PDF_SYSTEM,
  p,
  width = 12,
  height = max(8, nrow(df_system_top) * 0.5),
  dpi = 600
)

print(p)

cat("\nSystem enrichment plot generated successfully\n")

# version latest
library(dplyr)
library(ggplot2)

# READ DATA FROM CSV
red_subsystem_raw <- read.csv(file.path(DATA_DIR, "Red_Subsystem_Intersection.csv"), stringsAsFactors = FALSE)

# Check column names
cat("\nColumn names in your CSV:\n")
print(colnames(red_subsystem_raw))

# GET TOTAL UNIQUE REDUCED GENES
# Adjust column name based on what you have (CHO GENE SYMBOL or CHO.GENE.SYMBOL)
gene_col <- colnames(red_subsystem_raw)[grep("CHO|GENE|SYMBOL", colnames(red_subsystem_raw), ignore.case = TRUE)][1]
p_col <- colnames(red_subsystem_raw)[grep("p.val|p-val|pval", colnames(red_subsystem_raw), ignore.case = TRUE)][1]

cat("\nUsing gene column:", gene_col, "\n")
cat("Using p-value column:", p_col, "\n\n")

total_genes <- n_distinct(red_subsystem_raw[[gene_col]])
cat("Total unique genes:", total_genes, "\n")

# COLLAPSE GENE-LEVEL -> SUBSYSTEM-LEVEL
red_subsystem_enrich <- red_subsystem_raw %>%
  group_by(Subsystem) %>%
  summarise(
    Count = n_distinct(.data[[gene_col]]),
    p.adjust = min(.data[[p_col]], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(GeneRatio = Count / total_genes) %>%
  arrange(p.adjust)

cat("\nUnique subsystems in your data:\n")
print(unique(red_subsystem_raw$Subsystem))

# SELECT TOP 15 BY P-VALUE
red_subsystem_top <- red_subsystem_enrich %>%
  slice_head(n = 15) %>%
  mutate(
    Subsystem = factor(
      Subsystem,
      levels = rev(Subsystem[order(GeneRatio)])
    )
  )

# SUBSYSTEM ENRICHMENT DOTPLOT
p <- ggplot(
  red_subsystem_top,
  aes(
    x = GeneRatio,
    y = Subsystem,
    fill = p.adjust,
    size = Count
  )
) +
  geom_point(shape = 21, color = "black", stroke = 0.6, alpha = 0.8) +
  
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    name = "p.adj"
  ) +
  
  scale_size_continuous(
    name = "Gene count",
    range = c(3, 8)
  ) +
  
  guides(fill = guide_colorbar(order = 1),
         size = guide_legend(order = 2, override.aes = list(fill = "black"))) +
  
  labs(
    title = "Subsystem Enrichment (Reduced Genes)",
    x = "Gene Ratio",
    y = NULL
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold"),
    panel.grid.major.x = element_line(linetype = "dashed", color = "gray80"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "right"
  )

print(p)

# SAVE FIGURE (600 DPI)
ggsave(
  "Reduced_Subsystem_Enrichment.pdf",
  p,
  width = 11,
  height = max(8, nrow(red_subsystem_top) * 0.5),
  dpi = 600
)

cat("\nPlot saved as 'Reduced_Subsystem_Enrichment.pdf'\n")
