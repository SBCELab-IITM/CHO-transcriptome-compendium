# Paths come from config.R in the repository root; run this script from the repository root.
source("config.R")
# Gained / lost / reduced genes in CHO vs hamster tissues, for the 4 CHO subsets
# Subset1 = CHO-S, Subset2 = K1 + S, Subset3 = K1 + DXB11, Subset4 = DG44
# Gained/lost: expressed = value > 1 in >= 85% of samples.  Reduced: limma, logFC < -1 and adj.P < 0.05
# Inputs are in data/gain_loss_inputs/. Prints a summary table and writes results/Gain_loss_summary.csv

library(limma)

input_dir <- file.path(DATA_DIR, "gain_loss_inputs")

read_counts <- function(file, sep) {
  x <- read.delim(file.path(input_dir, file), sep = sep, check.names = FALSE)
  m <- as.matrix(x[, -1])
  rownames(m) <- x[[1]]
  m
}
base_name <- function(g) sub("_[0-9]+$", "", g)   # gene names carry _1, _2 ... suffixes

# Subset1 file has CHO + tissue columns, all other files have only the 16 CHO columns
both <- read_counts("Subset1_CHO_Tissue_counts.txt", "\t")
cho1 <- colnames(read_counts("Subset1_v1_Count.csv", ","))
tissue_cols <- setdiff(colnames(both), cho1)
ref_cho <- read_counts("Subset4_Count.csv", ",")      # same gene order as 'both'

classify <- function(mat, cutoff = 1, prop = 0.85) {
  n <- ncol(mat)
  ifelse(rowSums(mat > cutoff) / n >= prop, "Expressed",
         ifelse(rowSums(mat <= cutoff) / n >= prop, "Non-expressed", "Indeterminate"))
}

# put the rows of a CHO-only file into the same order as the tissue block
# (repeated gene names are matched by closest mean expression)
match_rows <- function(cho) {
  cho_key <- base_name(rownames(cho))
  tis_key <- base_name(rownames(both))
  ref <- rowMeans(log2(ref_cho + 1))
  mine <- rowMeans(log2(cho + 1))
  idx <- integer(nrow(cho))
  for (g in unique(cho_key)) {
    a <- which(cho_key == g)
    b <- which(tis_key == g)
    for (i in a) {
      k <- b[which.min(abs(ref[b] - mine[i]))]
      idx[i] <- k
      b <- setdiff(b, k)
    }
  }
  idx
}

run_subset <- function(cho, tissue) {
  expr <- cbind(cho, tissue)
  group <- factor(c(rep("CHO", ncol(cho)), rep("Tissue", ncol(tissue))))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  fit <- lmFit(expr, design)
  fit <- eBayes(contrasts.fit(fit, makeContrasts(CHOvsTissue = CHO - Tissue, levels = design)))
  deg <- topTable(fit, coef = "CHOvsTissue", number = Inf)
  cho_class <- classify(cho)
  tis_class <- classify(tissue)
  c(Gained = sum(cho_class == "Expressed" & tis_class == "Non-expressed"),
    Lost = sum(cho_class == "Non-expressed" & tis_class == "Expressed"),
    Reduced = sum(deg$logFC < -1 & deg$adj.P.Val < 0.05))
}

files <- c(Subset1 = "Subset1_v1_Count.csv", Subset2 = "Subset2_Count.csv",
           Subset3 = "Subset3_Count.csv", Subset4 = "Subset4_Count.csv")
results <- NULL
for (s in names(files)) {
  cho <- read_counts(files[s], ",")
  if (s %in% c("Subset1", "Subset4")) {
    # same row order as the tissue file
    tissue <- both[, tissue_cols]
  } else {
    idx <- match_rows(cho)
    tissue <- both[idx, tissue_cols]
  }
  results <- rbind(results, run_subset(cho, tissue))
}
rownames(results) <- names(files)
results <- cbind(results, Workbook_Gained = c(12, 19, 58, 60), Workbook_Lost = c(1962, 1360, 1419, 619),
                 Workbook_Reduced = c(5090, 4174, 5140, 4420))
print(results)
write.csv(results, file.path(OUT_DIR, "Gain_loss_summary.csv"))
