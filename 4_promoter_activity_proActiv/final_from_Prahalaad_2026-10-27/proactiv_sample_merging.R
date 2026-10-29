# Paths come from config.R in the repository root; run this script from the repository root.
source("config.R")
#######################################
# CONVERT AND MERGE PROACTIV RDS FILES
#######################################
#
# PURPOSE
#   Convert per-sample proActiv RDS files into CSV files and create one
#   merged promoter-count table for downstream analysis.
#
# INPUT
#   - One RDS file per sample, containing a `promoterCounts` assay and
#     promoter annotations in `rowData`.
#   - The input directory is set by `input_dir`.
#
# PROCESS
#   1. Read each RDS file and extract its promoter counts and annotations.
#   2. Align counts with annotations and save a separate CSV for each sample.
#   3. Check that promoter annotations agree across samples.
#   4. Merge sample counts by `promoterId.
#
# OUTPUT
#   - One CSV per sample in `individual_dir`.
#   - `merged_promoter_raw_counts.csv` in `output_dir`, with annotation
#     columns and one `<sample_id>.raw_counts` column per sample.
#
# NOTE
#   Promoters absent from a sample are represented as missing values (NA)
#   in the merged table. Review this missing-value handling before using
#   the merged counts for downstream analysis.
#
# SETUP
#   Update `input_dir` and `output_dir` for your environment. The input
#   directory must contain the RDS files to convert.
#
#######################################

library(SummarizedExperiment)

input_dir <- DATA_DIR
output_dir <- file.path(OUT_DIR, "csv_output")
individual_dir <- file.path(output_dir, "individual_samples")

dir.create(individual_dir, recursive = TRUE, showWarnings = FALSE)

rds_files <- list.files(
  input_dir,
  pattern = "\\.rds$",
  full.names = TRUE
)

# Test with up to the first 3 RDS files
#rds_files <- head(rds_files, 3)

if (length(rds_files) == 0) {
  stop("No .rds files found in input_dir. Check that the path contains the RDS files.")
}

sample_tables <- lapply(rds_files, function(rdsfile) {
  obj <- readRDS(rdsfile)
  
  counts <- SummarizedExperiment::assay(obj, "promoterCounts")
  rowdata <- as.data.frame(SummarizedExperiment::rowData(obj))
  
  if (is.null(rownames(counts)) || is.null(rownames(rowdata))) {
    stop("Missing row names needed to align promoter counts and annotations: ", rdsfile)
  }
  
  if (anyDuplicated(rownames(counts)) || anyDuplicated(rownames(rowdata))) {
    stop("Duplicate row names found; cannot safely align counts and annotations: ", rdsfile)
  }
  
  # Align annotation rows to the count matrix.
  row_index <- match(rownames(counts), rownames(rowdata))
  if (anyNA(row_index) || nrow(counts) != nrow(rowdata)) {
    stop("Promoter rows do not match between counts and annotations: ", rdsfile)
  }
  rowdata <- rowdata[row_index, , drop = FALSE]
  
  # Use promoterId if present; otherwise use the aligned row names.
  if (!"promoterId" %in% names(rowdata)) {
    rowdata$promoterId <- rownames(rowdata)
  }
  if (anyNA(rowdata$promoterId) || anyDuplicated(rowdata$promoterId)) {
    stop("promoterId values must be present and unique: ", rdsfile)
  }
  
  # Convert list-valued annotations to CSV-friendly strings.
  list_cols <- vapply(rowdata, is.list, logical(1))
  rowdata[list_cols] <- lapply(rowdata[list_cols], function(column) {
    vapply(column, function(value) {
      paste(as.character(unlist(value)), collapse = ";")
    }, character(1))
  })
  
  if (ncol(counts) != 1) {
    stop("Expected one sample per RDS file, but found ", ncol(counts), ": ", rdsfile)
  }
  
  sample_id <- tools::file_path_sans_ext(basename(rdsfile))
  count_col <- paste0(sample_id, ".raw_counts")
  rowdata[[count_col]] <- as.numeric(counts[, 1])
  
  sample_csv <- file.path(individual_dir, paste0(sample_id, ".csv"))
  write.csv(rowdata, sample_csv, row.names = FALSE)
  
  message("Saved per-sample CSV: ", sample_csv)
  rowdata
})

# Exclude sample-specific proActiv fields (e.g. .mean, .gene.mean, .class) from shared annotations.
sample_ids <- tools::file_path_sans_ext(basename(rds_files))
sample_prefixes <- paste0(sample_ids, ".")

get_annotation_cols <- function(table) {
  is_sample_specific <- vapply(
    names(table),
    function(column) any(startsWith(column, sample_prefixes)),
    logical(1)
  )
  
  names(table)[!is_sample_specific]
}

annotation_cols <- get_annotation_cols(sample_tables[[1]])

for (i in seq_along(sample_tables)) {
  current_annotations <- get_annotation_cols(sample_tables[[i]])
  
  if (!setequal(annotation_cols, current_annotations)) {
    stop(
      "Shared annotation columns differ in ",
      basename(rds_files[[i]]),
      ". Check the columns before merging."
    )
  }
}

all_annotations <- do.call(rbind, lapply(sample_tables, function(table) {
  table[, annotation_cols, drop = FALSE]
}))

# Check that a promoter's annotations agree across samples.
rows_by_promoter <- split(seq_len(nrow(all_annotations)), all_annotations$promoterId)
has_conflicts <- vapply(rows_by_promoter, function(rows) {
  nrow(unique(all_annotations[rows, , drop = FALSE])) > 1
}, logical(1))

if (any(has_conflicts)) {
  stop("Annotation values differ between samples for one or more promoterId values.")
}

merged <- all_annotations[!duplicated(all_annotations$promoterId), , drop = FALSE]
promoter_order <- merged$promoterId

# Add each sample's counts, matching rows by promoterId.
for (i in seq_along(sample_tables)) {
  sample_id <- tools::file_path_sans_ext(basename(rds_files[[i]]))
  count_col <- paste0(sample_id, ".raw_counts")
  count_data <- sample_tables[[i]][, c("promoterId", count_col), drop = FALSE]
  
  merged <- merge(
    merged,
    count_data,
    by = "promoterId",
    all = TRUE,
    sort = FALSE
  )
}

# Preserve the original promoter order and put annotations before sample counts.
merged <- merged[match(promoter_order, merged$promoterId), , drop = FALSE]
count_cols <- paste0(
  tools::file_path_sans_ext(basename(rds_files)),
  ".raw_counts"
)
merged <- merged[, c(annotation_cols, count_cols), drop = FALSE]

merged_csv <- file.path(output_dir, "merged_promoter_raw_counts.csv")
write.csv(merged, merged_csv, row.names = FALSE)

message("Saved merged CSV: ", merged_csv)
message("Merged table dimensions: ", nrow(merged), " rows x ", ncol(merged), " columns")
