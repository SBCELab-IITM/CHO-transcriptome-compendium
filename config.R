# Central paths for every R script. Edit here (or set the environment variables) - nowhere else.
# Run every script from the repository root, e.g.  Rscript 3_cellline_DEG_vs_K1/Cellline_combined.R
DATA_DIR <- Sys.getenv("CHO_DATA_DIR", "data")      # input files (counts, metadata, gene lists)
OUT_DIR  <- Sys.getenv("CHO_OUT_DIR",  "results")   # everything the scripts write
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
