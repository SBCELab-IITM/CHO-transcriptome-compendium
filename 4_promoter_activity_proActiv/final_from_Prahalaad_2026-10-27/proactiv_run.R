# Paths come from config.R in the repository root; run this script from the repository root.
source("config.R")
#######################################
# PROACTIV PROMOTER-COUNT GENERATION
#######################################
#
# PURPOSE
#   Generate promoter-level counts from BAM files using proActiv and the
#   Chinese hamster promoter annotation.
#
# INPUTS
#   - BAM files from the directory specified by `bam_files`
#   - A GTF file used to prepare the promoter annotation
#   - The BSgenome.Cgriseus.NCBI.CriGriPICR genome package
#
# PROCESS
#   1. Prepare the promoter annotation from the GTF file.
#   2. Run proActiv separately on each BAM file.
#   3. Remove rows with missing promoter counts.
#   4. Save each sample's result as an RDS file.
#
# OUTPUTS
#   One `<sample_name>.rds` file per BAM file, saved in `output_dir`.
#   These RDS files are inputs to the separate script that exports
#   per-sample CSVs and creates the merged promoter-count table.
#
# SETUP
#   Update the GTF, BAM input, and RDS output paths for your environment.
#   Ensure the output directory exists before running the script.
#
#######################################


# Load required library
library(proActiv)
library(BSgenome.Cgriseus.NCBI.CriGriPICR)

# Step 0:
# Add Cgriseus genomeinfo file to BSgenome
# Keep BSgenome.Cgriseus.NCBI.CriGriPICR folder in the working dir


# Step 1: Prepare the promoter annotation (object needs to be created only once, can use for all runs)
gtf.file <- file.path(DATA_DIR, "GCF_003668045.1_CriGri-PICR_genomic.gtf")
promoterAnnotation <- preparePromoterAnnotation(file = gtf.file, species = 'Cricetulus_griseus') # Chinese Hamster

# Step 2: Define file paths
bam_files <- list.files(
  file.path(DATA_DIR, "CHO_bam_files"),
  full.names = TRUE
)


# Step 3: Define output directory
output_dir <- file.path(OUT_DIR, "proactiv_rds_outputs")


# Step 4: Initialize counter
countfile <- 0

# Step 5: Loop through BAM files and process each
for (file in bam_files) {
  # Extract sample name (base name without extension)
  sample_name <- tools::file_path_sans_ext(basename(file))
  
  countfile <- countfile + 1
  
  # Set condition as the sample name (just to provide a name in the condition param while running proactiv)
  condition <- sample_name
  
  # Run the proActiv pipeline
  result <- proActiv(
    files = file,
    promoterAnnotation = promoterAnnotation,
    condition = condition,
    genome = 'BSgenome.Cgriseus.NCBI.CriGriPICR'
  )
  
  # Remove rows with missing promoter counts
  result <- result[complete.cases(assays(result)$promoterCounts), ]
  
  # Save the result as an .rds file
  output_file <- file.path(output_dir, paste0(sample_name, ".rds"))
  saveRDS(result, output_file)
  
  # Print progress
  cat("Processed and saved:", sample_name, countfile, "\n")
}

cat("All files processed and saved in", output_dir, "\n")


