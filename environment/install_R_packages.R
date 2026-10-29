# Installs the R packages used across the scripts. Run once: Rscript environment/install_R_packages.R
# Record exact versions afterwards with renv::snapshot() (recommended before the Zenodo release).
cran <- c("tidyverse","readr","readxl","dplyr","tidyr","stringr","ggplot2","ggraph","magrittr",
          "MASS","Matrix","RColorBrewer","Rtsne","irlba","uwot","factoextra","VennDiagram","openxlsx","writexl","BiocManager")
cran <- cran[nzchar(cran)]
bioc <- c("DESeq2","edgeR","limma","SummarizedExperiment","clusterProfiler","enrichplot","AnnotationHub",
          "biomaRt","rtracklayer","proActiv")
install.packages(setdiff(cran, rownames(installed.packages())))
BiocManager::install(setdiff(bioc, rownames(installed.packages())), ask = FALSE)
# Genome / annotation packages used by proActiv and promoter scripts (see script headers):
# BSgenome.Cgriseus.NCBI.CriGriPICR, EnsDb.Cgriseus.v1 (install if available from your Bioconductor release).
