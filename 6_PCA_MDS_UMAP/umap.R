rm(list = ls())
# Paths come from config.R in the repository root; run this script from the repository root.
source("config.R")
set.seed(42)   # UMAP is stochastic; fixed for reproducibility
library(readr)
library(RColorBrewer)
data <- as.matrix(read.table(file.path(DATA_DIR, "Norm_file.csv"), header = TRUE,sep = ",",row.names = 1, as.is = TRUE))

library(readr)
library(uwot)
library("Matrix")
library("irlba")

Counts_t<-t(log2(data+1))
umap_results <- umap(Counts_t, min_dist = 0.1, n_neighbors = 15, n_epochs = 1000, metric = "euclidean", learning_rate = 0.5, init = "spectral", verbose = TRUE, n_threads = 12)
leg <- as.matrix(read.table(file.path(DATA_DIR, "metadata_new.txt"), header = TRUE,sep = "\t", row.names = 1, as.is = TRUE))
umap_results <- as.data.frame(cbind(umap_results, leg))
colnames(umap_results) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7")
head(umap_results)
#umap_results$V1=as.numeric(levels(umap_results$V1))[umap_results$V1]
#umap_results$V2=as.numeric(levels(umap_results$V2))[umap_results$V2]
umap_results$V1=as.numeric(umap_results$V1)
umap_results$V2=as.numeric(umap_results$V2)
library(ggplot2)
 
ggplot(umap_results, aes(V1, V2)) + geom_point(aes(colour = factor(V7)),alpha = 0.6, size=2) + scale_color_manual(labels = c('K1', 'K1 (GS)','DG44','S','Brain','Kidney','Liver','Lung','Ovary/Uterus','Spleen','Reproductive Tract','Heart','Muscle','Pancreas','DXB11'), values = c('#9f4b39','#a77c86','#fd8d3c','#6a51a3','#238b45','#cb181d','#2171b5','#f6db2f','#bdb76b','#98bdf0','black','violet','purple','pink','#FF7373')) + labs(x = paste0("UMAP1"), y = paste0("UMAP2"), colour = "Batch\n", shape = "Cell line\n") + theme_bw()  + guides(col = guide_legend(nrow = 16)) 

write.csv(umap_results,file.path(OUT_DIR, "UMAP_results.csv"))



#tsne


library(Rtsne)
tsne_results <- Rtsne(Counts_t, is_distance = FALSE, dim = 2, pca = TRUE, perplexity = 30, max_iter = 1000, num_threads = 12, verbose = TRUE)
tsne_results <- as.data.frame(cbind(tsne_results$Y, leg))
#tsne_results$V1 = as.numeric(levels(tsne_results$V1))[tsne_results$V1]
#tsne_results$V2 = as.numeric(levels(tsne_results$V2))[tsne_results$V2]
ggplot(tsne_results, aes(V1, V2)) + geom_point(aes(colour = factor(Cell_line)), size=2, alpha = 0.3)
                                                 