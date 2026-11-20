# TMM normalisation, MDS, UMAP, heatmap and variancePartition for CHO samples
# Final analysis script (SBCE lab, IIT Madras), taken unchanged from the lab's combined script.
# NOTE: input paths are hard-coded to the original machine; edit them to your local copies of the
# count matrices and metadata (Supplementary Table 1) before running.

rm(list = ls())
library(readr)
library(edgeR)
Counts <- as.matrix(read.table("C:/Users/Meiyappan/OneDrive/Documents/CHO_rna_new/Latest/Count_cell_line_Updated_06092025.csv", header = TRUE,sep = ",", row.names = 1, as.is = TRUE))
coldata <- as.data.frame(read.table("C:/Users/Meiyappan/OneDrive/Documents/CHO_rna_new/Latest/Final_Metadata_Cell_line_06092025.txt", header = TRUE,sep = "\t", row.names = 1, as.is = TRUE))
coldata$Study <- factor(coldata$Study)
coldata$Layout <- factor(coldata$Layout)
coldata$Cell_line <- factor(coldata$Cell_line)
coldata$Producer <- factor(coldata$Producer)
coldata$Product <- factor(coldata$Product)
coldata$Batch <- factor(coldata$Batch)
coldata$Media <- factor(coldata$Media)
coldata$Library <- factor(coldata$Library)
design.mat <- model.matrix(~Cell_line+Study+Producer+Product+Batch+Media+Layout+Library,coldata)
dgList <- DGEList(counts=Counts, genes=rownames(Counts))
tcc <- calcNormFactors(dgList, method="TMM", refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75)
tcc <- estimateCommonDisp(tcc, design.mat)
norm_counts <- tcc[["pseudo.counts"]]
library(magrittr)
library(readr)
t_Counts <- t(log2(norm_counts+1))
library("factoextra")
dist.cor <- get_dist(t_Counts, method = "pearson")
fit <- cmdscale(dist.cor,eig=TRUE, k=50)
eigvlaues <- as.data.frame(fit[["eig"]])
pc1<-paste0(round(eigvlaues[1,1]/sum(eigvlaues)*100,0),"%")
pc2<-paste0(round(eigvlaues[2,1]/sum(eigvlaues)*100,0),"%")
pc3<-paste0(round(eigvlaues[3,1]/sum(eigvlaues)*100,0),"%")
library(MASS)
mds_result <- isoMDS(dist.cor, y = cmdscale(dist.cor, 50), k = 50, maxit = 350, trace = TRUE, tol = 1e-3, p = 2)
mds_fit <- as.data.frame(mds_result[["points"]])
leg <- as.matrix(read.table("C:/Users/Meiyappan/OneDrive/Documents/CHO_rna_new/Latest/Final_Metadata_Cell_line_06092025.txt", header = TRUE,sep = "\t", row.names = 1, as.is = TRUE))
mds_results <- as.data.frame(cbind(mds_fit$V1,mds_fit$V2,mds_fit$V3,leg))
mds_results$V1=as.numeric(mds_results$V1)
mds_results$V2=as.numeric(mds_results$V2)
mds_results$V3=as.numeric(mds_results$V3)
library(ggplot2)
ggplot(mds_results, aes(mds_results$V1, mds_results$V2)) + geom_point(aes(shape = factor(mds_results$Cell_line), colour = factor(mds_results$Media)), size=3) + labs(x = paste0("PC1 (",pc1,")"), y = paste0("PC2 (",pc2,")"), shape = "Cell line\n", colour = "Media\n") + theme_bw()
ggplot(mds_results, aes(mds_results$V1, mds_results$V2)) + geom_point() + geom_text(label=mds_results$Batch)

library(uwot)
library(RColorBrewer)
t_Counts <- t(log2(norm_counts+1))
umap_results <- umap(t_Counts, min_dist = 0.1, n_neighbors = 70, n_epochs = 1000, metric = "euclidean", learning_rate = 0.5, init = "spectral", verbose = TRUE, n_threads = 12)
leg <- as.matrix(read.table("C:/Users/Meiyappan/OneDrive/Documents/CHO_rna_new/Latest/Final_Metadata_Cell_line_06092025.txt", header = TRUE,sep = "\t", row.names = 1, as.is = TRUE))
umap_results <- as.data.frame(cbind(umap_results, leg))
umap_results$V1=as.numeric(umap_results$V1)
umap_results$V2=as.numeric(umap_results$V2)
library(ggplot2)
ggplot(umap_results, aes(umap_results$V1, umap_results$V2)) + geom_point(aes(shape = factor(umap_results$Cell_line), colour = factor(umap_results$Media)), size=3) + labs(x = paste0("UMAP1"), y = paste0("UMAP2"), shape = "Cell line\n", colour = "Media\n") + theme_bw()
ggplot(umap_results, aes(umap_results$V1, umap_results$V2)) + geom_point() + geom_text(label=umap_results$Batch)

library(readr)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
t_Counts <- log2(norm_counts+1)
makeColorRampPalette <- function(colors, cutoff.fraction1, cutoff.fraction2, cutoff.fraction3, cutoff.fraction4, num.colors.in.palette)
{
  stopifnot(length(colors) == 8)
  ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction1)
  ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (cutoff.fraction2 - cutoff.fraction1))
  ramp3 <- colorRampPalette(colors[4:5])(num.colors.in.palette * (cutoff.fraction3 - cutoff.fraction2))
  ramp4 <- colorRampPalette(colors[5:6])(num.colors.in.palette * (cutoff.fraction4 - cutoff.fraction3))
  ramp5 <- colorRampPalette(colors[6:7])(num.colors.in.palette * (1 - cutoff.fraction4))
  return(c(ramp1, ramp2, ramp3, ramp4, ramp5))
}
cutoff.distance1 <- 1
cutoff.distance2 <- 5
cutoff.distance3 <- 7
cutoff.distance4 <- 11
cols <- makeColorRampPalette(c("white","#c7e8ec","#c7e8ec","#45b4c2","#45b4c2","#16507c","#16507c","#051440"), cutoff.distance1 / 21, cutoff.distance2 / 21, cutoff.distance3 / 21, cutoff.distance4 / 21, 100)
my_sample_col <- as.data.frame(umap_results$Cell_line)
row.names(my_sample_col) <- colnames(t_Counts)
my_colour = list(sample = c(DG44 = "#4269d0", DXB11 = "#6cc5b0", K1 = "#efb118", S = "#ff725c"))
pheatmap(t_Counts, annotation_col = my_sample_col, annotation_colors = my_colour, na_col = "#f2f2f2", show_rownames = F, show_colnames = F, border_color = "black",  color = cols)

library(readr)
library(variancePartition)
library(limma)
library(edgeR)
t_Counts <- log2(norm_counts+1)
myinfo <- read_csv("C:/Users/Meiyappan/OneDrive/Documents/CHO_rna_new/Latest/Final_Metadata_Cell_line_06092025.csv")
form <- ~ Batch + (1|Phase) + (1|Cell_line) + (1|Producer) + (1|Media)  + (1|Library) + (1|Layout)
isexpr <- rowSums(cpm(t_Counts)>1) >= 0.01 * ncol(t_Counts)
gExpr <- t_Counts[isexpr,]
varPart <- fitExtractVarPartModel(gExpr, form, myinfo)
vp <- sortCols(varPart)
plotVarPart(vp)

i <- which.max(varPart$Cell_line)
i <- 8954
GE <- data.frame(Expression = gExpr[i, ], Phase = myinfo$Phase)
plotStratify(Expression ~ Phase, GE, main = rownames(gExpr)[i])
write.csv(varPart,"C:/Users/Meiyappan/OneDrive/Documents/CHO_rna_new/Latest/varpartnew.csv")

C <- canCorPairs(form, myinfo)
plotCorrMatrix(C)

