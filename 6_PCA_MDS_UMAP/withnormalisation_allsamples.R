rm(list = ls())
# Paths come from config.R in the repository root; run this script from the repository root.
source("config.R")
library(readr)
library(edgeR)
Counts <- as.matrix(read.table(file.path(DATA_DIR, "merged_processed1.csv"), header = TRUE,sep = ",", row.names = 1, as.is = TRUE))

coldata <- data.frame(Cell=c("K1", "K1", "K1", "K1", "K1", "K1", "K1", "K1", "K1", "K1", "K1", "K1", "K1", "K1", "K1", "K1", "K1", "K1", "K1", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "K1", "K1", "K1", "K1", "K1", "K1", "K1", "K1", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "K1", "K1", "K1", "K1", "K1", "K1", "Brain", "Brain", "ReproductiveTract", "Heart", "Kidney", "Liver", "Lung", "Muscle", "Pancreas", "Spleen", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "K1", "K1", "DG44", "DG44", "DG44", "DG44", "DG44", "DG44", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1 (GS)", "K1", "K1", "Brain", "Kidney", "Liver", "Lung", "Ovary", "Spleen", "K1", "K1", "K1", "K1", "K1", "K1", "K1", "K1", "K1", "K1", "K1", "K1", "K1", "K1", "K1" , "DXB11", "DXB11" , "DXB11" , "DXB11" , "DXB11" , "DXB11" , "DXB11" , "DXB11" , "DXB11" , "DXB11" , "DXB11" , "DXB11" , "DXB11" , "DXB11" , " K1" , "K1" , "K1" , "K1" , "K1" , "K1" , "K1" , "K1"),  Batch = c("3", "3", "3", "3", "3", "3", "3", "3", "3", "10", "10", "10", "10", "10", "10", "10", "10", "10", "10", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "14", "7", "7", "7", "7", "7", "7", "7", "7", "7", "7", "7", "7", "7", "7", "7", "8", "8", "8", "8", "8", "8", "8", "8", "8", "8", "8", "8", "8", "8", "8", "8", "8", "8", "8", "8", "8", "8", "1", "1", "1", "1", "1", "1", "1", "1", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "13", "13", "13", "13", "13", "13", "12", "12", "12", "12", "12", "12", "12", "12", "12", "12", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "4", "4", "4", "4", "4", "4", "4", "4", "9", "9", "9", "9", "9", "9", "9", "9", "9", "9", "9", "9", "11", "11", "11", "11", "11", "11", "11", "11", "15", "15", "15", "15", "15", "15", "15", "15", "15", "15", "15", "15", "15", "15", "15", "16" , "16" , "16" , "16" , "16" , "16" , "16" , "16" , "16" , "16" , "16" , "16" , "16" , "16", "17" , "17" , "17" , "17" , "17" , "17" , "17" , "17")) 
#Media = c("A", "A", "A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B", "C", "C", "C", "C", "C", "C", "C", "C", "C", "D", "D", "E", "E", "E", "E", "E", "E", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "F", "F", "G", "G", "G", "G", "G", "G", "C", "C" , "C" , "C" , "C" , "C" , "C" , "C" , "C" , "C" , "C" , "C" , "C" , "C" ,"A" , "A" , "A" , "A" , "A" , "A" , "A" , "A")
design.mat <- model.matrix(~Cell+Batch,coldata)
dgList <- DGEList(counts=Counts, genes=rownames(Counts))

tcc <- calcNormFactors(dgList, method="TMM", refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75)
tcc <- estimateCommonDisp(tcc, design.mat)
norm_counts <- tcc[["pseudo.counts"]]
write.csv(norm_counts,file.path(OUT_DIR, "Norm_clean1.csv"))
                     
Counts_t<-t(log2(norm_counts+1))
alpha<-scale(Counts_t)
alphadot<-t(alpha)
pc<-prcomp(na.omit(alphadot), retx = TRUE)
pcs<-pc$rotation
pc1_MA<-paste0(round(summary(pc)$importance[2,1]*100,0),"%")
pc2_MA<-paste0(round(summary(pc)$importance[2,2]*100,0),"%")
pc3_MA<-paste0(round(summary(pc)$importance[2,3]*100,0),"%")
summary(pc)
cumpro <- cumsum(pc$sdev^2 / sum(pc$sdev^2))
plot(cumpro[0:15], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
write.csv(pcs,file.path(OUT_DIR, "mergedcounts_cleaned_pcs.csv"));
                      
library(readr)
pcs <- read_csv(file.path(OUT_DIR, "mergedcounts_cleaned_pcs.csv"))
library(ggplot2)
ggplot(pcs, aes(pcs$PC1, pcs$PC2)) + geom_point(aes(colour = factor(pcs$Media), shape = factor(pcs$Cell)), size=2) + scale_color_brewer(palette="Set2") + labs(x = paste0("PC1 (",pc1_MA,")"), y = paste0("PC2 (",pc2_MA,")"), colour = "Media\n", shape = "Cell line\n") + theme_bw()
ggplot(pcs, aes(pcs$PC1, pcs$PC3)) + geom_point(aes(colour = factor(pcs$Media), shape = factor(pcs$Cell)), size=2) + labs(x = paste0("PC1 (",pc1_MA,")"), y = paste0("PC2 (",pc3_MA,")"), colour = "Media\n", shape = "Cell line\n") + theme_bw()
ggplot(pcs, aes(pcs$PC1, pcs$PC2)) + geom_point(aes(colour = factor(pcs$Batch), shape = factor(pcs$Cell)), size=2) + labs(x = paste0("PC1 (",pc1_MA,")"), y = paste0("PC2 (",pc2_MA,")"), colour = "Batch\n", shape = "Cell line\n") + theme_bw()  + guides(col = guide_legend(nrow = 4)) 
ggplot(pcs, aes(pcs$PC1, pcs$PC3)) + geom_point(aes(colour = factor(pcs$Batch), shape = factor(pcs$Cell)), size=2) + labs(x = paste0("PC1 (",pc1_MA,")"), y = paste0("PC2 (",pc3_MA,")"), colour = "Batch\n", shape = "Cell line\n") + theme_bw()  + guides(col = guide_legend(nrow = 4))
ggplot(pcs, aes(pcs$PC1, pcs$PC2)) + geom_point(aes(colour = factor(pcs$Cell)), size=2) + scale_color_manual(labels = c('K1', 'K1 (GS)','DG44','S','Brain','Kidney','Liver','Lung','Ovary/Utreus','Spleen','Reproductive Tract','Heart','Muscle','Pancreas','DXB11'), values = c('#9f4b39','#a77c86','#fd8d3c','#6a51a3','#238b45','#cb181d','#2171b5','#f6db2f','#bdb76b','#98bdf0','black','violet','purple','pink','#FF7373')) + labs(x = paste0("PC1 (",pc1_MA,")"), y = paste0("PC2 (",pc2_MA,")"), colour = "Batch\n", shape = "Cell line\n") + theme_bw()  + guides(col = guide_legend(nrow = 16)) 

data <- as.matrix(read.table(file.path(DATA_DIR, "Norm_file.csv"), header = TRUE, sep = ",", row.names = 1, as.is = TRUE))
norm_counts <- data
library(magrittr)
library(readr)
t_Counts <- t(log2(norm_counts+1))
install.packages("factoextra")
library("factoextra")
dist.cor <- get_dist(Counts_t, method = "pearson")
fit <- cmdscale(dist.cor,eig=TRUE, k=50)
eigvlaues <- as.data.frame(fit[["eig"]])
pc1<-paste0(round(eigvlaues[1,1]/sum(eigvlaues)*100,0),"%")
pc2<-paste0(round(eigvlaues[2,1]/sum(eigvlaues)*100,0),"%")
pc3<-paste0(round(eigvlaues[3,1]/sum(eigvlaues)*100,0),"%")
library(MASS)
mds_result <- isoMDS(dist.cor, y = cmdscale(dist.cor, 50), k = 50, maxit = 100, trace = TRUE, tol = 1e-3, p = 2)
mds_fit <- as.data.frame(mds_result[["points"]])
pcs <- read_csv(file.path(OUT_DIR, "mergedcounts_cleaned_pcs.csv"))
library(ggplot2)
pcs1 <- as.data.frame(cbind(mds_fit$V1,mds_fit$V2,pcs$Media,pcs$Cell,pcs$Prod,pcs$Phase,pcs$Batch))
pcs1$V1 = as.numeric(levels(pcs1$V1))[pcs1$V1]
pcs1$V2 = as.numeric(levels(pcs1$V2))[pcs1$V2]
pcs1$V1=as.numeric(pcs1$V1)
pcs1$V2=as.numeric(pcs1$V2)
ggplot(pcs1, aes(V1, V2)) + geom_point(aes(colour = factor(V4)), size=1)  + labs(x = paste0("Coordinate 1 (",pc1,")"), y = paste0("Coordinate 2 (",pc2,")"), colour = "Cell line\n") + theme_bw()
#+ scale_color_manual(labels = c('K1', 'K1 (GS)','DG44','S','Brain','Kidney','Liver','Lung','Ovary/Utreus','Spleen','Reproductive Tract','Heart','Muscle','Pancreas'), values = c('#9f4b39','#a77c86','#fd8d3c','#6a51a3','#238b45','#cb181d','#2171b5','#f6db2f','#bdb76b','#98bdf0','black','violet','purple','pink'))
ggplot(pcs1, aes(V1, V2)) + geom_point(aes(colour = factor(V3), shape = factor(V4)), size=2) + labs(x = paste0("Coordinate 1 (",pc1,")"), y = paste0("Coordinate 2 (",pc2,")"), colour = "Media\n", shape = "Cell line\n") + theme_bw()
cumpro <- cumsum(eigvlaues / sum(eigvlaues))
plot(cumpro[1:15,], xlab = "Coordinate #", ylab = "Amount of explained variance", main = "Cumulative variance plot")

mds_fit <- as.data.frame(cbind(mds_fit$V1,mds_fit$V2,mds_fit$V3))
write.csv(mds_fit,file.path(OUT_DIR, "MDS_results.csv"))



  