#MDS plot
rm(list = ls())
# Paths come from config.R in the repository root; run this script from the repository root.
source("config.R")
library(readr)
library(edgeR)

data<- as.matrix(read.table(file.path(DATA_DIR, "Combined_Filtered_New_Counts_nonzero.csv"), header = TRUE,sep = ",", row.names = 1, as.is = TRUE))
library(magrittr)
library(readr)
Counts <- as.matrix(read.table(file.path(DATA_DIR, "Final_Metadata.csv"), header = TRUE,sep = ",", row.names = 1, as.is = TRUE))
Counts_t<-t(log2(data+1))
alpha<-t(Counts_t)
library("factoextra")
dist.cor <- get_dist(Counts_t, method = "pearson")
fit <- cmdscale(dist.cor,eig=TRUE, k=34)
eigvlaues <- as.data.frame(fit[["eig"]])
pc1<-paste0(round(eigvlaues[1,1]/sum(eigvlaues)*100,0),"%")
pc2<-paste0(round(eigvlaues[2,1]/sum(eigvlaues)*100,0),"%")
pc3<-paste0(round(eigvlaues[3,1]/sum(eigvlaues)*100,0),"%")
library(MASS)
initial_config <- cmdscale(dist.cor, k = 34)
mds_result <- isoMDS(dist.cor, y = cmdscale(dist.cor, 34), k = 34, maxit = 200, trace = TRUE, tol = 1e-3, p = 2)
mds_fit <- as.data.frame(mds_result[["points"]])
pcs <- read_csv(file.path(DATA_DIR, "new_count_tmm_pcs.csv"))
library(ggplot2)
head(pcs)
pcs1 <- as.data.frame(cbind(mds_fit$V1,mds_fit$V2,pcs$Media,pcs$Cel_Line,pcs$Prod,pcs$Phase))
#pcs1$V1 = as.numeric(levels(pcs1$V1))[pcs1$V1]
#pcs1$V2 = as.numeric(levels(pcs1$V2))[pcs1$V2]
pcs1$V1=as.numeric(pcs1$V1)
pcs1$V2=as.numeric(pcs1$V2)

str(pcs) #checking the string of pcs
colnames(pcs) #checking the colnames of pcs


# Corrected version with a missing comma added
p <- ggplot(pcs, aes(x = PC1, y = PC2)) + 
  geom_point(aes(colour = factor(Batch)), size = 2) + 
  labs(x = "Principal Component 1", 
       y = "Principal Component 2", 
       colour = "Batch") + 
  theme_bw() +
  scale_color_manual(labels = c('1' ,'2','3' ,'4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'),
                     values = c('#f6db2f', '#2171b5','#fd8d3c','#6a51a3','violet','black','pink','#98bdf0','#009999',"#0000FF","#ABB065","#ACA4E2","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
                                "#A65628","#F781BF","#999999","#66C2A5","#FFD92F"))

print(p)


p <- ggplot(pcs, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = Media), size = 2) +
  labs(x = "Principal Component 1",
       y = "Principal Component 2",
       colour = "Media") +
  theme_bw() +
  scale_color_manual(
    values = c("CD-CHO-SFM" = "#f6db2f",
               "EXCELL CD-CHO (SAFC)" = "#2171b5",
               "CD-CHO (Gibco)" = "#fd8d3c",
               "CDM4CHO (HyClone)" = "#6a51a3",
               "PowerCHO (Lonza)" = "#cb181d",
               "CDFortiCHO" ="violet",
               "G13 medium" ="black",
               "DA1 medium" ="pink",
               "ActiPro (Gibco)" ="#98bdf0",
               "NM" = "#009999",
               "80%DMEM+10%serum+20%CD-CHO" ="#0000FF",
               "CHO-S-SFM II (Gibco)" ="#ABB065",
               "CD-CHO (Gibco) + PF-CHO (HyClone)" = "#ACA4E2"
               )
  )

print(p)

p <- ggplot(pcs, aes(x = PC1, y = PC2)) + 
  geom_point(aes(colour = factor(Batch)), size = 2) + 
  labs(x = "Principal Component 1", 
       y = "Principal Component 2", 
       colour = "Batch") + 
  theme_bw() +
  scale_color_manual(
    labels = c('1','2','3','4','5','6','7','8','9','10',
               '11','12','13','14','15','16','17','18','19','20','21','22'),
    values = c('#f6db2f', '#2171b5','#fd8d3c','#6a51a3','violet','black','pink',
               '#98bdf0','#009999',"#0000FF","#ABB065","#ACA4E2",
               "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
               "#A65628","#F781BF","#999999","#66C2A5","#FFD92F")
  )

print(p)


unique(pcs$Media)
table(pcs$Media)





p <- ggplot(pcs, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = factor(Batch)), size = 2) +
  labs(x = "Principal Component 1",
       y = "Principal Component 2",
       colour = "Batch") +
  theme_bw() +
  scale_color_manual(
    labels = c('1','2','3','4','5','6','7','8','9','10',
               '11','12','13','14','15','16','17','18','19','20','21','22'),
    values = c('#f6db2f', '#2171b5','#fd8d3c','#6a51a3','violet','black','pink',
               '#98bdf0','#009999',"#0000FF","#ABB065","#ACA4E2",
               "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
               "#A65628","#F781BF","#999999","#66C2A5","#FFD92F")
  )

print(p)












#Plotting K1 based on Media
p <- ggplot(pcs, aes(x = PC1, y = PC2)) + 
  geom_point(aes(colour = factor(Media)), size = 2) + 
  labs(x = "Principal Component 1", 
       y = "Principal Component 2", 
       colour = "Media") + 
  theme_bw()
#+ scale_color_manual(labels = c( 'CHO_S_SFM_II', 'ExCELL_CD-CHO' ,'CDFortiCHO' ,'CD_CHO', 'CDM4CHO' , 'ActiPro','G13' ,'DA1' ), values = c('#fd8d3c','#6a51a3''#bdb76b','#98bdf0','black','violet','purple','pink'))









ggplot(pcs1, aes(V1, V2)) + geom_point(aes(colour = factor(V4)), size=2)  + labs(x = paste0("Coordinate 1 (",pc1,")"), y = paste0("Coordinate 2 (",pc2,")"), colour = "Cell line\n") + theme_bw()
#+ scale_color_manual(labels = c('K1', 'P' ,'NP' ), values = c('#9f4b39' , '#a77c86' ,'#fd8d3c'))
ggplot(pcs, aes(V1, V2)) + geom_point(aes(colour = factor(V7), shape = factor(V5)), size=2) + labs(x = paste0("Coordinate 1 (",pc1,")"), y = paste0("Coordinate 2 (",pc2,")"), colour = "Prod\n", shape = "Cell line\n") + theme_bw()
cumpro <- cumsum(eigvlaues / sum(eigvlaues))
plot(cumpro[1:15,], xlab = "Coordinate #", ylab = "Amount of explained variance", main = "Cumulative variance plot")

mds_fit <- as.data.frame(cbind(mds_fit$V1,mds_fit$V2,mds_fit$V3))
write.csv(mds_fit,file.path(OUT_DIR, "MDS_results.csv"))
