setwd("~/data/CD63GFP/WD/")
output <- "../WD/"

source("../scripts/Functions/20170228_Standard_Bioinformatics_functions.R")
source("../scripts/Functions/20160428_tSNE_Stouffers_functions.R")
source("../scripts/Functions/20170228_tSNE_Graphing_functions.R")
source("../scripts/Functions/20170228_WGCNA_SCDE.R")
library(Rtsne)
library(fields)
library(Matrix)
library(igraph)
library(scales)
library(colorspace)
library(gplots)

rpkmdata<-read.table("../WD/S_30m_1x_2x_all_batch_correction_norm_express.txt", header=T)
data <- t(rpkmdata)

SpleenGenes <- read.csv("../WD/SpleenBcell100GenesEnsembl.txt")
Spleen <- rep(0, ncol(data))
names(Spleen) <- colnames(data)
Spleen[which(names(Spleen) %in% SpleenGenes$ensembl_gene_id)] <- 100
SpleenBcell <- Spleen

SpleenGenes <- read.csv("../WD/SpleenMacrophage100GenesEnsembl.txt")
Spleen <- rep(0, ncol(data))
names(Spleen) <- colnames(data)
Spleen[which(names(Spleen) %in% SpleenGenes$ensembl_gene_id)] <- 100
SpleenMacrophage <- Spleen

SpleenGenes <- read.csv("../WD/SpleenTcell100GenesEnsembl.txt")
Spleen <- rep(0, ncol(data))
names(Spleen) <- colnames(data)
Spleen[which(names(Spleen) %in% SpleenGenes$ensembl_gene_id)] <- 100
SpleenTcell <- Spleen

MitoGenes <- read.table("../WD/mouse_mitochondrial_genes_list_annotation.txt", sep = "\t", header = T)
Spleen <- rep(0, ncol(data))
names(Spleen) <- colnames(data)
Spleen[which(names(Spleen) %in% MitoGenes$ensembl_gene_id)] <- 100
Mitochondrial <- Spleen

MetaSpleen <- as.data.frame(matrix(NA, 4, ncol(data)))
MetaSpleen[1,] <- as.numeric(SpleenBcell)
MetaSpleen[2,] <- as.numeric(SpleenMacrophage)
MetaSpleen[3,] <- as.numeric(SpleenTcell)
MetaSpleen[4,] <- as.numeric(Mitochondrial)
colnames(MetaSpleen) <- colnames(data)
rownames(MetaSpleen) <- c("Bcell","Macrophage","Tcell","Mitochondrial")
MetaSpleen <- t(MetaSpleen)

################################################

rpkmdata<-read.table("../WD/S_30m_1x_2x_all_batch_correction_norm_express.txt", header=T)

S_HEK_30m_1x <- rpkmdata[,grep('^S_HEK_30m_1x', colnames(rpkmdata))]
#S_HEK_30m_2x <- rpkmdata[,grep('^S_HEK_30m_2x', colnames(rpkmdata))]
#S_HEK_3h_2x <- rpkmdata[,grep('^S_HEK_3h_2x', colnames(rpkmdata))]
S_CAP_30m_1x <- rpkmdata[,grep('^S_CAP_30m_1x', colnames(rpkmdata))]
#S_CAP_30m_2x <- rpkmdata[,grep('^S_CAP_30m_2x', colnames(rpkmdata))]
#S_CAP_3h_2x <- rpkmdata[,grep('^S_CAP_3h_2x', colnames(rpkmdata))]

rpkmdata <- cbind(S_HEK_30m_1x, S_CAP_30m_1x)
data <- cbind(MetaSpleen, S_HEK_30m_1x, S_CAP_30m_1x)

data <- t(as.data.frame(data))
data[data==0] <- 0.0001
data <- log2(data)

#var.rpkms<-var.genes(data)
#data <- subset(data, select = colnames(data) %in% var.rpkms)

genes <- read.csv("../WD/Spleen100TabulaGenesEnsembl.txt", header = F)
#genes <- c(genes$ensembl_gene_id, MitoGenes$ensembl_gene_id)
data <- subset(data, select = colnames(data) %in% genes$V1)

CellsCor=(cor(t(data)))
matrix<-as.matrix(CellsCor)

SpleenBcell_cor <- CellsCor[1,(5:616)]

#SpleenBcell.level <- cRamp.orig(SpleenBcell_cor, colors = c("white","white","red","darkred"))
SpleenBcell.level <- as.factor(findInterval(SpleenBcell_cor, c(min(SpleenBcell_cor),0.05,0.1,0.2,max(SpleenBcell_cor)), rightmost.closed = FALSE, all.inside = FALSE))
SpleenBcell.level <- c("white","pink","red","darkred")[SpleenBcell.level] #c("white","greenyellow", "green3", "darkgreen")
SpleenBcell.level.dark <- col2rgb(SpleenBcell.level)
SpleenBcell.level.dark <- rgb2hsv(SpleenBcell.level.dark)
SpleenBcell.level.dark[3,] <- SpleenBcell.level.dark[3,]*0.5
SpleenBcell.level.dark <- hsv(SpleenBcell.level.dark[1,], SpleenBcell.level.dark[2,], SpleenBcell.level.dark[3,], 1)
SpleenBcell.level.dark[which(is.na(SpleenBcell_cor))] <- NA

SpleenMacrophage_cor <- CellsCor[2,(5:616)]

#SpleenMacrophage.level <- cRamp.orig(SpleenMacrophage_cor, colors = c("white","white","red","darkred"))
SpleenMacrophage.level <- as.factor(findInterval(SpleenMacrophage_cor, c(min(SpleenMacrophage_cor),0.05,0.1,0.2,max(SpleenMacrophage_cor)), rightmost.closed = FALSE, all.inside = FALSE))
SpleenMacrophage.level <- c("white","pink","red","darkred")[SpleenMacrophage.level] #c("white","greenyellow", "green3", "darkgreen")
SpleenMacrophage.level.dark <- col2rgb(SpleenMacrophage.level)
SpleenMacrophage.level.dark <- rgb2hsv(SpleenMacrophage.level.dark)
SpleenMacrophage.level.dark[3,] <- SpleenMacrophage.level.dark[3,]*0.5
SpleenMacrophage.level.dark <- hsv(SpleenMacrophage.level.dark[1,], SpleenMacrophage.level.dark[2,], SpleenMacrophage.level.dark[3,], 1)
SpleenMacrophage.level.dark[which(is.na(SpleenMacrophage_cor))] <- NA

SpleenTcell_cor <- CellsCor[3,(5:616)]

#SpleenTcell.level <- cRamp.orig(SpleenTcell_cor, colors = c("white","white","red","darkred"))
SpleenTcell.level <- as.factor(findInterval(SpleenTcell_cor, c(min(SpleenTcell_cor),0.05,0.1,0.2,max(SpleenTcell_cor)), rightmost.closed = FALSE, all.inside = FALSE))
SpleenTcell.level <- c("white","pink","red","darkred")[SpleenTcell.level] #c("white","greenyellow", "green3", "darkgreen")
SpleenTcell.level.dark <- col2rgb(SpleenTcell.level)
SpleenTcell.level.dark <- rgb2hsv(SpleenTcell.level.dark)
SpleenTcell.level.dark[3,] <- SpleenTcell.level.dark[3,]*0.5
SpleenTcell.level.dark <- hsv(SpleenTcell.level.dark[1,], SpleenTcell.level.dark[2,], SpleenTcell.level.dark[3,], 1)
SpleenTcell.level.dark[which(is.na(SpleenTcell_cor))] <- NA

Mitochondrial_cor <- CellsCor[4,(5:616)]

#Mitochondrial.level <- cRamp.orig(Mitochondrial_cor, colors = c("white","white","red","darkred"))
Mitochondrial.level <- as.factor(findInterval(Mitochondrial_cor, c(min(Mitochondrial_cor),0.05,0.1,0.2,max(Mitochondrial_cor)), rightmost.closed = FALSE, all.inside = FALSE))
Mitochondrial.level <- c("white","pink","red","darkred")[Mitochondrial.level] #c("white","greenyellow", "green3", "darkgreen")
Mitochondrial.level.dark <- col2rgb(Mitochondrial.level)
Mitochondrial.level.dark <- rgb2hsv(Mitochondrial.level.dark)
Mitochondrial.level.dark[3,] <- Mitochondrial.level.dark[3,]*0.5
Mitochondrial.level.dark <- hsv(Mitochondrial.level.dark[1,], Mitochondrial.level.dark[2,], Mitochondrial.level.dark[3,], 1)
Mitochondrial.level.dark[which(is.na(Mitochondrial_cor))] <- NA

SpleenCors <- as.data.frame(cbind(SpleenBcell_cor, SpleenMacrophage_cor, SpleenTcell_cor))
SpleenCors$Max<-pmax(SpleenCors$SpleenBcell_cor,SpleenCors$SpleenMacrophage_cor, SpleenCors$SpleenTcell_cor)
SpleenBcells <- rownames(SpleenCors)[which(SpleenCors$SpleenBcell_cor == SpleenCors$Max)]
SpleenMacrophages <- rownames(SpleenCors)[which(SpleenCors$SpleenMacrophage_cor == SpleenCors$Max)]
SpleenTcells <- rownames(SpleenCors)[which(SpleenCors$SpleenTcell_cor == SpleenCors$Max)]

rpkmdata <- t(rpkmdata)
CorCellType <- rep(NA, nrow(rpkmdata))
names(CorCellType) <- rownames(rpkmdata)
CorCellType[which(names(CorCellType) %in% SpleenBcells)] <- "#00A08A"
CorCellType[which(names(CorCellType) %in% SpleenMacrophages)] <- "#EBCC2A"
CorCellType[which(names(CorCellType) %in% SpleenTcells)] <- "#F21A00"

CorCellType.dark <- col2rgb(CorCellType)
CorCellType.dark <- rgb2hsv(CorCellType.dark)
CorCellType.dark[3,] <- CorCellType.dark[3,]*0.5
CorCellType.dark <- hsv(CorCellType.dark[1,], CorCellType.dark[2,], CorCellType.dark[3,], 1)
CorCellType.dark[which(is.na(CorCellType))] <- NA

load("../WD/QCCells_Features.RData")
S_HEK_30m_1xFeat <- CellFeat[grep('^S_HEK_30m_1x', rownames(CellFeat)),]
S_CAP_30m_1xFeat <- CellFeat[grep('^S_CAP_30m_1x', rownames(CellFeat)),]
CellFeat <- rbind(S_HEK_30m_1xFeat, S_CAP_30m_1xFeat)

RandomForest <- rep(NA, nrow(rpkmdata))
names(RandomForest) <- rownames(rpkmdata)
RandomForest[which(names(RandomForest) %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#78B7C5")])] <- "#00A08A"
RandomForest[which(names(RandomForest) %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#E1AF00")])] <- "#EBCC2A"
RandomForest[which(names(RandomForest) %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#7294D4")])] <- "#F21A00"

RandomForest.dark <- col2rgb(RandomForest)
RandomForest.dark <- rgb2hsv(RandomForest.dark)
RandomForest.dark[3,] <- RandomForest.dark[3,]*0.5
RandomForest.dark <- hsv(RandomForest.dark[1,], RandomForest.dark[2,], RandomForest.dark[3,], 1)
RandomForest.dark[which(is.na(RandomForest))] <- NA

HEK <- length(which(CellFeat$CellType == "#CCCCCC"))
CAP <- length(which(CellFeat$CellType == "#666666"))
Bcells <- length(SpleenBcells)
aiBcells <- length(which(CellFeat$RandomForest == "#78B7C5"))
BcellsHEK <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% SpleenBcells))
BcellsCAP <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% SpleenBcells))
BcellsHEKai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#78B7C5")]))
BcellsCAPai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#78B7C5")]))
Macrophages <- length(SpleenMacrophages)
aiMacrophages <- length(which(CellFeat$RandomForest == "#E1AF00"))
MacrophagesHEK <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% SpleenMacrophages))
MacrophagesCAP <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% SpleenMacrophages))
MacrophagesHEKai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#E1AF00")]))
MacrophagesCAPai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#E1AF00")]))
Tcells <- length(SpleenTcells)
aiTcells <- length(which(CellFeat$RandomForest == "#7294D4"))
TcellsHEK <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% SpleenTcells))
TcellsCAP <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% SpleenTcells))
TcellsHEKai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#7294D4")]))
TcellsCAPai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#7294D4")]))
Numbers <- as.data.frame(c(HEK,CAP, Bcells, BcellsHEK, BcellsCAP, Macrophages, MacrophagesHEK, MacrophagesCAP, Tcells, TcellsHEK, TcellsCAP))
aiNumbers <- as.data.frame(c(HEK,CAP, aiBcells, BcellsHEKai, BcellsCAPai, aiMacrophages, MacrophagesHEKai, MacrophagesCAPai, aiTcells, TcellsHEKai, TcellsCAPai))
rownames(Numbers) <- c("HEK","CAP", "Bcells", "BcellsHEK", "BcellsCAP", "Macrophages", "MacrophagesHEK", "MacrophagesCAP", "Tcells", "TcellsHEK", "TcellsCAP")
rownames(aiNumbers) <- c("HEK","CAP", "Bcells", "BcellsHEK", "BcellsCAP", "Macrophages", "MacrophagesHEK", "MacrophagesCAP", "Tcells", "TcellsHEK", "TcellsCAP")
write.csv(Numbers, file="SpleenNumbers.csv")
write.csv(aiNumbers, file="aiSpleenNumbers.csv")

###################################################################

load("../WD/CD63spleen1x30mCellsQC_tSNE-100noid-3nod-200non_299100TabulaGenes_tSNEresults.Rdata")
Cluster <- rep(NA, nrow(rpkmdata))
names(Cluster) <- rownames(rpkmdata)
Cluster1 <- names(which(tSNEresults$g.info[[1]] == "grey"))
Cluster2 <- names(which(tSNEresults$g.info[[1]] == "orange"))
Cluster3 <- names(which(tSNEresults$g.info[[1]] == "yellow"))
#Cluster4 <- names(which(tSNEresults$g.info[[1]] == "red"))
Cluster[which(names(Cluster) %in% Cluster1)] <- "#00A08A"
Cluster[which(names(Cluster) %in% Cluster2)] <- "#EBCC2A"
Cluster[which(names(Cluster) %in% Cluster3)] <- "#F21A00"
#Cluster[which(names(Cluster) %in% Cluster4)] <- "#7294D4"

Cluster.dark <- col2rgb(Cluster)
Cluster.dark <- rgb2hsv(Cluster.dark)
Cluster.dark[3,] <- Cluster.dark[3,]*0.5
Cluster.dark <- hsv(Cluster.dark[1,], Cluster.dark[2,], Cluster.dark[3,], 1)
Cluster.dark[which(is.na(Cluster))] <- NA

### Plot data
pdf("Spleen30mMetaCells.pdf", width=10, height=10)
par(mfrow=c(4,4), mar=c(1,1,1,1))

load("../WD/CD63spleen1x30mCellsQC_tSNE-100noid-3nod-20non_299100TabulaGenes_tSNEresults.Rdata")
g <- tSNEresults$g
g.info <- tSNEresults$g.info
g.layout <- tSNEresults$g.layout

#adjacency.plot2(g, main = sprintf("Organ %s", ncol(data)), CellFeat$Organ, CellFeat$Organ.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("CellType %s", ncol(data)), CellFeat$CellType, CellFeat$CellType.dark, g.layout,6)
#adjacency.plot2(g, main = sprintf("TimePoint %s", ncol(data)), CellFeat$TimePoint, CellFeat$TimePoint.dark, g.layout,6)
#adjacency.plot2(g, main = sprintf("Dose %s", ncol(data)), CellFeat$Dose, CellFeat$Dose.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("DetGenes %s", ncol(data)), CellFeat$DetGenes, CellFeat$DetGenes.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("ReadCount %s", ncol(data)), CellFeat$ReadCount, CellFeat$ReadCount.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("clusters %s", ncol(data)), CellFeat$clusters, CellFeat$clusters.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("Infomap %s", ncol(data)), g.info[[1]], g.info[[2]], g.layout,6)
adjacency.plot2(g, main = sprintf("Cluster %s", ncol(data)), Cluster, Cluster.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("CorCellType %s", ncol(data)), CorCellType, CorCellType.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("RandomForest %s", ncol(data)), CellFeat$RandomForest, CellFeat$RandomForest.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("RandomForest %s", ncol(data)), RandomForest, RandomForest.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("Batch %s", ncol(data)), CellFeat$Batch, CellFeat$Batch.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("SpleenBcellMetaCell %s", ncol(data)), SpleenBcell.level, SpleenBcell.level.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("SpleenMacrophageMetaCell %s", ncol(data)), SpleenMacrophage.level, SpleenMacrophage.level.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("SpleenTcellMetaCell %s", ncol(data)), SpleenTcell.level, SpleenTcell.level.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("Mitochondrial %s", ncol(data)), Mitochondrial.level, Mitochondrial.level.dark, g.layout,6)

### Plot expression data

rpkmdata<-read.table("../WD/S_30m_1x_2x_all_batch_correction_norm_express.txt", header=T)

S_HEK_30m_1x <- rpkmdata[,grep('^S_HEK_30m_1x', colnames(rpkmdata))]
#S_HEK_30m_2x <- rpkmdata[,grep('^S_HEK_30m_2x', colnames(rpkmdata))]
#S_HEK_3h_2x <- rpkmdata[,grep('^S_HEK_3h_2x', colnames(rpkmdata))]
S_CAP_30m_1x <- rpkmdata[,grep('^S_CAP_30m_1x', colnames(rpkmdata))]
#S_CAP_30m_2x <- rpkmdata[,grep('^S_CAP_30m_2x', colnames(rpkmdata))]
#S_CAP_3h_2x <- rpkmdata[,grep('^S_CAP_3h_2x', colnames(rpkmdata))]

rpkmdata <- cbind(S_HEK_30m_1x, S_CAP_30m_1x)

data <- t(as.data.frame(data))
data[data==0] <- 0.0001
data <- log2(data)

w<-c("ENSMUSG00000069516", "ENSMUSG00000030724","ENSMUSG00000044206","ENSMUSG00000020884")#LYZ2, CD19, VSIG4, ASGR1

tempdata <- t(rpkmdata)
tempdata[tempdata==0] <- 0.0001
tempdata <- log2(tempdata)
tempdata <- tempdata[,w[which(w %in% colnames(tempdata))]]


for (i in colnames(tempdata)){
  exp.level <- as.factor(findInterval(tempdata[,i], c(min(rpkmdata),1,4,6.65,9,11,max(rpkmdata)), rightmost.closed = FALSE, all.inside = FALSE))
  exp.level <- c("white", "darkseagreen1", "green3", "green4", "darkgreen", "black")[exp.level]
  exp.level <- col2rgb(exp.level)
  exp.level <- rgb2hsv(exp.level)
  exp.level <- hsv(exp.level[1,], exp.level[2,], exp.level[3,], 1)
  exp.level[which(exp.level == "#FFFFFFFF")] <- rep(makeTransparent("white", 0.5), sum(exp.level == "#FFFFFFFF"))
  
  exp.level.dark <- col2rgb(exp.level)
  exp.level.dark <- rgb2hsv(exp.level.dark)
  exp.level.dark[3,] <- exp.level.dark[3,]*0.5
  exp.level.dark <- hsv(exp.level.dark[1,], exp.level.dark[2,], exp.level.dark[3,], 1)
  exp.level.dark[which(is.na(exp.level))] <- NA
  
  adjacency.plot2(g, main = i, exp.level, exp.level.dark, g.layout,6)
}

dev.off()

#################################

#SpleenClusters <- as.data.frame(cbind(Cluster, SpleenBcell_cor, SpleenMacrophage_cor, SpleenTcell_cor))
#pleenCluster1 <- SpleenClusters[which(SpleenClusters$Cluster == "#00A08A"),]
#SpleenCluster2 <- SpleenClusters[which(SpleenClusters$Cluster == "#EBCC2A"),]
#SpleenCluster3 <- SpleenClusters[which(SpleenClusters$Cluster == "#F21A00"),]
#SpleenCluster4 <- SpleenClusters[which(SpleenClusters$Cluster == "#7294D4"),]
#write.csv(SpleenCluster1, "SpleenClusterGreen.csv")
#write.csv(SpleenCluster2, "SpleenClusterYellow.csv")
#write.csv(SpleenCluster3, "SpleenClusterRed.csv")
#write.csv(SpleenCluster4, "SpleenClusterPurple.csv")
#write.csv(SpleenClusters, "SpleenClusters.csv")


