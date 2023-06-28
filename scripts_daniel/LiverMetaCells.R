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

rpkmdata<-read.table("../WD/LI_30m_1x_2x_all_2_batch_correction_norm_express.txt", header=T)
data <- t(rpkmdata)

LiverGenes <- read.csv("../WD/LiverBcell100GenesEnsembl.txt")
Liver <- rep(0, ncol(data))
names(Liver) <- colnames(data)
Liver[which(names(Liver) %in% LiverGenes$ensembl_gene_id)] <- 100
LiverBcell <- Liver

LiverGenes <- read.csv("../WD/LiverHepatocyte100GenesEnsembl.txt")
Liver <- rep(0, ncol(data))
names(Liver) <- colnames(data)
Liver[which(names(Liver) %in% LiverGenes$ensembl_gene_id)] <- 100
LiverHepatocyte <- Liver

LiverGenes <- read.csv("../WD/LiverKupffer100GenesEnsembl.txt")
Liver <- rep(0, ncol(data))
names(Liver) <- colnames(data)
Liver[which(names(Liver) %in% LiverGenes$ensembl_gene_id)] <- 100
LiverKupffer <- Liver

LiverGenes <- read.csv("../WD/LiverNKcell100GenesEnsembl.txt")
Liver <- rep(0, ncol(data))
names(Liver) <- colnames(data)
Liver[which(names(Liver) %in% LiverGenes$ensembl_gene_id)] <- 100
LiverNKcell <- Liver

LiverGenes <- read.csv("../WD/LiverEndothelial100GenesEnsembl.txt")
Liver <- rep(0, ncol(data))
names(Liver) <- colnames(data)
Liver[which(names(Liver) %in% LiverGenes$ensembl_gene_id)] <- 100
LiverEndothelial <- Liver

MitoGenes <- read.table("../WD/mouse_mitochondrial_genes_list_annotation.txt", sep = "\t", header = T)
Liver <- rep(0, ncol(data))
names(Liver) <- colnames(data)
Liver[which(names(Liver) %in% MitoGenes$ensembl_gene_id)] <- 100
Mitochondrial <- Liver

MetaLiver <- as.data.frame(matrix(NA, 5, ncol(data)))
MetaLiver[1,] <- as.numeric(LiverBcell)
MetaLiver[2,] <- as.numeric(LiverEndothelial)
MetaLiver[3,] <- as.numeric(LiverHepatocyte)
MetaLiver[4,] <- as.numeric(LiverKupffer)
MetaLiver[5,] <- as.numeric(LiverNKcell)
MetaLiver[6,] <- as.numeric(Mitochondrial)
colnames(MetaLiver) <- colnames(data)
rownames(MetaLiver) <- c("Bcell","Endothelial","Hepatocyte","Kupffer","NKcell","Mitochondrial")
MetaLiver <- t(MetaLiver)

################################################

rpkmdata<-read.table("../WD/LI_30m_1x_2x_all_2_batch_correction_norm_express.txt", header=T)

LI_HEK_30m_1x <- rpkmdata[,grep('^LI_HEK_30m_1x', colnames(rpkmdata))]
#LI_HEK_30m_2x <- rpkmdata[,grep('^LI_HEK_30m_2x', colnames(rpkmdata))]
#LI_HEK_3h_2x <- rpkmdata[,grep('^LI_HEK_3h_2x', colnames(rpkmdata))]
LI_CAP_30m_1x <- rpkmdata[,grep('^LI_CAP_30m_1x', colnames(rpkmdata))]
#LI_CAP_30m_2x <- rpkmdata[,grep('^LI_CAP_30m_2x', colnames(rpkmdata))]
#LI_CAP_3h_2x <- rpkmdata[,grep('^LI_CAP_3h_2x', colnames(rpkmdata))]

rpkmdata <- cbind(LI_HEK_30m_1x, LI_CAP_30m_1x)
data <- cbind(MetaLiver, LI_HEK_30m_1x, LI_CAP_30m_1x)

data <- t(as.data.frame(data))
data[data==0] <- 0.0001
data <- log2(data)

#var.rpkms<-var.genes(data)
#data <- subset(data, select = colnames(data) %in% var.rpkms)

genes <- read.csv("../WD/Liver100TabulaGenesEnsembl.txt", header=F)
#genes <- c(genes$ensembl_gene_id, MitoGenes$ensembl_gene_id)
data <- subset(data, select = colnames(data) %in% genes$V1)

CellsCor=(cor(t(data)))
matrix<-as.matrix(CellsCor)

LiverBcell_cor <- CellsCor[1,(7:655)]

#LiverBcell.level <- cRamp.orig(LiverBcell_cor, colors = c("white","white","red","darkred"))
LiverBcell.level <- as.factor(findInterval(LiverBcell_cor, c(min(LiverBcell_cor),0.05,0.1,0.2,max(LiverBcell_cor)), rightmost.closed = FALSE, all.inside = FALSE))
LiverBcell.level <- c("white","pink","red","darkred")[LiverBcell.level]
LiverBcell.level.dark <- col2rgb(LiverBcell.level)
LiverBcell.level.dark <- rgb2hsv(LiverBcell.level.dark)
LiverBcell.level.dark[3,] <- LiverBcell.level.dark[3,]*0.5
LiverBcell.level.dark <- hsv(LiverBcell.level.dark[1,], LiverBcell.level.dark[2,], LiverBcell.level.dark[3,], 1)
LiverBcell.level.dark[which(is.na(LiverBcell_cor))] <- NA

LiverEndothelial_cor <- CellsCor[2,(7:655)]

#LiverEndothelial.level <- cRamp.orig(LiverEndothelial_cor, colors = c("white","white","red","darkred"))
LiverEndothelial.level <- as.factor(findInterval(LiverEndothelial_cor, c(min(LiverEndothelial_cor),0.05,max(LiverEndothelial_cor)), rightmost.closed = FALSE, all.inside = FALSE))
LiverEndothelial.level <- c("white","pink")[LiverEndothelial.level]
LiverEndothelial.level.dark <- col2rgb(LiverEndothelial.level)
LiverEndothelial.level.dark <- rgb2hsv(LiverEndothelial.level.dark)
LiverEndothelial.level.dark[3,] <- LiverEndothelial.level.dark[3,]*0.5
LiverEndothelial.level.dark <- hsv(LiverEndothelial.level.dark[1,], LiverEndothelial.level.dark[2,], LiverEndothelial.level.dark[3,], 1)
LiverEndothelial.level.dark[which(is.na(LiverEndothelial_cor))] <- NA

LiverHepatocyte_cor <- (CellsCor[3,(7:655)]+.1)*4

#LiverHepatocyte.level <- cRamp.orig(LiverHepatocyte_cor, colors = c("white","white","red","darkred"))
LiverHepatocyte.level <- as.factor(findInterval(LiverHepatocyte_cor, c(min(LiverHepatocyte_cor),0.05,0.1,0.2,max(LiverHepatocyte_cor)), rightmost.closed = FALSE, all.inside = FALSE))
LiverHepatocyte.level <- c("white","pink","red","darkred")[LiverHepatocyte.level]
LiverHepatocyte.level.dark <- col2rgb(LiverHepatocyte.level)
LiverHepatocyte.level.dark <- rgb2hsv(LiverHepatocyte.level.dark)
LiverHepatocyte.level.dark[3,] <- LiverHepatocyte.level.dark[3,]*0.5
LiverHepatocyte.level.dark <- hsv(LiverHepatocyte.level.dark[1,], LiverHepatocyte.level.dark[2,], LiverHepatocyte.level.dark[3,], 1)
LiverHepatocyte.level.dark[which(is.na(LiverHepatocyte_cor))] <- NA

LiverKupffer_cor <- CellsCor[4,(7:655)]

#LiverKupffer.level <- cRamp.orig(LiverKupffer_cor, colors = c("white","white","red","darkred"))
LiverKupffer.level <- as.factor(findInterval(LiverKupffer_cor, c(min(LiverKupffer_cor),0.05,0.1,0.2,max(LiverKupffer_cor)), rightmost.closed = FALSE, all.inside = FALSE))
LiverKupffer.level <- c("white","pink","red","darkred")[LiverKupffer.level]
LiverKupffer.level.dark <- col2rgb(LiverKupffer.level)
LiverKupffer.level.dark <- rgb2hsv(LiverKupffer.level.dark)
LiverKupffer.level.dark[3,] <- LiverKupffer.level.dark[3,]*0.5
LiverKupffer.level.dark <- hsv(LiverKupffer.level.dark[1,], LiverKupffer.level.dark[2,], LiverKupffer.level.dark[3,], 1)
LiverKupffer.level.dark[which(is.na(LiverKupffer_cor))] <- NA

LiverNKcell_cor <- CellsCor[5,(7:655)]

#LiverNKcell.level <- cRamp.orig(LiverNKcell_cor, colors = c("white","white","red","darkred"))
LiverNKcell.level <- as.factor(findInterval(LiverNKcell_cor, c(min(LiverNKcell_cor),0.05,0.1,0.2,max(LiverNKcell_cor)), rightmost.closed = FALSE, all.inside = FALSE))
LiverNKcell.level <- c("white","pink","red","darkred")[LiverNKcell.level]
LiverNKcell.level.dark <- col2rgb(LiverNKcell.level)
LiverNKcell.level.dark <- rgb2hsv(LiverNKcell.level.dark)
LiverNKcell.level.dark[3,] <- LiverNKcell.level.dark[3,]*0.5
LiverNKcell.level.dark <- hsv(LiverNKcell.level.dark[1,], LiverNKcell.level.dark[2,], LiverNKcell.level.dark[3,], 1)
LiverNKcell.level.dark[which(is.na(LiverNKcell_cor))] <- NA

Mitochondrial_cor <- CellsCor[6,(7:655)]

Mitochondrial.level <- cRamp.orig(Mitochondrial_cor, colors = c("white","white","red","darkred"))
Mitochondrial.level <- as.factor(findInterval(Mitochondrial_cor, c(min(Mitochondrial_cor),0.05,0.1,0.2,max(Mitochondrial_cor)), rightmost.closed = FALSE, all.inside = FALSE))
Mitochondrial.level <- c("white", "pink","red","darkred")[Mitochondrial.level]
Mitochondrial.level.dark <- col2rgb(Mitochondrial.level)
Mitochondrial.level.dark <- rgb2hsv(Mitochondrial.level.dark)
Mitochondrial.level.dark[3,] <- Mitochondrial.level.dark[3,]*0.5
Mitochondrial.level.dark <- hsv(Mitochondrial.level.dark[1,], Mitochondrial.level.dark[2,], Mitochondrial.level.dark[3,], 1)
Mitochondrial.level.dark[which(is.na(Mitochondrial_cor))] <- NA

LiverCors <- as.data.frame(cbind(LiverBcell_cor, LiverEndothelial_cor, LiverHepatocyte_cor, LiverKupffer_cor, LiverNKcell_cor))
LiverCors$Max<-pmax(LiverCors$LiverBcell_cor,LiverCors$LiverEndothelial_cor, LiverCors$LiverHepatocyte_cor, LiverCors$LiverKupffer_cor, LiverCors$LiverNKcell_cor)
LiverBcells <- rownames(LiverCors)[which(LiverCors$LiverBcell_cor == LiverCors$Max)]
LiverEndothelials <- rownames(LiverCors)[which(LiverCors$LiverEndothelial_cor == LiverCors$Max)]
LiverHepatocytes <- rownames(LiverCors)[which(LiverCors$LiverHepatocyte_cor == LiverCors$Max)]
LiverKupffers <- rownames(LiverCors)[which(LiverCors$LiverKupffer_cor == LiverCors$Max)]
LiverNKcells <- rownames(LiverCors)[which(LiverCors$LiverNKcell_cor == LiverCors$Max)]

rpkmdata <- t(rpkmdata)
CorCellType <- rep(NA, nrow(rpkmdata))
names(CorCellType) <- rownames(rpkmdata)
CorCellType[which(names(CorCellType) %in% LiverBcells)] <- "#00A08A"
CorCellType[which(names(CorCellType) %in% LiverEndothelials)] <- "#7294D4"
CorCellType[which(names(CorCellType) %in% LiverHepatocytes)] <- "#F21A00"
CorCellType[which(names(CorCellType) %in% LiverKupffers)] <- "#EBCC2A"
CorCellType[which(names(CorCellType) %in% LiverNKcells)] <- "#3B9AB2"

CorCellType.dark <- col2rgb(CorCellType)
CorCellType.dark <- rgb2hsv(CorCellType.dark)
CorCellType.dark[3,] <- CorCellType.dark[3,]*0.5
CorCellType.dark <- hsv(CorCellType.dark[1,], CorCellType.dark[2,], CorCellType.dark[3,], 1)
CorCellType.dark[which(is.na(CorCellType))] <- NA

load("../WD/QCCells_Features.RData")
LI_HEK_30m_1xFeat <- CellFeat[grep('^LI_HEK_30m_1x', rownames(CellFeat)),]
LI_CAP_30m_1xFeat <- CellFeat[grep('^LI_CAP_30m_1x', rownames(CellFeat)),]
CellFeat <- rbind(LI_HEK_30m_1xFeat, LI_CAP_30m_1xFeat)

RandomForest <- rep(NA, nrow(rpkmdata))
names(RandomForest) <- rownames(rpkmdata)
RandomForest[which(names(RandomForest) %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#78B7C5")])] <- "#00A08A"
RandomForest[which(names(RandomForest) %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#00A08A")])] <- "#F21A00"
RandomForest[which(names(RandomForest) %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#E1AF00")])] <- "#EBCC2A"
RandomForest[which(names(RandomForest) %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#D8A499")])] <- "#3B9AB2"

RandomForest.dark <- col2rgb(RandomForest)
RandomForest.dark <- rgb2hsv(RandomForest.dark)
RandomForest.dark[3,] <- RandomForest.dark[3,]*0.5
RandomForest.dark <- hsv(RandomForest.dark[1,], RandomForest.dark[2,], RandomForest.dark[3,], 1)
RandomForest.dark[which(is.na(RandomForest))] <- NA

HEK <- length(which(CellFeat$CellType == "#CCCCCC"))
CAP <- length(which(CellFeat$CellType == "#666666"))
Bcells <- length(LiverBcells)
aiBcells <- length(which(CellFeat$RandomForest == "#78B7C5"))
BcellsHEK <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% LiverBcells))
BcellsCAP <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% LiverBcells))
BcellsHEKai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#78B7C5")]))
BcellsCAPai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#78B7C5")]))
Endothelials <- length(LiverEndothelials)
aiEndothelials <- length(which(CellFeat$RandomForest == "#"))
EndothelialsHEK <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% LiverEndothelials))
EndothelialsCAP <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% LiverEndothelials))
EndothelialsHEKai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#")]))
EndothelialsCAPai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#")]))
Hepatocytes <- length(LiverHepatocytes)
aiHepatocytes <- length(which(CellFeat$RandomForest == "#00A08A"))
HepatocytesHEK <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% LiverHepatocytes))
HepatocytesCAP <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% LiverHepatocytes))
HepatocytesHEKai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#00A08A")]))
HepatocytesCAPai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#00A08A")]))
Kupffers <- length(LiverKupffers)
aiKupffers <- length(which(CellFeat$RandomForest == "#E1AF00"))
KupffersHEK <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% LiverKupffers))
KupffersCAP <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% LiverKupffers))
KupffersHEKai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#E1AF00")]))
KupffersCAPai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#E1AF00")]))
NKcells <- length(LiverNKcells)
aiNKcells <- length(which(CellFeat$RandomForest == "#D8A499"))
NKcellsHEK <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% LiverNKcells))
NKcellsCAP <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% LiverNKcells))
NKcellsHEKai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#D8A499")]))
NKcellsCAPai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#D8A499")]))
Numbers <- as.data.frame(c(HEK,CAP, Bcells, BcellsHEK, BcellsCAP, Endothelials, EndothelialsHEK, EndothelialsCAP, Hepatocytes, HepatocytesHEK, HepatocytesCAP,
                           Kupffers, KupffersHEK, KupffersCAP, NKcells, NKcellsHEK, NKcellsCAP))
aiNumbers <- as.data.frame(c(HEK,CAP, aiBcells, BcellsHEKai, BcellsCAPai, aiEndothelials, EndothelialsHEKai, EndothelialsCAPai, aiHepatocytes, HepatocytesHEKai, HepatocytesCAPai,
                             aiKupffers, KupffersHEKai, KupffersCAPai, aiNKcells, NKcellsHEKai, NKcellsCAPai))
rownames(Numbers) <- c("HEK","CAP", "Bcells", "BcellsHEK", "BcellsCAP", "Endothelials", "EndothelialsHEK", "EndothelialsCAP", "Hepatocytes", "HepatocytesHEK", "HepatocytesCAP", 
                       "Kupffers", "KupffersHEK", "KupffersCAP", "NKcells", "NKcellsHEK", "NKcellsCAP")
rownames(aiNumbers) <- c("HEK","CAP", "Bcells", "BcellsHEK", "BcellsCAP", "Endothelials", "EndothelialsHEK", "EndothelialsCAP", "Hepatocytes", "HepatocytesHEK", "HepatocytesCAP", 
                         "Kupffers", "KupffersHEK", "KupffersCAP", "NKcells", "NKcellsHEK", "NKcellsCAP")
write.csv(Numbers, file="LiverNumbers.csv")
write.csv(aiNumbers, file="aiLiverNumbers.csv")

###################################################################

load("../WD/CD63liver1x30mCellsQC_tSNE-100noid-3nod-200non_497Tabula100Genes_tSNEresults.Rdata")
Cluster <- rep(NA, nrow(rpkmdata))
names(Cluster) <- rownames(rpkmdata)
Cluster1 <- names(which(tSNEresults$g.info[[1]] == "grey"))
Cluster2 <- names(which(tSNEresults$g.info[[1]] == "orange"))
Cluster3 <- names(which(tSNEresults$g.info[[1]] == "yellow"))

Cluster[which(names(Cluster) %in% Cluster1)] <- "#3B9AB2"
Cluster[which(names(Cluster) %in% Cluster2)] <- "#EBCC2A"
Cluster[which(names(Cluster) %in% Cluster3)] <- "#F21A00"

Cluster.dark <- col2rgb(Cluster)
Cluster.dark <- rgb2hsv(Cluster.dark)
Cluster.dark[3,] <- Cluster.dark[3,]*0.5
Cluster.dark <- hsv(Cluster.dark[1,], Cluster.dark[2,], Cluster.dark[3,], 1)
Cluster.dark[which(is.na(Cluster))] <- NA

### Plot data
pdf("Liver30mMetaCells.pdf", width=10, height=10)
par(mfrow=c(4,4), mar=c(1,1,1,1))

load("../WD/CD63liver1x30mCellsQC_tSNE-100noid-3nod-20non_497Tabula100Genes_tSNEresults.Rdata")
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
adjacency.plot2(g, main = sprintf("CorCellType %s", ncol(data)), CorCellType, CorCellType.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("Cluster %s", ncol(data)), Cluster, Cluster.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("RandomForest %s", ncol(data)), CellFeat$RandomForest, CellFeat$RandomForest.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("RandomForest %s", ncol(data)), RandomForest, RandomForest.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("Batch %s", ncol(data)), CellFeat$Batch, CellFeat$Batch.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("LiverBcellMetaCell %s", ncol(data)), LiverBcell.level, LiverBcell.level.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("LiverEndothelialMetaCell %s", ncol(data)), LiverEndothelial.level, LiverEndothelial.level.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("LiverHepatocyteMetaCell %s", ncol(data)), LiverHepatocyte.level, LiverHepatocyte.level.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("LiverKupfferMetaCell %s", ncol(data)), LiverKupffer.level, LiverKupffer.level.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("LiverNKcellMetaCell %s", ncol(data)), LiverNKcell.level, LiverNKcell.level.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("Mitochondrial %s", ncol(data)), Mitochondrial.level, Mitochondrial.level.dark, g.layout,6)

### Plot expression data
rpkmdata<-read.table("../WD/LI_30m_1x_2x_all_2_batch_correction_norm_express.txt", header=T)

LI_HEK_30m_1x <- rpkmdata[,grep('^LI_HEK_30m_1x', colnames(rpkmdata))]
#LI_HEK_30m_2x <- rpkmdata[,grep('^LI_HEK_30m_2x', colnames(rpkmdata))]
#LI_HEK_3h_2x <- rpkmdata[,grep('^LI_HEK_3h_2x', colnames(rpkmdata))]
LI_CAP_30m_1x <- rpkmdata[,grep('^LI_CAP_30m_1x', colnames(rpkmdata))]
#LI_CAP_30m_2x <- rpkmdata[,grep('^LI_CAP_30m_2x', colnames(rpkmdata))]
#LI_CAP_3h_2x <- rpkmdata[,grep('^LI_CAP_3h_2x', colnames(rpkmdata))]

rpkmdata <- cbind(LI_HEK_30m_1x, LI_CAP_30m_1x)

data <- t(rpkmdata)
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

#LiverClusters <- as.data.frame(cbind(Cluster, LiverBcell_cor, LiverEndothelial_cor, LiverHepatocyte_cor, LiverKupffer_cor, LiverNKcell_cor))
#LiverCluster1 <- LiverClusters[which(LiverClusters$Cluster == "#3B9AB2"),]
##LiverCluster2 <- LiverClusters[which(LiverClusters$Cluster == "#EBCC2A"),]
#LiverCluster3 <- LiverClusters[which(LiverClusters$Cluster == "#F21A00"),]
#write.csv(LiverCluster1, "LiverClusterBlue.csv")
#write.csv(LiverCluster2, "LiverClusterYellow.csv")
#write.csv(LiverCluster3, "LiverClusterRed.csv")
