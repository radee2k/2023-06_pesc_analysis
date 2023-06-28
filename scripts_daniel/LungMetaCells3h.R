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

rpkmdata<-read.table("../WD/Lu_3h_2x_all_batch_correction_norm_express.txt", header=T)
data <- t(rpkmdata)

LungGenes <- read.csv("../WD/LungBcell100GenesEnsembl.txt")
Lung <- rep(0, ncol(data))
names(Lung) <- colnames(data)
Lung[which(names(Lung) %in% LungGenes$ensembl_gene_id)] <- 100
LungBcell <- Lung

LungGenes <- read.csv("../WD/LungCiliated100GenesEnsembl.txt")
Lung <- rep(0, ncol(data))
names(Lung) <- colnames(data)
Lung[which(names(Lung) %in% LungGenes$ensembl_gene_id)] <- 100
LungCiliated <- Lung

LungGenes <- read.csv("../WD/LungClasMono100GenesEnsembl.txt")
Lung <- rep(0, ncol(data))
names(Lung) <- colnames(data)
Lung[which(names(Lung) %in% LungGenes$ensembl_gene_id)] <- 100
LungClassicMonocyte <- Lung

LungGenes <- read.csv("../WD/LungEndothelial100GenesEnsembl.txt")
Lung <- rep(0, ncol(data))
names(Lung) <- colnames(data)
Lung[which(names(Lung) %in% LungGenes$ensembl_gene_id)] <- 100
LungEndothelial <- Lung

LungGenes <- read.csv("../WD/LungEpithelial100GenesEnsembl.txt")
Lung <- rep(0, ncol(data))
names(Lung) <- colnames(data)
Lung[which(names(Lung) %in% LungGenes$ensembl_gene_id)] <- 100
LungEpithelial <- Lung

LungGenes <- read.csv("../WD/LungLeukocyte100GenesEnsembl.txt")
Lung <- rep(0, ncol(data))
names(Lung) <- colnames(data)
Lung[which(names(Lung) %in% LungGenes$ensembl_gene_id)] <- 100
LungLeukocyte <- Lung

LungGenes <- read.csv("../WD/LungMonocyte100GenesEnsembl.txt")
Lung <- rep(0, ncol(data))
names(Lung) <- colnames(data)
Lung[which(names(Lung) %in% LungGenes$ensembl_gene_id)] <- 100
LungMonocyte <- Lung

LungGenes <- read.csv("../WD/LungMyeloid100GenesEnsembl.txt")
Lung <- rep(0, ncol(data))
names(Lung) <- colnames(data)
Lung[which(names(Lung) %in% LungGenes$ensembl_gene_id)] <- 100
LungMyeloid <- Lung

LungGenes <- read.csv("../WD/LungNKcell100GenesEnsembl.txt")
Lung <- rep(0, ncol(data))
names(Lung) <- colnames(data)
Lung[which(names(Lung) %in% LungGenes$ensembl_gene_id)] <- 100
LungNKcell <- Lung

LungGenes <- read.csv("../WD/LungStromal100GenesEnsembl.txt")
Lung <- rep(0, ncol(data))
names(Lung) <- colnames(data)
Lung[which(names(Lung) %in% LungGenes$ensembl_gene_id)] <- 100
LungStromal <- Lung

LungGenes <- read.csv("../WD/LungTcell100GenesEnsembl.txt")
Lung <- rep(0, ncol(data))
names(Lung) <- colnames(data)
Lung[which(names(Lung) %in% LungGenes$ensembl_gene_id)] <- 100
LungTcell <- Lung

MitoGenes <- read.table("../WD/mouse_mitochondrial_genes_list_annotation.txt", sep = "\t", header = T)
Lung <- rep(0, ncol(data))
names(Lung) <- colnames(data)
Lung[which(names(Lung) %in% MitoGenes$ensembl_gene_id)] <- 100
Mitochondrial <- Lung

MetaLung <- as.data.frame(matrix(NA, 12, ncol(data)))
MetaLung[1,] <- as.numeric(LungBcell)
MetaLung[2,] <- as.numeric(LungCiliated)
MetaLung[3,] <- as.numeric(LungClassicMonocyte)
MetaLung[4,] <- as.numeric(LungEndothelial)
MetaLung[5,] <- as.numeric(LungEpithelial)
MetaLung[6,] <- as.numeric(LungLeukocyte)
MetaLung[7,] <- as.numeric(LungMonocyte)
MetaLung[8,] <- as.numeric(LungMyeloid)
MetaLung[9,] <- as.numeric(LungNKcell)
MetaLung[10,] <- as.numeric(LungStromal)
MetaLung[11,] <- as.numeric(LungTcell)
MetaLung[12,] <- as.numeric(Mitochondrial)
colnames(MetaLung) <- colnames(data)
rownames(MetaLung) <- c("Bcell","Ciliated","ClassicMonocyte","Endothelial","Epithelial","Leukocyte","Monocyte","Myeloid","NKcell","Stromal","Tcell","Mitochondrial")
MetaLung <- t(MetaLung)

################################################

rpkmdata<-read.table("../WD/Lu_3h_2x_all_batch_correction_norm_express.txt", header=T)

#LU_HEK_30m_1x <- rpkmdata[,grep('^LU_HEK_30m_1x', colnames(rpkmdata))]
#LI_HEK_30m_2x <- rpkmdata[,grep('^LI_HEK_30m_2x', colnames(rpkmdata))]
LU_HEK_3h_2x <- rpkmdata[,grep('^LU_HEK_3h_2x', colnames(rpkmdata))]
#LU_CAP_30m_1x <- rpkmdata[,grep('^LU_CAP_30m_1x', colnames(rpkmdata))]
#LI_CAP_30m_2x <- rpkmdata[,grep('^LI_CAP_30m_2x', colnames(rpkmdata))]
LU_CAP_3h_2x <- rpkmdata[,grep('^LU_CAP_3h_2x', colnames(rpkmdata))]

rpkmdata <- cbind(LU_HEK_3h_2x, LU_CAP_3h_2x)
data <- cbind(MetaLung, LU_HEK_3h_2x, LU_CAP_3h_2x)

data <- t(as.data.frame(data))
data[data==0] <- 0.0001
data <- log2(data)

#var.rpkms<-var.genes(data)
#data <- subset(data, select = colnames(data) %in% var.rpkms)

genes <- read.csv("../WD/Lung100TabulaGenesEnsembl.txt", header=F)
#genes <- c(genes$ensembl_gene_id, MitoGenes$ensembl_gene_id)
data <- subset(data, select = colnames(data) %in% genes$V1)

CellsCor=(cor(t(data)))
matrix<-as.matrix(CellsCor)

LungBcell_cor <- CellsCor[1,(13:446)]

#LungBcell.level <- cRamp.orig(LungBcell_cor, colors = c("white","white","red","darkred"))
LungBcell.level <- as.factor(findInterval(LungBcell_cor, c(min(LungBcell_cor),0.05,0.1,0.2,max(LungBcell_cor)), rightmost.closed = FALSE, all.inside = FALSE))
LungBcell.level <- c("white","pink","red","darkred")[LungBcell.level]
LungBcell.level.dark <- col2rgb(LungBcell.level)
LungBcell.level.dark <- rgb2hsv(LungBcell.level.dark)
LungBcell.level.dark[3,] <- LungBcell.level.dark[3,]*0.5
LungBcell.level.dark <- hsv(LungBcell.level.dark[1,], LungBcell.level.dark[2,], LungBcell.level.dark[3,], 1)
LungBcell.level.dark[which(is.na(LungBcell_cor))] <- NA

LungCiliated_cor <- CellsCor[2,(13:446)]

#LungCiliated.level <- cRamp.orig(LungCiliated_cor, colors = c("white","white","red","darkred"))
LungCiliated.level <- as.factor(findInterval(LungCiliated_cor, c(min(LungCiliated_cor),max(LungCiliated_cor)), rightmost.closed = FALSE, all.inside = FALSE))
LungCiliated.level <- c("white")[LungCiliated.level]
LungCiliated.level.dark <- col2rgb(LungCiliated.level)
LungCiliated.level.dark <- rgb2hsv(LungCiliated.level.dark)
LungCiliated.level.dark[3,] <- LungCiliated.level.dark[3,]*0.5
LungCiliated.level.dark <- hsv(LungCiliated.level.dark[1,], LungCiliated.level.dark[2,], LungCiliated.level.dark[3,], 1)
LungCiliated.level.dark[which(is.na(LungCiliated_cor))] <- NA

LungClassicMonocyte_cor <- CellsCor[3,(13:446)]

#LungClassicMonocyte.level <- cRamp.orig(LungClassicMonocyte_cor, colors = c("white","white","red","darkred"))
LungClassicMonocyte.level <- as.factor(findInterval(LungClassicMonocyte_cor, c(min(LungClassicMonocyte_cor),0.05,0.1,0.2,max(LungClassicMonocyte_cor)), rightmost.closed = FALSE, all.inside = FALSE))
LungClassicMonocyte.level <- c("white","pink","red","darkred")[LungClassicMonocyte.level]
LungClassicMonocyte.level.dark <- col2rgb(LungClassicMonocyte.level)
LungClassicMonocyte.level.dark <- rgb2hsv(LungClassicMonocyte.level.dark)
LungClassicMonocyte.level.dark[3,] <- LungClassicMonocyte.level.dark[3,]*0.5
LungClassicMonocyte.level.dark <- hsv(LungClassicMonocyte.level.dark[1,], LungClassicMonocyte.level.dark[2,], LungClassicMonocyte.level.dark[3,], 1)
LungClassicMonocyte.level.dark[which(is.na(LungClassicMonocyte_cor))] <- NA

LungEndothelial_cor <- CellsCor[4,(13:446)]

#LungEndothelial.level <- cRamp.orig(LungEndothelial_cor, colors = c("white","white","red","darkred"))
LungEndothelial.level <- as.factor(findInterval(LungEndothelial_cor, c(min(LungEndothelial_cor),max(LungEndothelial_cor)), rightmost.closed = FALSE, all.inside = FALSE))
LungEndothelial.level <- c("white")[LungEndothelial.level]
LungEndothelial.level.dark <- col2rgb(LungEndothelial.level)
LungEndothelial.level.dark <- rgb2hsv(LungEndothelial.level.dark)
LungEndothelial.level.dark[3,] <- LungEndothelial.level.dark[3,]*0.5
LungEndothelial.level.dark <- hsv(LungEndothelial.level.dark[1,], LungEndothelial.level.dark[2,], LungEndothelial.level.dark[3,], 1)
LungEndothelial.level.dark[which(is.na(LungEndothelial_cor))] <- NA

LungEpithelial_cor <- CellsCor[5,(13:446)]

#LungEpithelial.level <- cRamp.orig(LungEpithelial_cor, colors = c("white","white","red","darkred"))
LungEpithelial.level <- as.factor(findInterval(LungEpithelial_cor, c(min(LungEpithelial_cor),0.05,0.1,0.2,max(LungEpithelial_cor)), rightmost.closed = FALSE, all.inside = FALSE))
LungEpithelial.level <- c("white","pink","red","darkred")[LungEpithelial.level]
LungEpithelial.level.dark <- col2rgb(LungEpithelial.level)
LungEpithelial.level.dark <- rgb2hsv(LungEpithelial.level.dark)
LungEpithelial.level.dark[3,] <- LungEpithelial.level.dark[3,]*0.5
LungEpithelial.level.dark <- hsv(LungEpithelial.level.dark[1,], LungEpithelial.level.dark[2,], LungEpithelial.level.dark[3,], 1)
LungEpithelial.level.dark[which(is.na(LungEpithelial_cor))] <- NA

LungLeukocyte_cor <- CellsCor[6,(13:446)]

#LungLeukocyte.level <- cRamp.orig(LungLeukocyte_cor, colors = c("white","white","red","darkred"))
LungLeukocyte.level <- as.factor(findInterval(LungLeukocyte_cor, c(min(LungLeukocyte_cor),0.05,0.1,0.2,max(LungLeukocyte_cor)), rightmost.closed = FALSE, all.inside = FALSE))
LungLeukocyte.level <- c("white","pink","red","darkred")[LungLeukocyte.level]
LungLeukocyte.level.dark <- col2rgb(LungLeukocyte.level)
LungLeukocyte.level.dark <- rgb2hsv(LungLeukocyte.level.dark)
LungLeukocyte.level.dark[3,] <- LungLeukocyte.level.dark[3,]*0.5
LungLeukocyte.level.dark <- hsv(LungLeukocyte.level.dark[1,], LungLeukocyte.level.dark[2,], LungLeukocyte.level.dark[3,], 1)
LungLeukocyte.level.dark[which(is.na(LungLeukocyte_cor))] <- NA

LungMonocyte_cor <- CellsCor[7,(13:446)]

#LungMonocyte.level <- cRamp.orig(LungMonocyte_cor, colors = c("white","white","red","darkred"))
LungMonocyte.level <- as.factor(findInterval(LungMonocyte_cor, c(min(LungMonocyte_cor),0.05,0.1,0.2,max(LungMonocyte_cor)), rightmost.closed = FALSE, all.inside = FALSE))
LungMonocyte.level <- c("white","pink","red","darkred")[LungMonocyte.level]
LungMonocyte.level.dark <- col2rgb(LungMonocyte.level)
LungMonocyte.level.dark <- rgb2hsv(LungMonocyte.level.dark)
LungMonocyte.level.dark[3,] <- LungMonocyte.level.dark[3,]*0.5
LungMonocyte.level.dark <- hsv(LungMonocyte.level.dark[1,], LungMonocyte.level.dark[2,], LungMonocyte.level.dark[3,], 1)
LungMonocyte.level.dark[which(is.na(LungMonocyte_cor))] <- NA

LungMyeloid_cor <- CellsCor[8,(13:446)]

#LungMyeloid.level <- cRamp.orig(LungMyeloid_cor, colors = c("white","white","red","darkred"))
LungMyeloid.level <- as.factor(findInterval(LungMyeloid_cor, c(min(LungMyeloid_cor),0.05,0.1,max(LungMyeloid_cor)), rightmost.closed = FALSE, all.inside = FALSE))
LungMyeloid.level <- c("white","pink","red")[LungMyeloid.level]
LungMyeloid.level.dark <- col2rgb(LungMyeloid.level)
LungMyeloid.level.dark <- rgb2hsv(LungMyeloid.level.dark)
LungMyeloid.level.dark[3,] <- LungMyeloid.level.dark[3,]*0.5
LungMyeloid.level.dark <- hsv(LungMyeloid.level.dark[1,], LungMyeloid.level.dark[2,], LungMyeloid.level.dark[3,], 1)
LungMyeloid.level.dark[which(is.na(LungMyeloid_cor))] <- NA

LungNKcell_cor <- CellsCor[9,(13:446)]

#LungNKcell.level <- cRamp.orig(LungNKcell_cor, colors = c("white","white","red","darkred"))
LungNKcell.level <- as.factor(findInterval(LungNKcell_cor, c(min(LungNKcell_cor),0.05,0.1,max(LungNKcell_cor)), rightmost.closed = FALSE, all.inside = FALSE))
LungNKcell.level <- c("white","pink","red")[LungNKcell.level]
LungNKcell.level.dark <- col2rgb(LungNKcell.level)
LungNKcell.level.dark <- rgb2hsv(LungNKcell.level.dark)
LungNKcell.level.dark[3,] <- LungNKcell.level.dark[3,]*0.5
LungNKcell.level.dark <- hsv(LungNKcell.level.dark[1,], LungNKcell.level.dark[2,], LungNKcell.level.dark[3,], 1)
LungNKcell.level.dark[which(is.na(LungNKcell_cor))] <- NA

LungStromal_cor <- CellsCor[10,(13:446)]

#LungStromal.level <- cRamp.orig(LungStromal_cor, colors = c("white","white","red","darkred"))
LungStromal.level <- as.factor(findInterval(LungStromal_cor, c(min(LungStromal_cor),max(LungStromal_cor)), rightmost.closed = FALSE, all.inside = FALSE))
LungStromal.level <- c("white")[LungStromal.level]
LungStromal.level.dark <- col2rgb(LungStromal.level)
LungStromal.level.dark <- rgb2hsv(LungStromal.level.dark)
LungStromal.level.dark[3,] <- LungStromal.level.dark[3,]*0.5
LungStromal.level.dark <- hsv(LungStromal.level.dark[1,], LungStromal.level.dark[2,], LungStromal.level.dark[3,], 1)
LungStromal.level.dark[which(is.na(LungStromal_cor))] <- NA

LungTcell_cor <- CellsCor[11,(13:446)]

#LungTcell.level <- cRamp.orig(LungTcell_cor, colors = c("white","white","red","darkred"))
LungTcell.level <- as.factor(findInterval(LungTcell_cor, c(min(LungTcell_cor),0.05,0.1,0.2,max(LungTcell_cor)), rightmost.closed = FALSE, all.inside = FALSE))
LungTcell.level <- c("white","pink","red","darkred")[LungTcell.level]
LungTcell.level.dark <- col2rgb(LungTcell.level)
LungTcell.level.dark <- rgb2hsv(LungTcell.level.dark)
LungTcell.level.dark[3,] <- LungTcell.level.dark[3,]*0.5
LungTcell.level.dark <- hsv(LungTcell.level.dark[1,], LungTcell.level.dark[2,], LungTcell.level.dark[3,], 1)
LungTcell.level.dark[which(is.na(LungTcell_cor))] <- NA

Mitochondrial_cor <- CellsCor[12,(13:446)]

#Mitochondrial.level <- cRamp.orig(Mitochondrial_cor, colors = c("white","white","red","darkred"))
Mitochondrial.level <- as.factor(findInterval(Mitochondrial_cor, c(min(Mitochondrial_cor),max(Mitochondrial_cor)), rightmost.closed = FALSE, all.inside = FALSE))
Mitochondrial.level <- c("white")[Mitochondrial.level]
Mitochondrial.level.dark <- col2rgb(Mitochondrial.level)
Mitochondrial.level.dark <- rgb2hsv(Mitochondrial.level.dark)
Mitochondrial.level.dark[3,] <- Mitochondrial.level.dark[3,]*0.5
Mitochondrial.level.dark <- hsv(Mitochondrial.level.dark[1,], Mitochondrial.level.dark[2,], Mitochondrial.level.dark[3,], 1)
Mitochondrial.level.dark[which(is.na(Mitochondrial_cor))] <- NA

LungCors <- as.data.frame(cbind(LungBcell_cor, LungCiliated_cor, LungClassicMonocyte_cor, LungEndothelial_cor, LungEpithelial_cor, LungLeukocyte_cor, LungMonocyte_cor, 
                                LungMyeloid_cor, LungNKcell_cor, LungStromal_cor, LungTcell_cor))
LungCors$Max<-pmax(LungCors$LungBcell_cor, LungCors$LungCiliated_cor, LungCors$LungClassicMonocyte_cor, LungCors$LungEndothelial_cor, LungCors$LungEpithelial_cor, 
                   LungCors$LungLeukocyte_cor, LungCors$LungMonocyte_cor, LungCors$LungMyeloid_cor, LungCors$LungNKcell_cor, LungCors$LungStromal_cor, LungCors$LungTcell_cor)
LungBcells <- rownames(LungCors)[which(LungCors$LungBcell_cor == LungCors$Max)]
LungCiliateds <- rownames(LungCors)[which(LungCors$LungCiliated_cor == LungCors$Max)]
LungClassicMonocytes <- rownames(LungCors)[which(LungCors$LungClassicMonocyte_cor == LungCors$Max)]
LungEndothelials <- rownames(LungCors)[which(LungCors$LungEndothelial_cor == LungCors$Max)]
LungEpithelials <- rownames(LungCors)[which(LungCors$LungEpithelial_cor == LungCors$Max)]
LungLeukocytes <- rownames(LungCors)[which(LungCors$LungLeukocyte_cor == LungCors$Max)]
LungMonocytes <- rownames(LungCors)[which(LungCors$LungMonocyte_cor == LungCors$Max)]
LungMyeloids <- rownames(LungCors)[which(LungCors$LungMyeloid_cor == LungCors$Max)]
LungNKcells <- rownames(LungCors)[which(LungCors$LungNKcell_cor == LungCors$Max)]
LungStromals <- rownames(LungCors)[which(LungCors$LungStromal_cor == LungCors$Max)]
LungTcells <- rownames(LungCors)[which(LungCors$LungTcell_cor == LungCors$Max)]

rpkmdata <- t(rpkmdata)
CorCellType <- rep(NA, nrow(rpkmdata))
names(CorCellType) <- rownames(rpkmdata)
CorCellType[which(names(CorCellType) %in% LungBcells)] <- "#00A08A"
CorCellType[which(names(CorCellType) %in% LungCiliateds)] <- "#F21A00"
CorCellType[which(names(CorCellType) %in% LungClassicMonocytes)] <- "#EBCC2A"
CorCellType[which(names(CorCellType) %in% LungEndothelials)] <- "#78B7C5"
CorCellType[which(names(CorCellType) %in% LungEpithelials)] <- "#7294D4"
CorCellType[which(names(CorCellType) %in% LungLeukocytes)] <- "#EBCC2A"
CorCellType[which(names(CorCellType) %in% LungMonocytes)] <- "#EBCC2A"
CorCellType[which(names(CorCellType) %in% LungMyeloids)] <- "#EBCC2A"
CorCellType[which(names(CorCellType) %in% LungNKcells)] <- "#3B9AB2"
CorCellType[which(names(CorCellType) %in% LungStromals)] <- "#D8A499"
CorCellType[which(names(CorCellType) %in% LungTcells)] <- "#F21A00"

CorCellType.dark <- col2rgb(CorCellType)
CorCellType.dark <- rgb2hsv(CorCellType.dark)
CorCellType.dark[3,] <- CorCellType.dark[3,]*0.5
CorCellType.dark <- hsv(CorCellType.dark[1,], CorCellType.dark[2,], CorCellType.dark[3,], 1)
CorCellType.dark[which(is.na(CorCellType))] <- NA

load("../WD/QCCells_Features.RData")
LU_HEK_3h_2xFeat <- CellFeat[grep('^LU_HEK_3h_2x', rownames(CellFeat)),]
LU_CAP_3h_2xFeat <- CellFeat[grep('^LU_CAP_3h_2x', rownames(CellFeat)),]
CellFeat <- rbind(LU_HEK_3h_2xFeat, LU_CAP_3h_2xFeat)

RandomForest <- rep(NA, nrow(rpkmdata))
names(RandomForest) <- rownames(rpkmdata)
RandomForest[which(names(RandomForest) %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#78B7C5")])] <- "#00A08A"
RandomForest[which(names(RandomForest) %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#EBCC2A")])] <- "#EBCC2A"
RandomForest[which(names(RandomForest) %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#F21A00")])] <- "#78B7C5"
RandomForest[which(names(RandomForest) %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#00A08A")])] <- "#7294D4"
RandomForest[which(names(RandomForest) %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#3B9AB2")])] <- "#EBCC2A"
RandomForest[which(names(RandomForest) %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#BABABA")])] <- "#EBCC2A"
RandomForest[which(names(RandomForest) %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#E1AF00")])] <- "#EBCC2A"
RandomForest[which(names(RandomForest) %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#000000")])] <- "#D8A499"
RandomForest[which(names(RandomForest) %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#7294D4")])] <- "#F21A00"

RandomForest.dark <- col2rgb(RandomForest)
RandomForest.dark <- rgb2hsv(RandomForest.dark)
RandomForest.dark[3,] <- RandomForest.dark[3,]*0.5
RandomForest.dark <- hsv(RandomForest.dark[1,], RandomForest.dark[2,], RandomForest.dark[3,], 1)
RandomForest.dark[which(is.na(RandomForest))] <- NA

HEK <- length(which(CellFeat$CellType == "#CCCCCC"))
CAP <- length(which(CellFeat$CellType == "#666666"))

Bcells <- length(LungBcells)
aiBcells <- length(which(CellFeat$RandomForest == "#78B7C5"))
BcellsHEK <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% LungBcells))
BcellsCAP <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% LungBcells))
BcellsHEKai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#78B7C5")]))
BcellsCAPai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#78B7C5")]))

Ciliateds <- length(LungCiliateds)
aiCiliateds <- length(which(CellFeat$RandomForest == "#"))
CiliatedsHEK <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% LungCiliateds))
CiliatedsCAP <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% LungCiliateds))
CiliatedsHEKai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#")]))
CiliatedsCAPai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#")]))

ClassicMonocytes <- length(LungClassicMonocytes)
aiClassicMonocytes <- length(which(CellFeat$RandomForest == "#EBCC2A"))
ClassicMonocytesHEK <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% LungClassicMonocytes))
ClassicMonocytesCAP <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% LungClassicMonocytes))
ClassicMonocytesHEKai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#EBCC2A")]))
ClassicMonocytesCAPai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#EBCC2A")]))

Endothelials <- length(LungEndothelials)
aiEndothelials <- length(which(CellFeat$RandomForest == "#F21A00"))
EndothelialsHEK <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% LungEndothelials))
EndothelialsCAP <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% LungEndothelials))
EndothelialsHEKai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#F21A00")]))
EndothelialsCAPai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#F21A00")]))

Epithelials <- length(LungEpithelials)
aiEpithelials <- length(which(CellFeat$RandomForest == "#00A08A"))
EpithelialsHEK <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% LungEpithelials))
EpithelialsCAP <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% LungEpithelials))
EpithelialsHEKai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#00A08A")]))
EpithelialsCAPai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#00A08A")]))

Leukocytes <- length(LungLeukocytes)
aiLeukocytes <- length(which(CellFeat$RandomForest == "#3B9AB2"))
LeukocytesHEK <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% LungLeukocytes))
LeukocytesCAP <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% LungLeukocytes))
LeukocytesHEKai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#3B9AB2")]))
LeukocytesCAPai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#3B9AB2")]))

Monocytes <- length(LungMonocytes)
aiMonocytes <- length(which(CellFeat$RandomForest == "#BABABA"))
MonocytesHEK <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% LungMonocytes))
MonocytesCAP <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% LungMonocytes))
MonocytesHEKai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#BABABA")]))
MonocytesCAPai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#BABABA")]))

Myeloids <- length(LungMyeloids)
aiMyeloids <- length(which(CellFeat$RandomForest == "#E1AF00"))
MyeloidsHEK <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% LungMyeloids))
MyeloidsCAP <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% LungMyeloids))
MyeloidsHEKai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#E1AF00")]))
MyeloidsCAPai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#E1AF00")]))

NKcells <- length(LungNKcells)
aiNKcells <- length(which(CellFeat$RandomForest == "#"))
NKcellsHEK <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% LungNKcells))
NKcellsCAP <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% LungNKcells))
NKcellsHEKai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#")]))
NKcellsCAPai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#")]))

Stromals <- length(LungStromals)
aiStromals <- length(which(CellFeat$RandomForest == "#000000"))
StromalsHEK <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% LungStromals))
StromalsCAP <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% LungStromals))
StromalsHEKai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#000000")]))
StromalsCAPai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#000000")]))

Tcells <- length(LungTcells)
aiTcells <- length(which(CellFeat$RandomForest == "#7294D4"))
TcellsHEK <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% LungTcells))
TcellsCAP <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% LungTcells))
TcellsHEKai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#CCCCCC")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#7294D4")]))
TcellsCAPai <- length(which(rownames(CellFeat)[which(CellFeat$CellType == "#666666")] %in% rownames(CellFeat)[which(CellFeat$RandomForest == "#7294D4")]))

Numbers <- as.data.frame(c(HEK,CAP, Bcells, BcellsHEK, BcellsCAP, Ciliateds, CiliatedsHEK, CiliatedsCAP, ClassicMonocytes, ClassicMonocytesHEK, ClassicMonocytesCAP, Endothelials, 
                           EndothelialsHEK, EndothelialsCAP, Epithelials, EpithelialsHEK, EpithelialsCAP, Leukocytes, LeukocytesHEK, LeukocytesCAP, Monocytes, MonocytesHEK, 
                           MonocytesCAP, Myeloids, MyeloidsHEK, MyeloidsCAP, NKcells, NKcellsHEK, NKcellsCAP, Stromals, StromalsHEK, StromalsCAP, Tcells, TcellsHEK, TcellsCAP))
aiNumbers <- as.data.frame(c(HEK,CAP, aiBcells, BcellsHEKai, BcellsCAPai, aiCiliateds, CiliatedsHEKai, CiliatedsCAPai, aiClassicMonocytes, ClassicMonocytesHEKai, ClassicMonocytesCAPai, aiEndothelials, 
                             EndothelialsHEKai, EndothelialsCAPai, aiEpithelials, EpithelialsHEKai, EpithelialsCAPai, aiLeukocytes, LeukocytesHEKai, LeukocytesCAPai, aiMonocytes, MonocytesHEKai, 
                             MonocytesCAPai, aiMyeloids, MyeloidsHEKai, MyeloidsCAPai, aiNKcells, NKcellsHEKai, NKcellsCAPai, aiStromals, StromalsHEKai, StromalsCAPai, aiTcells, TcellsHEKai, TcellsCAPai))
rownames(Numbers) <- c("HEK","CAP", "Bcells", "BcellsHEK", "BcellsCAP", "Ciliateds", "CiliatedsHEK", "CiliatedsCAP", "ClassicMonocytes", "ClassicMonocytesHEK", "ClassicMonocytesCAP", "Endothelials", 
                       "EndothelialsHEK", "EndothelialsCAP", "Epithelials", "EpithelialsHEK", "EpithelialsCAP", "Leukocytes", "LeukocytesHEK", "LeukocytesCAP", "Monocytes", "MonocytesHEK", 
                       "MonocytesCAP", "Myeloids", "MyeloidsHEK", "MyeloidsCAP", "NKcells", "NKcellsHEK", "NKcellsCAP", "Stromals", "StromalsHEK", "StromalsCAP", "Tcells", "TcellsHEK", "TcellsCAP")
rownames(aiNumbers) <- c("HEK","CAP", "Bcells", "BcellsHEK", "BcellsCAP", "Ciliateds", "CiliatedsHEK", "CiliatedsCAP", "ClassicMonocytes", "ClassicMonocytesHEK", "ClassicMonocytesCAP", "Endothelials", 
                         "EndothelialsHEK", "EndothelialsCAP", "Epithelials", "EpithelialsHEK", "EpithelialsCAP", "Leukocytes", "LeukocytesHEK", "LeukocytesCAP", "Monocytes", "MonocytesHEK", 
                         "MonocytesCAP", "Myeloids", "MyeloidsHEK", "MyeloidsCAP", "NKcells", "NKcellsHEK", "NKcellsCAP", "Stromals", "StromalsHEK", "StromalsCAP", "Tcells", "TcellsHEK", "TcellsCAP")
write.csv(Numbers, file="LungNumbers3h.csv")
write.csv(aiNumbers, file="aiLungNumbers3h.csv")

###################################################################

load("../WD/CD63lung3hCellsQC_tSNE-100noid-3nod-200non_983100TabulaGenes_tSNEresults.Rdata")
Cluster <- rep(NA, nrow(rpkmdata))
names(Cluster) <- rownames(rpkmdata)
Cluster1 <- names(which(tSNEresults$g.info[[1]] == "grey"))
#Cluster2 <- names(which(tSNEresults$g.info[[1]] == "orange"))
Cluster3 <- names(which(tSNEresults$g.info[[1]] == "yellow"))

Cluster[which(names(Cluster) %in% Cluster1)] <- "#00A08A"
#Cluster[which(names(Cluster) %in% Cluster2)] <- "#EBCC2A"
Cluster[which(names(Cluster) %in% Cluster3)] <- "#7294D4"

Cluster.dark <- col2rgb(Cluster)
Cluster.dark <- rgb2hsv(Cluster.dark)
Cluster.dark[3,] <- Cluster.dark[3,]*0.5
Cluster.dark <- hsv(Cluster.dark[1,], Cluster.dark[2,], Cluster.dark[3,], 1)
Cluster.dark[which(is.na(Cluster))] <- NA

### Plot data
pdf("Lung3hMetaCells.pdf", width=10, height=10)
par(mfrow=c(4,4), mar=c(1,1,1,1))

load("../WD/CD63lung3hCellsQC_tSNE-100noid-3nod-40non_983100TabulaGenes_tSNEresults.Rdata")
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
adjacency.plot2(g, main = sprintf("LungBcellMetaCell %s", ncol(data)), LungBcell.level, LungBcell.level.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("LungCiliatedMetaCell %s", ncol(data)), LungCiliated.level, LungCiliated.level.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("LungClassicMonocyteMetaCell %s", ncol(data)), LungClassicMonocyte.level, LungClassicMonocyte.level.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("LungEndothelialMetaCell %s", ncol(data)), LungEndothelial.level, LungEndothelial.level.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("LungEpithelialMetaCell %s", ncol(data)), LungEpithelial.level, LungEpithelial.level.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("LungLeukocyteMetaCell %s", ncol(data)), LungLeukocyte.level, LungLeukocyte.level.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("LungMonocyteMetaCell %s", ncol(data)), LungMonocyte.level, LungMonocyte.level.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("LungMyeloidMetaCell %s", ncol(data)), LungMyeloid.level, LungMyeloid.level.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("LungNKcellMetaCell %s", ncol(data)), LungNKcell.level, LungNKcell.level.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("LungStromalMetaCell %s", ncol(data)), LungStromal.level, LungStromal.level.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("LungTcellMetaCell %s", ncol(data)), LungTcell.level, LungTcell.level.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("Mitochondrial %s", ncol(data)), Mitochondrial.level, Mitochondrial.level.dark, g.layout,6)


### Plot expression data
rpkmdata<-read.table("../WD/Lu_3h_2x_all_batch_correction_norm_express.txt", header=T)

#LU_HEK_30m_1x <- rpkmdata[,grep('^LU_HEK_30m_1x', colnames(rpkmdata))]
#LI_HEK_30m_2x <- rpkmdata[,grep('^LI_HEK_30m_2x', colnames(rpkmdata))]
LU_HEK_3h_2x <- rpkmdata[,grep('^LU_HEK_3h_2x', colnames(rpkmdata))]
#LU_CAP_30m_1x <- rpkmdata[,grep('^LU_CAP_30m_1x', colnames(rpkmdata))]
#LI_CAP_30m_2x <- rpkmdata[,grep('^LI_CAP_30m_2x', colnames(rpkmdata))]
LU_CAP_3h_2x <- rpkmdata[,grep('^LU_CAP_3h_2x', colnames(rpkmdata))]

rpkmdata <- cbind(LU_HEK_3h_2x, LU_CAP_3h_2x)

data <- t(rpkmdata)
data[data==0] <- 0.0001
data <- log2(data)

w<-c("ENSMUSG00000069516", "ENSMUSG00000030724","ENSMUSG00000044206","ENSMUSG00000020884","ENSMUSG00000024653")#LYZ2, CD19, VSIG4, ASGR1

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

#LungClusters <- as.data.frame(cbind(Cluster, LungBcell_cor, LungCiliated_cor, LungClassicMonocyte_cor, LungEndothelial_cor, LungEpithelial_cor,
#                                    LungLeukocyte_cor, LungMonocyte_cor, LungMyeloid_cor, LungNKcell_cor, LungStromal_cor, LungTcell_cor))
#LungCluster1 <- LungClusters[which(LungClusters$Cluster == "#00A08A"),]
#LungCluster2 <- LungClusters[which(LungClusters$Cluster == "#EBCC2A"),]
#LungCluster3 <- LungClusters[which(LungClusters$Cluster == "#7294D4"),]
#write.csv(LungCluster1, "LungClusterGreen.csv")
#write.csv(LungCluster2, "LungClusterYellow.csv")
#write.csv(LungCluster3, "LungClusterPurple.csv")
