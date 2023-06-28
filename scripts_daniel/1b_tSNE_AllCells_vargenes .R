setwd("~/data/CD63GFP/WD/")
output <- "../WD/"

args<-commandArgs(TRUE)

#BiocManager::install('WGCNA')
#BiocManager::install('Rtsne')
#BiocManager::install('fields')
#BiocManager::install('Matrix')
#BiocManager::install('igraph')
#BiocManager::install('gplots')
#BiocManager::install('colorspace')
#BiocManager::install('psych')
#BiocManager::install('rdist')
#BiocManager::install('flashClust')
#BiocManager::install('Cairo')
#BiocManager::install('scde')

source("../scripts/Functions/20170228_PCA_functions.R")
source("../scripts/Functions/20170228_Standard_Bioinformatics_functions.R")
source("../scripts/Functions/20170228_tSNE_Graphing_functions.R")
source("../scripts/Functions/20170228_WGCNA_SCDE.R")
source("../scripts/Functions/20160428_tSNE_Stouffers_functions.R")

library(Rtsne)
library(fields)
library(Matrix)
library(igraph)
library(gplots)
library(colorspace)
library(psych)

load("../WD/QCCells_Features.RData")

#cen   <- 0.9
non   <- 200 #as.numeric(args[1])
nod   <- 3 # as.numeric(args[2])
noid  <- 100
plotgenes <- TRUE
colourPalette <- colorRampPalette(c("white", "cornsilk2", "darkgoldenrod1", "darkorange4", "darkred", "black"))(20)
window <- 10000

### load the data frame
load("../WD/SingleCells_ensemble_rpkm_QC.RData")
data <- t(rpkmdata)
data[data==0] <- 0.0001
data <- log2(data)

var.rpkms<-var.genes(data)
data <- subset(data, select = colnames(data) %in% var.rpkms)

name <- sprintf("CD63CellsQCharm_tSNE-%snoid-%snod-%snon_%sVargenes", noid, nod, non, ncol(data))

### tSNE Plotting
tsne_out <- Rtsne(data, dims= nod, initial_dims = noid, theta = 0.001, pca = FALSE)
rownames(tsne_out$Y) <- rownames(data)
orig.coord <- tsne_out$Y

d<-rdist(orig.coord)
colnames(d)<-rownames(orig.coord)
rownames(d)<-rownames(orig.coord)

#adj.mat <- connect.cen.max(d, cen, non)
adj.mat<-connect.neighbours(d, non)
rownames(adj.mat)<-colnames(adj.mat)

g <- graph.adjacency(adj.mat, mode = "directed", weighted=TRUE, diag=FALSE)
g.layout <- layout.fruchterman.reingold(g, weights=E(g)$weight) #coord = orig.coord
g.info <- InfoMap.Community.Colours(g)
g.walktrap.one <- Walktrap.Community.Colours(g, 1)
g.walktrap.two <- Walktrap.Community.Colours(g, 2)
g.walktrap.three <- Walktrap.Community.Colours(g, 3)
g.walktrap.four <- Walktrap.Community.Colours(g, 4)

### WGCNA analysis
WGCNA <- WGCNA.Workflow(data, "unsigned", 5, "p", "a", 5, 4, 0.1, colourPalette, TRUE)
tSNEresults <- list(data, tsne_out, orig.coord, adj.mat, g, g.layout, g.info, WGCNA)
names(tSNEresults) <- c("data", "tsne_out","orig.coord", "adj.mat", "g", "g.layout", "g.info", "WGCNAresults")
save(tSNEresults, file = sprintf("%s_tSNEresults.Rdata", name))

#load("../WD/AllCD63CellsQC_tSNE-100noid-3nod-200non_1573Vargenes_tSNEresults.Rdata")
#g <- tSNEresults$g
#g.info <- tSNEresults$g.info
#g.layout <- tSNEresults$g.layout

### Plot data
pdf(sprintf("%s.pdf", name), width=20, height=20)
par(mfrow=c(4,4), mar=c(1,1,1,1))

adjacency.plot2(g, main = sprintf("Organ %s", ncol(data)), CellFeat$Organ, CellFeat$Organ.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("CellType %s", ncol(data)), CellFeat$CellType, CellFeat$CellType.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("TimePoint %s", ncol(data)), CellFeat$TimePoint, CellFeat$TimePoint.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("Dose %s", ncol(data)), CellFeat$Dose, CellFeat$Dose.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("DetGenes %s", ncol(data)), CellFeat$DetGenes, CellFeat$DetGenes.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("ReadCount %s", ncol(data)), CellFeat$ReadCount, CellFeat$ReadCount.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("clusters %s", ncol(data)), CellFeat$clusters, CellFeat$clusters.dark, g.layout,6)
adjacency.plot2(g, main = sprintf("Infomap %s", ncol(data)), g.info[[1]], g.info[[2]], g.layout,6)
#adjacency.plot(g, main = sprintf("Walktrap1 %s", ncol(data)), g.walktrap.one[[1]], g.walktrap.one[[2]], g.layout)
#adjacency.plot(g, main = sprintf("Walktrap2 %s", ncol(data)), g.walktrap.two[[1]], g.walktrap.two[[2]], g.layout)
#adjacency.plot(g, main = sprintf("Walktrap3 %s", ncol(data)), g.walktrap.three[[1]], g.walktrap.three[[2]], g.layout)
#adjacency.plot(g, main = sprintf("Walktrap4 %s", ncol(data)), g.walktrap.four[[1]], g.walktrap.four[[2]], g.layout)

### Plot expression data
load("../WD/SingleCells_ensemble_rpkm_QC.RData")
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

### plot WGCNA
#AExpn <- WGCNA$MergeMEListExpression$averageExpr
#
#for (i in colnames(AExpn)){
#  exp.level <- cRamp.orig(AExpn[,i], c("white", "cornsilk2", "darkgoldenrod1", "darkorange4", "darkred", "black"))
#  exp.level <- col2rgb(exp.level)
#  exp.level <- rgb2hsv(exp.level)
#  exp.level <- hsv(exp.level[1,], exp.level[2,], exp.level[3,], 1)
#  exp.level[which(exp.level == "#FFFFFFFF")] <- rep(makeTransparent("white", 0.5), sum(exp.level == "#FFFFFFFF"))
#  
#  exp.level.dark <- col2rgb(exp.level)
#  exp.level.dark <- rgb2hsv(exp.level.dark)
#  exp.level.dark[3,] <- exp.level.dark[3,]*0.5
#  exp.level.dark <- hsv(exp.level.dark[1,], exp.level.dark[2,], exp.level.dark[3,], 1)
#  exp.level.dark[which(is.na(exp.level))] <- NA
#  
#  adjacency.plot(g, main = i, exp.level, exp.level.dark, g.layout)
#}

dev.off()

### do SCDE between similar maturation groups, followed by WGCNA
load("../WD/AllCD63Cells50kreads_tSNE-100noid-5nod-25non_1008Vargenes_tSNEresults.Rdata")
load("../WD/SingleCells_ensemblecounts_2rpkmin5cells_min50kGenes_1000-10000genes_ReadsGenesQC.RData")
load("../WD/SingleCells_ensemblerpkms_2rpkmin5cells_min50kGenes_1000-10000genes_ReadsGenesQC.RData")
ClusterPairsList <- PairsCompList(unlist(g.info[1]), CellFeat$NoDetGenes, window)
#ClusterPairsList1 <- PairsCompList(unlist(CellFeat$clusters), CellFeat$NoDetGenes, window)
ClusterPairsList1.4 <- ClusterPairsList[1:4]
ClusterPairsList5.8 <- ClusterPairsList[5:8]
ClusterPairsList9.12 <- ClusterPairsList[9:12]
ClusterPairsList13.16 <- ClusterPairsList[13:16]
ClusterPairsList17.20 <- ClusterPairsList[17:20]
ClusterPairsList21.24 <- ClusterPairsList[21:24]
ClusterPairsList25.28 <- ClusterPairsList[25:28]
ClusterPairsList29.32 <- ClusterPairsList[29:32]
ClusterPairsList33.36 <- ClusterPairsList[33:36]

InfomapSCDE1.4 <- SCDE.DiffExp.ChosenGroups(countdata, as.factor(unlist(g.info[[1]])), ClusterPairsList1.4, 10)
save(InfomapSCDE1.4, file="InfomapSCDE1.4.Rdata")
InfomapSCDE5.8 <- SCDE.DiffExp.ChosenGroups(countdata, as.factor(unlist(g.info[[1]])), ClusterPairsList5.8, 10)
save(InfomapSCDE5.8, file="InfomapSCDE5.8.Rdata")
InfomapSCDE9.12 <- SCDE.DiffExp.ChosenGroups(countdata, as.factor(unlist(g.info[[1]])), ClusterPairsList9.12, 10)
save(InfomapSCDE9.12, file="InfomapSCDE9.12.Rdata")
InfomapSCDE13.16 <- SCDE.DiffExp.ChosenGroups(countdata, as.factor(unlist(g.info[[1]])), ClusterPairsList13.16, 10)
save(InfomapSCDE13.16, file="InfomapSCDE13.16.Rdata")
InfomapSCDE17.20 <- SCDE.DiffExp.ChosenGroups(countdata, as.factor(unlist(g.info[[1]])), ClusterPairsList17.20, 10)
save(InfomapSCDE17.20, file="InfomapSCDE17.20.Rdata")
InfomapSCDE21.24 <- SCDE.DiffExp.ChosenGroups(countdata, as.factor(unlist(g.info[[1]])), ClusterPairsList21.24, 10)
save(InfomapSCDE21.24, file="InfomapSCDE21.24.Rdata")
InfomapSCDE25.28 <- SCDE.DiffExp.ChosenGroups(countdata, as.factor(unlist(g.info[[1]])), ClusterPairsList25.28, 10)
save(InfomapSCDE25.28, file="InfomapSCDE25.28.Rdata")
InfomapSCDE29.32 <- SCDE.DiffExp.ChosenGroups(countdata, as.factor(unlist(g.info[[1]])), ClusterPairsList29.32, 10)
save(InfomapSCDE29.32, file="InfomapSCDE29.32.Rdata")
InfomapSCDE33.36 <- SCDE.DiffExp.ChosenGroups(countdata, as.factor(unlist(g.info[[1]])), ClusterPairsList33.36, 10)
save(InfomapSCDE33.36, file="InfomapSCDE33.36.Rdata")

load("../WD/InfomapSCDE1.4.Rdata")
load("../WD/InfomapSCDE5.8.Rdata")
load("../WD/InfomapSCDE9.12.Rdata")
load("../WD/InfomapSCDE13.16.Rdata")
load("../WD/InfomapSCDE17.20.Rdata")
load("../WD/InfomapSCDE21.24.Rdata")
load("../WD/InfomapSCDE25.28.Rdata")
load("../WD/InfomapSCDE29.32.Rdata")
load("../WD/InfomapSCDE33.36.Rdata")

InfomapSCDE <- (c(InfomapSCDE1.4, InfomapSCDE5.8, InfomapSCDE9.12, InfomapSCDE13.16, 
                  InfomapSCDE17.20, InfomapSCDE21.24, InfomapSCDE21.24, InfomapSCDE25.28,
                  InfomapSCDE29.32, InfomapSCDE33.36))

StouffersComb <- StouffersComb(InfomapSCDE)
save(StouffersComb, file="StouffersComb.Rdata")

# WGCNA grouping of genes, removing grey genes, and plot.
data <- rpkmdata
data[data==0] <- 0.0001
data <- log2(t(data))
InfomapWGCNA <- WGCNA.Workflow(subset(data, select = colnames(data) %in% rownames(StouffersComb)[which(StouffersComb$Z > 15)]), "unsigned", 5, "p", "a", 5, 4, 0.1, TRUE, 10)

#AExpn <- InfomapWGCNA$MergeMEListExpression$averageExpr
load("../WD/AllCD63Cells50kreads_tSNE-100noid-5nod-25non_1008Vargenes_10000win-Infomap-SCDE-WGCNA-results.Rdata")

pdf(sprintf("%s_%swin-Infomap-SCDE-WGCNA-results.pdf", name, window), width=25, height=25)
par(mfrow=c(5,5), mar=c(1,1,1,1))

adjacency.plot2(g, main = sprintf("CellType %s", ncol(data)), CellFeat$CellType, CellFeat$CellType.dark, tSNEresults$g.layout,6)
adjacency.plot2(g, main = sprintf("DetGenes %s", ncol(data)), CellFeat$DetGenes, CellFeat$DetGenes.dark, tSNEresults$g.layout,6)
adjacency.plot2(g, main = sprintf("ReadCount %s", ncol(data)), CellFeat$ReadCount, CellFeat$ReadCount.dark, tSNEresults$g.layout,6)
adjacency.plot2(g, main = sprintf("Infomap %s", ncol(data)), tSNEresults$g.info[[1]], tSNEresults$g.info[[2]], tSNEresults$g.layout,6)

for (i in colnames(AExpn)){
  exp.level <- cRamp.orig(AExpn[,i], c("white", "cornsilk2", "darkgoldenrod1", "darkorange4", "darkred", "black"))
  exp.level <- col2rgb(exp.level)
  exp.level <- rgb2hsv(exp.level)
  exp.level <- hsv(exp.level[1,], exp.level[2,], exp.level[3,], 1)
  exp.level[which(exp.level == "#FFFFFFFF")] <- rep(makeTransparent("white", 0.5), sum(exp.level == "#FFFFFFFF"))
  
  exp.level.dark <- col2rgb(exp.level)
  exp.level.dark <- rgb2hsv(exp.level.dark)
  exp.level.dark[3,] <- exp.level.dark[3,]*0.5
  exp.level.dark <- hsv(exp.level.dark[1,], exp.level.dark[2,], exp.level.dark[3,], 1)
  exp.level.dark[which(is.na(exp.level))] <- NA
  
  adjacency.plot2(g, main = i, exp.level, exp.level.dark, tSNEresults$g.layout,6)
}

dev.off()

Infomapresults <- list(InfomapSCDE, InfomapWGCNA)
names(Infomapresults) <- c("InfomapSCDE", "InfomapWGCNA")
save(Infomapresults, file = sprintf("%s_%swin-Infomap-SCDE-WGCNA-results.Rdata", name, window))
write.csv(InfomapWGCNA$ResultsTable, file = sprintf("%s_%swin-Infomap-SCDE-WGCNA-results.csv", name, window))
