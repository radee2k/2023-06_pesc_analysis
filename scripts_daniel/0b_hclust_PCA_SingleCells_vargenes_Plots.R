setwd("~/Desktop/bioinformatics/R/CD63GFP/WD/")
output <- "../WD/"

source("../scripts/functions/20170228_Standard_Bioinformatics_functions.R")
library(gplots)

load("../WD/QCCells_Features.RData")
load("../WD/SingleCells_ensembleNorm_2rpkmin10cells_1sdDetgenes_1sdCellcor_QC.RData")
data <- t(normdata)
data[data==0] <- 0.0001
data <- log2(data)

var.rpkms<-var.genes(data)
data <- subset(data, select = colnames(data) %in% var.rpkms)

hc1cells<-hclust(as.dist(1-cor(t(data), method="spearman")), method="ward")

################### do some heirachial clustering
pdf(file = "AllCD63CellsQC_hclust_vargenes.pdf", width = 10, height = 10)

### color HCA by embstage and hcagroups

par(mfrow=c(4,1),cex=0.5, mar=c(5, 4, 3, 2), oma=c(3, 3, 3, 3))
colplot.clust(hc1cells,
              lab.col=CellFeat$CellType,
              cex=2,hang=0.01,
              main=sprintf(paste("Genes Spearman Clustering", paste(as.character(ncol(data)), "DE genes",sep=""),sep="")))

colplot.clust(hc1cells,
              lab.col=CellFeat$DetGenes,
              cex=2,hang=0.01,
              main=sprintf(paste("Genes Spearman Clustering", paste(as.character(ncol(data)), "DE genes",sep=""),sep="")))

colplot.clust(hc1cells,
              lab.col=CellFeat$ReadCount,
              cex=2,hang=0.01,
              main=sprintf(paste("Genes Spearman Clustering", paste(as.character(ncol(data)), "DE genes",sep=""),sep="")))

colplot.clust(hc1cells,
              lab.col=CellFeat$clusters,
              cex=2,hang=0.01,
              main=sprintf(paste("Genes Spearman Clustering", paste(as.character(ncol(data)), "DE genes",sep=""),sep="")))

##### heatmap2

CellsCor=(cor(t(data)))
colors = c(seq(0, 1,length=23))
my_palette <- colorRampPalette(c("white","red"))
matrix<-as.matrix(CellsCor)

heatmap.2(matrix,
          Rowv = as.dendrogram(hclust(as.dist(1-cor(t(data), method="spearman")), method="ward")),
          Colv = "Rowv",
          symm = TRUE,
          dendrogram="column",
          ColSideColors=CellFeat$CellType,
          trace="none",
          col=my_palette,
          key=TRUE,
          keysize=1,
          density.info="none")

heatmap.2(matrix,
          Rowv = as.dendrogram(hclust(as.dist(1-cor(t(data), method="spearman")), method="ward")),
          Colv = "Rowv",
          symm = TRUE,
          dendrogram="column",
          ColSideColors=CellFeat$clusters,
          trace="none",
          col=my_palette,
          key=TRUE,
          keysize=1,
          density.info="none")

dev.off()

pdf(file = "AllCD63CellsQC_PCA_vargenes.pdf", width = 10, height = 10)

pcadata.Data<-prcomp(data,
                     center = TRUE,
                     scale = TRUE)

pca.stage.data<-pcadata.Data$x[,c(1:5)]

pairs(pcadata.Data$x[,c(1:5)],
      line.main = 3,
      cex.labels = NULL,
      font.labels = 1,
      row1attop = TRUE,
      gap = 0.5,
      col=CellFeat$CellType, # fill color
      pch=16, # symbol type, 21= a filled circle
      cex=2) # symbol size

pairs(pcadata.Data$x[,c(1:5)],
      line.main = 3,
      cex.labels = NULL,
      font.labels = 1,
      row1attop = TRUE,
      gap = 0.5,
      col=CellFeat$DetGenes, # fill color
      pch=16, # symbol type, 21= a filled circle
      cex=2) # symbol size

pairs(pcadata.Data$x[,c(1:5)],
      line.main = 3,
      cex.labels = NULL,
      font.labels = 1,
      row1attop = TRUE,
      gap = 0.5,
      col=CellFeat$ReadCount, # fill color
      pch=16, # symbol type, 21= a filled circle
      cex=2) # symbol size

pairs(pcadata.Data$x[,c(1:5)],
      line.main = 3,
      cex.labels = NULL,
      font.labels = 1,
      row1attop = TRUE,
      gap = 0.5,
      col=CellFeat$clusters, # fill color
      pch=16, # symbol type, 21= a filled circle
      cex=2) # symbol size

#x<-1 # PCA component to select
#y<-2 # PCA component to select

#plot(pca.stage.data[,x], pca.stage.data[,y],
#     xlab = colnames(pca.stage.data)[x],
#     ylab = colnames(pca.stage.data)[y],
#     type="p", # symbol type = point
#     col="black", # line color
#     pch=21, # symbol type, 21= a filled circle
#     cex=.5,
#     lwd=1)
#identify(pca.stage.data[,x], pca.stage.data[,y], plot=TRUE)

#rownames(data)[20]

# calcualte the variance
vars<- pcadata.Data$sdev^2
vars<- vars/sum(vars)

# save the variance
par(mfcol=c(1,1), oma=c(2,2,2,2), mar=c(3,3,3,3))
barplot(vars[1:5],
        main = "PCA Explained Variance",
        names.arg=colnames(pcadata.Data$rotation)[1:5],
        cex.names=2,
        las=2)

# start a dataframe with the genes ordered according to the first component
PCA<-as.data.frame(pcadata.Data$rotation[,c(1:5)])
PCA.names<-rownames(PCA)
rownames(PCA)<-NULL
PCA.names<-cbind(as.data.frame(PCA.names), PCA)

pca.rank.name<-PCA.names[,c(1,2)]
order<-order(pca.rank.name[,2], decreasing=FALSE)
PCA.ordered<-pca.rank.name[order,]
rownames(PCA.ordered)<-NULL

# start a loop to fill the rest of the dataframe with 10 PCs altogether
for (i in 3:ncol(PCA.names))
{
  pca.rank.name<-PCA.names[,c(1,i)]
  order<-order(pca.rank.name[,2], decreasing=FALSE)
  pca.rank.name<-pca.rank.name[order,]
  rownames(pca.rank.name)<-NULL
  PCA.ordered<-cbind(PCA.ordered, pca.rank.name)
}

par(mfcol=c(2,6),mar=c(3,3,3,3), oma=c(3,3,3,3))
for (i in 2*(1:5)){
  
  barplot(rev(PCA.ordered[c(1:40),i]),
          horiz=TRUE,
          main =paste(as.character(colnames(PCA.ordered)[i]), as.character(round(pcadata.Data$sdev[i/2], 2)), sep=" "),
          cex.main=2,
          width=0.5,
          space=0.75,
          names.arg=rev(PCA.ordered[c(1:40),i-1]),
          cex.names=0.75,
          las=2)
  
  barplot(PCA.ordered[c(nrow(PCA.ordered):(nrow(PCA.ordered)-40)),i],
          horiz=TRUE,
          width=0.5,
          space=0.75,
          names.arg=PCA.ordered[c(nrow(PCA.ordered):(nrow(PCA.ordered)-40)),i-1],
          cex.names=0.75,
          las=2)
}

dev.off()

# save this dataframe
save(PCA.ordered, file=paste("../WD/AllCD63CellsQCPCAComponents_vargenes.Rdata" , sep=""))
write.csv(PCA.ordered, file=paste("../WD/AllCD63CellsQCPCAComponents_vargenes.csv" , sep=""))
