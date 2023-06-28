#!/bin/R




library(Rtsne)
library(org.Mm.eg.db)
library(DESeq2)
library(edgeR)
library(gplots)
library(plotrix)
library(monocle)


source("scrna_lib.r")

express_matrix <- read.csv("express_matrix.txt", header = T, sep = "\t")
gene_id <- read.csv("genes_title.txt", header = T, sep = "\t")
rownames(express_matrix) <- gene_id[,1]

######################## deseq2 normalization ########################

normalized_matrix <- deseq2_normalizaton(express_matrix, 2, ncol(express_matrix)-2)


######################## cal rpkm ########################

gene_length <- read.csv("genes_length.txt", header = T, sep = "\t")

d <- DGEList(counts=express_matrix)
#cpm(d)
#cpm(d,log=TRUE)

d$genes$Length <- gene_length[,1]

rpkm_express_matrix <- rpkm(d)
gene_id[,1] -> gene.name

RPKM <- rpkm_express_matrix # HFD and ND

NORM <- normalized_matrix
colnames(NORM) <- colnames(RPKM)
nS <- ncol(NORM) # number of samples
rownames(RPKM)<- gene.name


######################## plot detected genes ########################

ercc <- grep("ERCC",rownames(RPKM))
nD.ercc <- apply(RPKM[ercc,],2,detect.genes)

nD.tx <- apply(RPKM[-ercc,],2,detect.genes) ## rpkm > 1
ercc.ratio<-colSums(RPKM[-ercc,])/colSums(RPKM[ercc,])

par(mfrow=c(2,1),mar=c(3,3,1,1),cex=0.6)
hist(nD.tx,n=100,main="Detected transcripts")
hist(nD.ercc,n=100,main="Detected ERCC")


#################### remove the high mito expressed cells ################

cc <- colorRampPalette(c("yellow","red","black"))
col.nDet <- convert.to.color(nD.tx,cc)


mit <- read.csv("mouse_mitochondrial_genes_list.txt", header = T, sep = "\t")
match_index_mit <- match(mit[,1], rownames(RPKM))

nD.mit <- colSums(NORM[match_index_mit,]) / colSums(NORM[-ercc,])

plot(nD.tx, nD.mit, pch = 16, col = col.nDet$cols)

low_mito <- names(which(nD.mit <= 0.02))

RPKM <- RPKM[,low_mito]
NORM <- NORM[,low_mito]

ercc <- grep("ERCC",rownames(RPKM))
nD.ercc <- apply(RPKM[ercc,],2,detect.genes)

nD.tx <- apply(RPKM[-ercc,],2,detect.genes) ## rpkm > 1
ercc.ratio<-colSums(RPKM[-ercc,])/colSums(RPKM[ercc,])




par(mfrow=c(2,1),mar=c(3,3,1,1),cex=0.6)
hist(nD.tx,n=100,main="Detected transcripts")
hist(nD.ercc,n=100,main="Detected ERCC")


######################## cutoff at 6000 genes to filter out cells ##########################

nG<-nrow(RPKM)
nS<-ncol(RPKM)

passed<-which(nD.tx >= 100 & nD.tx<=1500) ############# need changes ##############

length(passed)


length(passed)/nS


NORM <- NORM[-ercc,passed]
RPKM <- RPKM[-ercc,passed]

##################### plot PCA ####################

cc <- colorRampPalette(c("yellow","red","black"))
col.nDet <- convert.to.color(nD.tx,cc)



PC <- run.pca(NORM)
pca.plot(PC,pch=16,main="PCA, color by detected genes",col=col.nDet$cols)
#text(PC$x[, c("PC1" ,"PC2")], labels = colnames(NORM), cex = 0.3, col = "magenta")

eigs <- PC$sdev^2
barplot(eigs / sum(eigs) * 100, main="Percentage of Principle Components", ylab = "Percentage", col = "magenta", border = F)





# why always get the number of PC = number of cells, instead of number of genes?
# In fact, prcomp does this truncation internally -- in stats:::prcomp.default, if (rank < ncol(x)) ..., that just truncates the svd
# In fact, prcomp does this truncation internally -- in stats:::prcomp.default, just truncate the top cells number of PC.

################## plot tSNE #####################

tsne2 <- Rtsne(t(log2(NORM+1)),  perplexity = 25) ############# need changes ##############

par(mfrow=c(2,2),mar=c(2,2,2,2),oma=c(1,1,1,1))

plot(tsne2$Y,col=col.nDet$cols,pch=16)

# test hclust on the tSNE matrix

rownames(tsne2$Y) <- substr(colnames(NORM), 1,5)
hc2 <- hclust(dist(tsne2$Y),method="ward.D2")
groups.cl <- cutree(hc2,5) ############# need changes ##############

hc2_cluster <- rect.hclust(hc2, k=5) # define groups ############# need changes ##############











clustcol_1 <- c( "black", "red", "green3", "blue","cyan", "magenta",  "yellow", "grey", "brown", "purple",  "pink",  "navy", "gold", "orange", "orangered",  "blueviolet", "beige", "deepskyblue", "darkgreen", "darkmagenta","bisque","maroon",  "turquoise", "chartreuse",   "violet",  "plum", "antiquewhite1","antiquewhite2","antiquewhite3","antiquewhite4","aquamarine1","aquamarine2","aquamarine3","aquamarine4","azure1","azure2","azure3","azure4","bisque1","bisque2","bisque3","bisque4","blue1","blue2","blue3","blue4","brown1","brown2","brown3","brown4","burlywood1","burlywood2","burlywood3","burlywood4","cadetblue1","cadetblue2","cadetblue3","cadetblue4","chartreuse1","chartreuse2","chartreuse3","chartreuse4","chocolate1","chocolate2","chocolate3","chocolate4","coral1","coral2","coral3","coral4","cornsilk1","cornsilk2","cornsilk3","cornsilk4","cyan1","cyan2","cyan3","cyan4","darkgoldenrod1","darkgoldenrod2","darkgoldenrod3","darkgoldenrod4","darkolivegreen1","darkolivegreen2","darkolivegreen3","darkolivegreen4","darkorange1","darkorange2","darkorange3","darkorange4","darkorchid1","darkorchid2","darkorchid3","darkorchid4","darkseagreen1","darkseagreen2","darkseagreen3","darkseagreen4","darkslategray1","darkslategray2","darkslategray3","darkslategray4","deeppink1","deeppink2","deeppink3","deeppink4","deepskyblue1","deepskyblue2","deepskyblue3","deepskyblue4")

clustcol_2 <- clustcol_1[1:length(hc2_cluster)]

groups.cl2 <- mgsub(c(paste0("^",seq(1:5),"$")), clustcol_2, groups.cl) ############### define groups

colplot.clust(hc2,lab.col=as.character(groups.cl2))

plot(tsne2$Y,col=as.character(groups.cl2),pch=16)


############# adjust to 7 clusters ###############

hc2 <- hclust(dist(tsne2$Y),method="ward.D2")
groups.cl <- cutree(hc2,7) ############### define groups

hc2_cluster <- rect.hclust(hc2, k=7) ############### define groups

#clustcol_2 <- palette()[1:length(hc2_cluster)]
clustcol_1 <- c( "black", "red", "green3", "blue","cyan", "magenta",  "yellow", "grey", "brown", "purple",  "pink",  "navy", "gold", "orange", "orangered",  "blueviolet", "beige", "deepskyblue", "darkgreen", "darkmagenta","bisque","maroon",  "turquoise", "chartreuse",   "violet",  "plum", "antiquewhite1","antiquewhite2","antiquewhite3","antiquewhite4","aquamarine1","aquamarine2","aquamarine3","aquamarine4","azure1","azure2","azure3","azure4","bisque1","bisque2","bisque3","bisque4","blue1","blue2","blue3","blue4","brown1","brown2","brown3","brown4","burlywood1","burlywood2","burlywood3","burlywood4","cadetblue1","cadetblue2","cadetblue3","cadetblue4","chartreuse1","chartreuse2","chartreuse3","chartreuse4","chocolate1","chocolate2","chocolate3","chocolate4","coral1","coral2","coral3","coral4","cornsilk1","cornsilk2","cornsilk3","cornsilk4","cyan1","cyan2","cyan3","cyan4","darkgoldenrod1","darkgoldenrod2","darkgoldenrod3","darkgoldenrod4","darkolivegreen1","darkolivegreen2","darkolivegreen3","darkolivegreen4","darkorange1","darkorange2","darkorange3","darkorange4","darkorchid1","darkorchid2","darkorchid3","darkorchid4","darkseagreen1","darkseagreen2","darkseagreen3","darkseagreen4","darkslategray1","darkslategray2","darkslategray3","darkslategray4","deeppink1","deeppink2","deeppink3","deeppink4","deepskyblue1","deepskyblue2","deepskyblue3","deepskyblue4")

clustcol_2 <- clustcol_1[1:length(hc2_cluster)]

groups.cl2 <- mgsub(c(paste0("^",seq(1:7),"$")), clustcol_2, groups.cl) ############### define groups

colplot.clust(hc2,lab.col=groups.cl2)



plot(tsne2$Y,col=as.character(groups.cl2),pch=16)
text(tsne2$Y, labels = colnames(RPKM), cex = 0.2)

rownames(tsne2$Y) <-  as.character(as.vector(seq(1, nrow(tsne2$Y))))
write.table(RPKM,"NBT_rpkm_express.txt", col.names = T, row.names = T, quote = F, sep = "\t")
write.table(groups.cl2,"NBT_groups.cl.txt", col.names = T, row.names = T, quote = F, sep = "\t")
write.table(tsne2$Y, "NBT_tsne2_y.txt", col.names = T, row.names = T, quote = F, sep = "\t")
