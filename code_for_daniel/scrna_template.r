#!/bin/R




library(Rtsne)
library(org.Mm.eg.db)
library(DESeq2)
library(edgeR)
library(gplots)
library(plotrix)



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


mit <- read.csv("mouse_mitochondrial_genes_list_annotation.txt", header = T, sep = "\t")
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


cc <- colorRampPalette(c("yellow","red","black"))
col.nDet <- convert.to.color(nD.tx,cc)


################## plot tSNE #####################

tsne2 <- Rtsne(t(log2(NORM+1)),  perplexity = 25) ############# need changes ##############

par(mfrow=c(2,2),mar=c(2,2,2,2),oma=c(1,1,1,1))

plot(tsne2$Y,col=col.nDet$cols,pch=16)

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


################## plot single gene tSNE #####################
#RPKM[46518:46610,] -> ercc_gfp

keys <- as.character(rownames(RPKM))
annotations <- select(org.Hs.eg.db, keys,  c("SYMBOL","GENENAME"), keytype="ENSEMBL")

match_index <- match(annotations[,1], keys)
RPKM <- RPKM[match_index,]

rownames(RPKM) <- annotations[,2]

candidate_genes <- c("Cd9", "Cd200r1", "Acta2", "Cldn5", "Pecam1", "Ptprc" ,"Podxl", "Pdgfra" ,"Pdgfrb", "Cdh5","Col1a1","Cd68" ,"Actb", "Gapdh", "Pecam1", "Pdgfrb", "Ptprc","Adgre1","Pecam1","Cd34","Cd3e", "Cd3d", "Cd3g","Cd19", "Cd79a", "Cd79b", "Cd22", "Ly6d", "Siglech", "Cd68", "Lyz2", "Cd14", "Cd55", "Cd5l", "Clec4f", "Cdh5", "Vsig4", "Clec1b", "Slc40a1", "Alb", "Rplp0", "Actb", "Adipoq", "Alb", "Cd3g", "Eng", "Igfbp7", "Cd19", "Cd68", "Cd3e", "Cd3d", "Cd3g","Cd19", "Cd79a", "Cd79b", "Cd22", "Ly6d", "Siglech", "Cd68", "Lyz2", "Cd14", "Cd55", "Cd5l", "Clec4f", "Cdh5", "Vsig4", "Clec1b", "Slc40a1", "Alb", "Rplp0", "Actb", "Adipoq", "Alpl", "Icam1", "Kit", "Nrp2", "Vegfc", "Vwf", "Pecam1", "Cdh5", "Cldn5",  "Kdr", "Tie1", "Flt4", "Prox1", "Bmp3", "Bmper", "Chad",  "Cxcl14", "Mkx", "Pdlim3", "Tmeff2", "Col1a1", "Col1a2",  "Dcn", "Lum", "Pdgfra", "Cd68", "Ptprc", "Acta2", "Pdgfrb", "Abcc9", "Kcnj8", "Cspg4", "Rgs5", "Vtn", "Actb", "Gapdh","Siglech", "Igfbp7", "Cd3g", "Cd79a", "Alb", "Lyz2","Plvap","Cldn5", "Slc1a1",  "Ptprc","Adgre1","Pecam1","Cd34","Cd3e", "Cd3d", "Cd3g","Cd19", "Cd79a", "Cd79b", "Cd22", "Ly6d", "Siglech", "Cd68", "Lyz2", "Cd14", "Cd55", "Cd5l", "Clec4f", "Cdh5", "Vsig4", "Clec1b", "Slc40a1", "Alb", "Rplp0", "Actb", "Adipoq", "Alb", "Cd3g", "Eng", "Igfbp7", "Cd19", "Cd68", "Cd3e", "Cd3d", "Cd3g","Cd19", "Cd79a", "Cd79b", "Cd22", "Ly6d", "Siglech", "Cd68", "Lyz2", "Cd14", "Cd55", "Cd5l", "Clec4f", "Cdh5", "Vsig4", "Clec1b", "Slc40a1", "Alb", "Rplp0", "Actb", "Adipoq")

#candidate_genes <- toupper(c("Acta2", "Cldn5", "Pecam1", "Ptprc" ,"Podxl", "Pdgfra" ,"Pdgfrb", "Cdh5","Col1a1","Cd68" ,"Actb", "Gapdh", "Pecam1", "Pdgfrb", "Ptprc","Adgre1","Pecam1","Cd34","Cd3e", "Cd3d", "Cd3g","Cd19", "Cd79a", "Cd79b", "Cd22", "Ly6d", "Cd68", "Cd14", "Cd55", "Cd5l", "Clec4f", "Cdh5", "Vsig4", "Clec1b", "Slc40a1", "Alb", "Rplp0", "Actb", "Adipoq", "Alb", "Cd3g", "Eng", "Igfbp7", "Cd19", "Cd68", "Cd3e", "Cd3d", "Cd3g","Cd19", "Cd79a", "Cd79b", "Cd22",   "Cd14", "Cd55", "Cd5l", "Clec4f", "Cdh5", "Vsig4", "Clec1b", "Slc40a1", "Alb", "Rplp0", "Actb", "Adipoq", "Alpl", "Icam1", "Kit", "Nrp2", "Vegfc", "Vwf", "Pecam1", "Cdh5", "Cldn5",  "Kdr", "Tie1", "Flt4", "Prox1", "Bmp3", "Bmper", "Chad",  "Cxcl14", "Mkx", "Pdlim3", "Tmeff2", "Col1a1", "Col1a2",  "Dcn", "Lum", "Pdgfra", "Cd68", "Ptprc", "Acta2", "Pdgfrb", "Abcc9", "Kcnj8", "Cspg4", "Rgs5", "Vtn", "Actb", "Gapdh", "Igfbp7", "Cd3g", "Cd79a", "Alb", "Plvap","Cldn5", "Slc1a1",  "Ptprc","Adgre1","Pecam1","Cd34","Cd3e", "Cd3d", "Cd3g","Cd19", "Cd79a", "Cd79b", "Cd22",  "Clec4f", "Cdh5", "Vsig4", "Clec1b", "Slc40a1", "Alb", "Rplp0", "Actb", "Adipoq", "Alb", "Cd3g", "Eng", "Igfbp7", "Cd19", "Cd68", "Cd3e", "Cd3d", "Cd3g","Cd19", "Cd79a", "Cd79b", "Cd22",  "Cd68",  "Cd14", "Cd55", "Cd5l", "Clec4f", "Cdh5", "Vsig4", "Clec1b", "Slc40a1", "Alb", "Rplp0", "Actb", "Adipoq"))

for (i in 1:length(candidate_genes))
{
    
    filename = paste(candidate_genes[i], '.pdf', sep = "") ###################### need change #########################
    
    pdf (file= filename)
    bla.cl <- convert.to.color(RPKM[candidate_genes[i],],cc)
    
    print (i)
    print (candidate_genes[i])
    plot(tsne2$Y,col=bla.cl$cols,pch=16, main = paste(candidate_genes[i],"")) ###################### need change #########################
    
    #dev.off()
    
}


