setwd("~/data/PEsc/WD/")
source("../scripts/Functions/20170228_Standard_Bioinformatics_functions.R")
output <- "../WD/"

#BiocManager::install('DESeq2')
library(DESeq2)

source("../scripts/deseq2_normalization.r")

vi <- "../WD/PEsc_matrix.txt"
vii <- "../WD/DC_matrix.txt"

##################################################################

SCcountdata <- read.csv(vi, header=T, sep="\t", row.names = NULL)
DCcountdata <- read.csv(vii, header=T, sep="\t", row.names = NULL)
genes <- read.csv("../WD/genes_title.txt")

rownames(SCcountdata) <- genes[,1]
rownames(DCcountdata) <- genes[,1]

SCcountdata <- t(SCcountdata)
DCcountdata <- t(DCcountdata)

######################## deseq2 normalization ########################

SCcountdata[is.na(SCcountdata)] <- 0
SCcountdata[1,] <- SCcountdata[1,] + 1
SCnormdata <- deseq2_normalizaton(SCcountdata, 2, ncol(SCcountdata)-2)

DCcountdata[1,] <- DCcountdata[1,] + 1
DCnormdata <- deseq2_normalizaton(DCcountdata, 2, ncol(DCcountdata)-2)

### detected genes

SCdetGenes <- apply(SCnormdata, 1, function(x) sum(x > 1))
SCavDetGenes <- mean(SCdetGenes)
SDDetGenes <- sd(DetGenes)
MinDetGenes <- AvDetGenes - 1.2*SDDetGenes
MaxDetGenes <- AvDetGenes + 2*SDDetGenes
MinDetCells <- names(which(DetGenes < MinDetGenes))
MaxDetCells <- names(which(DetGenes > MaxDetGenes))
BadDetCells <- union(MinDetCells, MaxDetCells)

DCdetGenes <- apply(DCnormdata, 1, function(x) sum(x > 1))
DCavDetGenes <- mean(DCdetGenes)
SDDetGenes <- sd(DetGenes)
MinDetGenes <- AvDetGenes - 1.2*SDDetGenes
MaxDetGenes <- AvDetGenes + 2*SDDetGenes
MinDetCells <- names(which(DetGenes < MinDetGenes))
MaxDetCells <- names(which(DetGenes > MaxDetGenes))
BadDetCells <- union(MinDetCells, MaxDetCells)

Mit <- read.csv("mouse_mitochondrial_genes_list_annotation.txt", header = T, sep = "\t")
MitGenes <- match(Mit[,1], colnames(countdata))
MitSums <- rowSums(countdata[,MitGenes])
MitRatio <- MitSums/CountSums
AvMitRatio <- mean(MitRatio)
SDMitRatio <- sd(MitRatio)
MinMitRatio <- AvMitRatio - 2*SDMitRatio
MaxMitRatio <- AvMitRatio + 2*SDMitRatio
MinMitCells <- names(which(MitRatio < MinMitRatio))
MaxMitCells <- names(which(MitRatio > MaxMitRatio))
BadMitCells <- union(MinMitCells, MaxMitCells)

BadCells <- union(BadDetCells, union(BadCountCells, BadMitCells))

normdata <- t(normdata)
countdata <- t(countdata)

normdata <- normdata[,-c(which(colnames(normdata) %in% BadCells))]
countdata <- countdata[,-c(which(colnames(countdata) %in% BadCells))]

### remove non-expressed genes and bad cells

normdata <- t(normdata)
countdata <- t(countdata)

rpkm <- 1 # rpkm cutoff
cells <- 20 # expressed above the rpkm cutoff in this many cells
data <- normdata[,which(apply(normdata[c(1:(nrow(normdata)/2)),], 2, function(x) length(which(x>rpkm)))>cells)]
data <- t(data)
ExGenes <- rownames(data)

normdata <- t(normdata)
countdata <- t(countdata)

normdata <- normdata[c(which(rownames(normdata) %in% ExGenes)),]
countdata <- countdata[c(which(rownames(countdata) %in% ExGenes)),]

############# SAVE THIS FILTERED DATA ##########################
save(rpkmdata, file="SingleCells_ensemble_norm_QC.RData")
save(countdata, file="SingleCells_ensemble_count_QC.RData")
################################################################