setwd("~/data/CD63GFP/WD/")
source("../scripts/Functions/20170228_Standard_Bioinformatics_functions.R")
output <- "../WD/"

#BiocManager::install('tximport')
#BiocManager::install('AnnotationHub')
BiocManager::install('ensembldb')
library(ensembldb)
library(AnnotationHub)
library(AcidGenomes)
library(tximport)

viia <- "../WD/mouse_total_raw_count_1.txt"
  vii <- "../WD/raw_count_mouse.txt"

##################################################################

countdata <- read.csv(vii, header=T, sep="\t", row.names = NULL)
  countdata1 <- read.csv(viia, header=T, sep="\t", row.names = NULL)
  countdata1 <- countdata1[,290:577]
countdata <- cbind(countdata, countdata1)

rownames(countdata) <- countdata[,1]
countdata <- countdata[,-1]

samples <- t(read.csv("../WD/samples.csv", header = T))
colnames(rpkmdata) <- make.unique(samples[2,])
colnames(countdata) <- make.unique(samples[2,])

countdata <- t(countdata)

tx2gene <- makeTx2GeneFromEnsembl("Mus musculus")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
rpkmdata <- 

S_countdata <- countdata[grep('^S_', rownames(countdata)),]
LI_countdata <- countdata[grep('^LI_', rownames(countdata)),]
LU_countdata <- countdata[grep('^LU_', rownames(countdata)),]
countdata <- rbind(S_countdata, LI_countdata, LU_countdata)

### detected genes

CountSums <- rowSums(countdata)
MinCountCells <- names(which(CountSums < 1000))
MaxCountCells <- names(which(CountSums > 1000000))
BadCountCells <- union(MinCountCells, MaxCountCells)

DetGenes <- apply(rpkmdata, 1, function(x) sum(x > 1))
AvDetGenes <- mean(DetGenes)
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

rpkmdata <- t(rpkmdata)
countdata <- t(countdata)

rpkmdata <- rpkmdata[,-c(which(colnames(rpkmdata) %in% BadCells))]
countdata <- countdata[,-c(which(colnames(countdata) %in% BadCells))]

### remove non-expressed genes and bad cells

rpkmdata <- t(rpkmdata)
countdata <- t(countdata)

rpkm <- 1 # rpkm cutoff
cells <- 20 # expressed above the rpkm cutoff in this many cells
data <- rpkmdata[,which(apply(rpkmdata[c(1:(nrow(rpkmdata)/2)),], 2, function(x) length(which(x>rpkm)))>cells)]
data <- t(data)
ExGenes <- rownames(data)

rpkmdata <- t(rpkmdata)
countdata <- t(countdata)

rpkmdata <- rpkmdata[c(which(rownames(rpkmdata) %in% ExGenes)),]
countdata <- countdata[c(which(rownames(countdata) %in% ExGenes)),]

############# SAVE THIS FILTERED DATA ##########################
save(rpkmdata, file="SingleCells_ensemble_rpkm_QC.RData")
save(countdata, file="SingleCells_ensemble_count_QC.RData")
################################################################