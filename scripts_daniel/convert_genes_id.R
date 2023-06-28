setwd("~/data/CD63GFP/WD/")
output <- "../WD/"

#BiocManager::install('remotes')
#iocManager::install('grimbough/biomaRt')
#BiocManager::install('BiocFileCache')
library(biomaRt)

#a <- "../WD/Lung50TabulaGenes.txt"
#data <- as.data.frame(read.csv(a, header=F))

#Liver <- read.csv("../WD/LiverMarkerGenes50cellsFCmin2.csv", header = T)
#Lung <- read.csv("../WD/LungMarkerGenes50cellsFCmin2.csv", header = T)
Spleen <- read.csv("../WD/SpleenMarkerGenes50cellsFCmin2.csv", header = T)

require("biomaRt")

mart = useMart("ensembl", dataset="mmusculus_gene_ensembl")

#keys <- data[,1]
keys <- Spleen$Tcell[1:100]

annotations <- getBM(
mart=mart,
attributes=c("external_gene_name","description", "ensembl_gene_id"),
filter="external_gene_name",
values=keys,
uniqueRows=FALSE)

write.csv(annotations, file="../WD/SpleenTcell100GenesEnsembl.txt")
