setwd("~/data/CD63GFP/WD/")
source("../scripts/Functions/20170228_Standard_Bioinformatics_functions.R")
source("../scripts/Functions/20170228_PCA_functions.R")
source("../scripts/Functions/20170228_tSNE_Graphing_functions.R")
source("../scripts/Functions/20170228_WGCNA_SCDE.R")
#BiocManager::install('psych')
library(psych)

### load the data frame
rpkmdata1<-read.table("../WD/LI_30m_1x_2x_all_2_batch_correction_norm_express.txt", header=T)
rpkmdata2<-read.table("../WD/LI_3h_2x_all_batch_correction_norm_express.txt", header=T)
rpkmdata3<-read.table("../WD/S_30m_1x_2x_all_batch_correction_norm_express.txt", header=T)
rpkmdata4<-read.table("../WD/S_3h_2x_all_batch_correction_norm_express.txt", header=T)
rpkmdata5<-read.table("../WD/Lu_30m_1x_2x_all_batch_correction_norm_express.txt", header=T)
rpkmdata5<-rpkmdata5[,1:946]
rpkmdata6<-read.table("../WD/Lu_3h_2x_all_batch_correction_norm_express.txt", header=T)
rpkmdata<-cbind(rpkmdata1, rpkmdata2, rpkmdata3, rpkmdata4, rpkmdata5, rpkmdata6)
data <- t(rpkmdata)

#data[data==0] <- 0.0001
#data <- log2(data)

#var.rpkms<-var.genes(data)
#data <- subset(data, select = colnames(data) %in% var.rpkms)

#plotCol(nearRcolor("#000000", "rgb", dist=50))
# black, orange, lightblue, green, yellow, darkblue, red, pink
#colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
#                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#pie(rep(1, 8), col = colorBlindBlack8)
#Zissou blue,lightblue,green,yellow,orange,red,pink,purple "#3B9AB2", "#78B7C5", "#00A08A", "#EBCC2A", "#E1AF00", "#F21A00", "#D8A499", "#7294D4"

### Derive Cell Ages, DetGenes, Reads

S_HEK_30m_1x <- data[grep('^S_HEK_30m_1x', rownames(data)),]
S_HEK_30m_2x <- data[grep('^S_HEK_30m_2x', rownames(data)),]
S_HEK_3h_2x <- data[grep('^S_HEK_3h_2x', rownames(data)),]
S_CAP_30m_1x <- data[grep('^S_CAP_30m_1x', rownames(data)),]
S_CAP_30m_2x <- data[grep('^S_CAP_30m_2x', rownames(data)),]
S_CAP_3h_2x <- data[grep('^S_CAP_3h_2x', rownames(data)),]

LI_HEK_30m_1x <- data[grep('^LI_HEK_30m_1x', rownames(data)),]
LI_HEK_30m_2x <- data[grep('^LI_HEK_30m_2x', rownames(data)),]
LI_HEK_3h_2x <- data[grep('^LI_HEK_3h_2x', rownames(data)),]
LI_CAP_30m_1x <- data[grep('^LI_CAP_30m_1x', rownames(data)),]
LI_CAP_30m_2x <- data[grep('^LI_CAP_30m_2x', rownames(data)),]
LI_CAP_3h_2x <- data[grep('^LI_CAP_3h_2x', rownames(data)),]

LU_HEK_30m_1x <- data[grep('^LU_HEK_30m_1x', rownames(data)),]
LU_HEK_30m_2x <- data[grep('^LU_HEK_30m_2x', rownames(data)),]
LU_HEK_3h_2x <- data[grep('^LU_HEK_3h_2x', rownames(data)),]
LU_CAP_30m_1x <- data[grep('^LU_CAP_30m_1x', rownames(data)),]
LU_CAP_30m_2x <- data[grep('^LU_CAP_30m_2x', rownames(data)),]
LU_CAP_3h_2x <- data[grep('^LU_CAP_3h_2x', rownames(data)),]

Organ <- rep(NA, nrow(data))
names(Organ) <- rownames(data)
S_ <- data[grep('^S_', rownames(data)),]
LI_ <- data[grep('^LI_', rownames(data)),]
LU_ <- data[grep('^LU_', rownames(data)),]
Organ[which(names(Organ) %in% rownames(S_))] <- "#56B4E9"
Organ[which(names(Organ) %in% rownames(LI_))] <- "#D55E00"
Organ[which(names(Organ) %in% rownames(LU_))] <- "#F0E442"

Organ.dark <- col2rgb(Organ)
Organ.dark <- rgb2hsv(Organ.dark)
Organ.dark[3,] <- Organ.dark[3,]*0.5
Organ.dark <- hsv(Organ.dark[1,], Organ.dark[2,], Organ.dark[3,], 1)
Organ.dark[which(is.na(Organ))] <- NA

CellType <- rep(NA, nrow(data))
names(CellType) <- rownames(data)
CellType[which(names(CellType) %in% rownames(S_HEK_30m_1x))] <- "#CCCCCC"
CellType[which(names(CellType) %in% rownames(S_HEK_30m_2x))] <- "#CCCCCC"
CellType[which(names(CellType) %in% rownames(S_HEK_3h_2x))] <- "#CCCCCC"
CellType[which(names(CellType) %in% rownames(S_CAP_30m_1x))] <- "#666666"
CellType[which(names(CellType) %in% rownames(S_CAP_30m_2x))] <- "#666666"
CellType[which(names(CellType) %in% rownames(S_CAP_3h_2x))] <- "#666666"
CellType[which(names(CellType) %in% rownames(LI_HEK_30m_1x))] <- "#CCCCCC"
CellType[which(names(CellType) %in% rownames(LI_HEK_30m_2x))] <- "#CCCCCC"
CellType[which(names(CellType) %in% rownames(LI_HEK_3h_2x))] <- "#CCCCCC"
CellType[which(names(CellType) %in% rownames(LI_CAP_30m_1x))] <- "#666666"
CellType[which(names(CellType) %in% rownames(LI_CAP_30m_2x))] <- "#666666"
CellType[which(names(CellType) %in% rownames(LI_CAP_3h_2x))] <- "#666666"
CellType[which(names(CellType) %in% rownames(LU_HEK_30m_1x))] <- "#CCCCCC"
CellType[which(names(CellType) %in% rownames(LU_HEK_30m_2x))] <- "#CCCCCC"
CellType[which(names(CellType) %in% rownames(LU_HEK_3h_2x))] <- "#CCCCCC"
CellType[which(names(CellType) %in% rownames(LU_CAP_30m_1x))] <- "#666666"
CellType[which(names(CellType) %in% rownames(LU_CAP_30m_2x))] <- "#666666"
CellType[which(names(CellType) %in% rownames(LU_CAP_3h_2x))] <- "#666666"

CellType.dark <- col2rgb(CellType)
CellType.dark <- rgb2hsv(CellType.dark)
CellType.dark[3,] <- CellType.dark[3,]*0.5
CellType.dark <- hsv(CellType.dark[1,], CellType.dark[2,], CellType.dark[3,], 1)
CellType.dark[which(is.na(CellType))] <- NA


TimePoint <- rep(NA, nrow(data))
names(TimePoint) <- rownames(data)
TimePoint[which(names(TimePoint) %in% rownames(S_HEK_30m_1x))] <- "#009E73"
TimePoint[which(names(TimePoint) %in% rownames(S_HEK_30m_2x))] <- "#009E73"
TimePoint[which(names(TimePoint) %in% rownames(S_HEK_3h_2x))] <- "#CC79A7"
TimePoint[which(names(TimePoint) %in% rownames(S_CAP_30m_1x))] <- "#009E73"
TimePoint[which(names(TimePoint) %in% rownames(S_CAP_30m_2x))] <- "#009E73"
TimePoint[which(names(TimePoint) %in% rownames(S_CAP_3h_2x))] <- "#CC79A7"
TimePoint[which(names(TimePoint) %in% rownames(LI_HEK_30m_1x))] <- "#009E73"
TimePoint[which(names(TimePoint) %in% rownames(LI_HEK_30m_2x))] <- "#009E73"
TimePoint[which(names(TimePoint) %in% rownames(LI_HEK_3h_2x))] <- "#CC79A7"
TimePoint[which(names(TimePoint) %in% rownames(LI_CAP_30m_1x))] <- "#009E73"
TimePoint[which(names(TimePoint) %in% rownames(LI_CAP_30m_2x))] <- "#009E73"
TimePoint[which(names(TimePoint) %in% rownames(LI_CAP_3h_2x))] <- "#CC79A7"
TimePoint[which(names(TimePoint) %in% rownames(LU_HEK_30m_1x))] <- "#009E73"
TimePoint[which(names(TimePoint) %in% rownames(LU_HEK_30m_2x))] <- "#009E73"
TimePoint[which(names(TimePoint) %in% rownames(LU_HEK_3h_2x))] <- "#CC79A7"
TimePoint[which(names(TimePoint) %in% rownames(LU_CAP_30m_1x))] <- "#009E73"
TimePoint[which(names(TimePoint) %in% rownames(LU_CAP_30m_2x))] <- "#009E73"
TimePoint[which(names(TimePoint) %in% rownames(LU_CAP_3h_2x))] <- "#CC79A7"

TimePoint.dark <- col2rgb(TimePoint)
TimePoint.dark <- rgb2hsv(TimePoint.dark)
TimePoint.dark[3,] <- TimePoint.dark[3,]*0.5
TimePoint.dark <- hsv(TimePoint.dark[1,], TimePoint.dark[2,], TimePoint.dark[3,], 1)
TimePoint.dark[which(is.na(TimePoint))] <- NA


Dose <- rep(NA, nrow(data))
names(Dose) <- rownames(data)
Dose[which(names(Dose) %in% rownames(S_HEK_30m_1x))] <- "#696969"
Dose[which(names(Dose) %in% rownames(S_HEK_30m_2x))] <- "#000000"
Dose[which(names(Dose) %in% rownames(S_HEK_3h_2x))] <- "#000000"
Dose[which(names(Dose) %in% rownames(S_CAP_30m_1x))] <- "#696969"
Dose[which(names(Dose) %in% rownames(S_CAP_30m_2x))] <- "#000000"
Dose[which(names(Dose) %in% rownames(S_CAP_3h_2x))] <- "#000000"
Dose[which(names(Dose) %in% rownames(LI_HEK_30m_1x))] <- "#696969"
Dose[which(names(Dose) %in% rownames(LI_HEK_30m_2x))] <- "#000000"
Dose[which(names(Dose) %in% rownames(LI_HEK_3h_2x))] <- "#000000"
Dose[which(names(Dose) %in% rownames(LI_CAP_30m_1x))] <- "#696969"
Dose[which(names(Dose) %in% rownames(LI_CAP_30m_2x))] <- "#000000"
Dose[which(names(Dose) %in% rownames(LI_CAP_3h_2x))] <- "#000000"
Dose[which(names(Dose) %in% rownames(LU_HEK_30m_1x))] <- "#696969"
Dose[which(names(Dose) %in% rownames(LU_HEK_30m_2x))] <- "#000000"
Dose[which(names(Dose) %in% rownames(LU_HEK_3h_2x))] <- "#000000"
Dose[which(names(Dose) %in% rownames(LU_CAP_30m_1x))] <- "#696969"
Dose[which(names(Dose) %in% rownames(LU_CAP_30m_2x))] <- "#000000"
Dose[which(names(Dose) %in% rownames(LU_CAP_3h_2x))] <- "#000000"

Dose.dark <- col2rgb(Dose)
Dose.dark <- rgb2hsv(Dose.dark)
Dose.dark[3,] <- Dose.dark[3,]*0.5
Dose.dark <- hsv(Dose.dark[1,], Dose.dark[2,], Dose.dark[3,], 1)
Dose.dark[which(is.na(Dose))] <- NA



hc1cells<-hclust(as.dist(1-cor(t(data), method="spearman")), method="ward")
cellgroups<-cutree(hc1cells,k=15)
clusters <- rep(NA, nrow(data))
names(clusters) <- rownames(data)
clust1 <- which(cellgroups == 1)
clust2 <- which(cellgroups == 2)
clust3 <- which(cellgroups == 3)
clust4 <- which(cellgroups == 4)
clust5 <- which(cellgroups == 5)
clust6 <- which(cellgroups == 6)
clust7 <- which(cellgroups == 7)
clust8 <- which(cellgroups == 8)
clust9 <- which(cellgroups == 9)
clust10 <- which(cellgroups == 10)
clust11 <- which(cellgroups == 11)
clust12 <- which(cellgroups == 12)
clust13 <- which(cellgroups == 13)
clust14 <- which(cellgroups == 14)
clust15 <- which(cellgroups == 15)
clusters[which(names(clusters) %in% names(clust1))] <- "blue"
clusters[which(names(clusters) %in% names(clust2))] <- "red"
clusters[which(names(clusters) %in% names(clust3))] <- "yellow"
clusters[which(names(clusters) %in% names(clust4))] <- "green"
clusters[which(names(clusters) %in% names(clust5))] <- "orange"
clusters[which(names(clusters) %in% names(clust6))] <- "purple"
clusters[which(names(clusters) %in% names(clust7))] <- "black"
clusters[which(names(clusters) %in% names(clust8))] <- "darkblue"
clusters[which(names(clusters) %in% names(clust9))] <- "darkred"
clusters[which(names(clusters) %in% names(clust10))] <- "yellow4"
clusters[which(names(clusters) %in% names(clust11))] <- "darkgreen"
clusters[which(names(clusters) %in% names(clust12))] <- "darkorange"
clusters[which(names(clusters) %in% names(clust13))] <- "purple4"
clusters[which(names(clusters) %in% names(clust14))] <- "grey66"
clusters[which(names(clusters) %in% names(clust15))] <- "grey33"

clusters.dark <- col2rgb(clusters)
clusters.dark <- rgb2hsv(clusters.dark)
clusters.dark[3,] <- clusters.dark[3,]*0.5
clusters.dark <- hsv(clusters.dark[1,], clusters.dark[2,], clusters.dark[3,], 1)
clusters.dark[which(is.na(clusters))] <- NA

DetGenes <- apply(data, 1, function(x) sum(x > log2(2)))
NoDetGenes <- DetGenes

DetGenes <- cRamp.orig(DetGenes, colors = c("green","yellow","orange","red","darkred"))
DetGenes.dark <- col2rgb(DetGenes)
DetGenes.dark <- rgb2hsv(DetGenes.dark)
DetGenes.dark[3,] <- DetGenes.dark[3,]*0.5
DetGenes.dark <- hsv(DetGenes.dark[1,], DetGenes.dark[2,], DetGenes.dark[3,], 1)
DetGenes.dark[which(is.na(DetGenes))] <- NA

ReadNumbers <- colSums(rpkmdata)

ReadCount <- cRamp.orig(ReadNumbers, colors = c("cyan", "green","yellow","orange","red","darkred"))
ReadCount.dark <- col2rgb(ReadCount)
ReadCount.dark <- rgb2hsv(ReadCount.dark)
ReadCount.dark[3,] <- ReadCount.dark[3,]*0.5
ReadCount.dark <- hsv(ReadCount.dark[1,], ReadCount.dark[2,], ReadCount.dark[3,], 1)
ReadCount.dark[which(is.na(ReadCount))] <- NA


RFliver30m <- read.table("../WD/LI_30m_1x_all_rf_prediction.txt", header=T)
RFliver3h <- read.table("../WD/LI_3h_2x_all_rf_prediction.txt", header=T)
RFspleen30m <- read.table("../WD/S_30m_1x_all_rf_prediction.txt", header=T)
RFspleen3h <- read.table("../WD/S_3h_2x_all_rf_prediction.txt", header=T)
RFlung30m <- read.table("../WD/LU_30m_1x_all_rf_prediction.txt", header=T)
RFlung3h <- read.table("../WD/LU_3h_2x_all_rf_prediction.txt", header=T)
RF <- rbind(RFliver30m, RFliver3h, RFspleen30m, RFspleen3h, RFlung30m, RFlung3h)
rownames(RF) <- RF$X0

B_cell <- RF[which(RF[,3] == "B_cell"),]
Hepatocyte <- RF[which(RF[,3] == "hepatocyte"),]
Kupffer <- RF[which(RF[,3] == "Kupffer_cell"),]
LungEndo <- RF[which(RF[,3] == "lung_endothelial_cell"),]
T_cell <- RF[which(RF[,3] == "T_cell"),]
Myeloid <- RF[which(RF[,3] == "myeloid_cell"),]
NK_cell <- RF[which(RF[,3] == "natural_killer_cell"),]
Macrophage <- RF[which(RF[,3] == "macrophage"),]
ClasMono <- RF[which(RF[,3] == "classical_monocyte"),]
Epithelial <- RF[which(RF[,3] == "epithelial_cell_of_lung"),]
Leukocyte <- RF[which(RF[,3] == "leukocyte"),]
Monocyte <- RF[which(RF[,3] == "monocyte"),]
Stromal <- RF[which(RF[,3] == "stromal_cell"),]

RandomForest <- rep(NA, nrow(data))
names(RandomForest) <- rownames(data)
RandomForest[which(names(RandomForest) %in% rownames(B_cell))] <- "#78B7C5"
RandomForest[which(names(RandomForest) %in% rownames(Hepatocyte))] <- "#00A08A"
RandomForest[which(names(RandomForest) %in% rownames(Kupffer))] <- "#E1AF00"
RandomForest[which(names(RandomForest) %in% rownames(LungEndo))] <- "#F21A00"
RandomForest[which(names(RandomForest) %in% rownames(T_cell))] <- "#7294D4"
RandomForest[which(names(RandomForest) %in% rownames(Myeloid))] <- "#E1AF00"
RandomForest[which(names(RandomForest) %in% rownames(NK_cell))] <- "#D8A499"
RandomForest[which(names(RandomForest) %in% rownames(Macrophage))] <- "#E1AF00"
RandomForest[which(names(RandomForest) %in% rownames(ClasMono))] <- "#EBCC2A"
RandomForest[which(names(RandomForest) %in% rownames(Epithelial))] <- "#00A08A"
RandomForest[which(names(RandomForest) %in% rownames(Leukocyte))] <- "#3B9AB2"
RandomForest[which(names(RandomForest) %in% rownames(Monocyte))] <- "#BABABA"
RandomForest[which(names(RandomForest) %in% rownames(Stromal))] <- "#000000"

RandomForest.dark <- col2rgb(RandomForest)
RandomForest.dark <- rgb2hsv(RandomForest.dark)
RandomForest.dark[3,] <- RandomForest.dark[3,]*0.5
RandomForest.dark <- hsv(RandomForest.dark[1,], RandomForest.dark[2,], RandomForest.dark[3,], 1)
RandomForest.dark[which(is.na(RandomForest))] <- NA


S_HEK_30m_1x_A <- data[grep('^S_HEK_30m_1x_A', rownames(data)),]
S_HEK_30m_2x_A <- data[grep('^S_HEK_30m_2x_A', rownames(data)),]
S_HEK_3h_2x_A <- data[grep('^S_HEK_3h_2x_A', rownames(data)),]
S_CAP_30m_1x_A <- data[grep('^S_CAP_30m_1x_A', rownames(data)),]
S_CAP_30m_2x_A <- data[grep('^S_CAP_30m_2x_A', rownames(data)),]
S_CAP_3h_2x_A <- data[grep('^S_CAP_3h_2x_A', rownames(data)),]
LI_HEK_30m_1x_A <- data[grep('^LI_HEK_30m_1x_A', rownames(data)),]
LI_HEK_30m_2x_A <- data[grep('^LI_HEK_30m_2x_A', rownames(data)),]
LI_HEK_3h_2x_A <- data[grep('^LI_HEK_3h_2x_A', rownames(data)),]
LI_CAP_30m_1x_A <- data[grep('^LI_CAP_30m_1x_A', rownames(data)),]
LI_CAP_30m_2x_A <- data[grep('^LI_CAP_30m_2x_A', rownames(data)),]
LI_CAP_3h_2x_A <- data[grep('^LI_CAP_3h_2x_A', rownames(data)),]
LU_HEK_30m_1x_A <- data[grep('^LU_HEK_30m_1x_A', rownames(data)),]
LU_HEK_30m_2x_A <- data[grep('^LU_HEK_30m_2x_A', rownames(data)),]
LU_HEK_3h_2x_A <- data[grep('^LU_HEK_3h_2x_A', rownames(data)),]
LU_CAP_30m_1x_A <- data[grep('^LU_CAP_30m_1x_A', rownames(data)),]
LU_CAP_30m_2x_A <- data[grep('^LU_CAP_30m_2x_A', rownames(data)),]
LU_CAP_3h_2x_A <- data[grep('^LU_CAP_3h_2x_A', rownames(data)),]
S_HEK_30m_1x_B <- data[grep('^S_HEK_30m_1x_B', rownames(data)),]
S_HEK_30m_2x_B <- data[grep('^S_HEK_30m_2x_B', rownames(data)),]
S_HEK_3h_2x_B <- data[grep('^S_HEK_3h_2x_B', rownames(data)),]
S_CAP_30m_1x_B <- data[grep('^S_CAP_30m_1x_B', rownames(data)),]
S_CAP_30m_2x_B <- data[grep('^S_CAP_30m_2x_B', rownames(data)),]
S_CAP_3h_2x_B <- data[grep('^S_CAP_3h_2x_B', rownames(data)),]
LI_HEK_30m_1x_B <- data[grep('^LI_HEK_30m_1x_B', rownames(data)),]
LI_HEK_30m_2x_B <- data[grep('^LI_HEK_30m_2x_B', rownames(data)),]
LI_HEK_3h_2x_B <- data[grep('^LI_HEK_3h_2x_B', rownames(data)),]
LI_CAP_30m_1x_B <- data[grep('^LI_CAP_30m_1x_B', rownames(data)),]
LI_CAP_30m_2x_B <- data[grep('^LI_CAP_30m_2x_B', rownames(data)),]
LI_CAP_3h_2x_B <- data[grep('^LI_CAP_3h_2x_B', rownames(data)),]
LU_HEK_30m_1x_B <- data[grep('^LU_HEK_30m_1x_B', rownames(data)),]
LU_HEK_30m_2x_B <- data[grep('^LU_HEK_30m_2x_B', rownames(data)),]
LU_HEK_3h_2x_B <- data[grep('^LU_HEK_3h_2x_B', rownames(data)),]
LU_CAP_30m_1x_B <- data[grep('^LU_CAP_30m_1x_B', rownames(data)),]
LU_CAP_30m_2x_B <- data[grep('^LU_CAP_30m_2x_B', rownames(data)),]
LU_CAP_3h_2x_B <- data[grep('^LU_CAP_3h_2x_B', rownames(data)),]

Batch <- rep(NA, nrow(data))
names(Batch) <- rownames(data)
Batch[which(names(Batch) %in% rownames(S_HEK_30m_1x_A))] <- "#9A8822"
Batch[which(names(Batch) %in% rownames(S_HEK_30m_2x_A))] <- "#9A8822"
Batch[which(names(Batch) %in% rownames(S_HEK_3h_2x_A))] <- "#9A8822"
Batch[which(names(Batch) %in% rownames(S_CAP_30m_1x_A))] <- "#9A8822"
Batch[which(names(Batch) %in% rownames(S_CAP_30m_2x_A))] <- "#9A8822"
Batch[which(names(Batch) %in% rownames(S_CAP_3h_2x_A))] <- "#9A8822"
Batch[which(names(Batch) %in% rownames(LI_HEK_30m_1x_A))] <- "#9A8822"
Batch[which(names(Batch) %in% rownames(LI_HEK_30m_2x_A))] <- "#9A8822"
Batch[which(names(Batch) %in% rownames(LI_HEK_3h_2x_A))] <- "#9A8822"
Batch[which(names(Batch) %in% rownames(LI_CAP_30m_1x_A))] <- "#9A8822"
Batch[which(names(Batch) %in% rownames(LI_CAP_30m_2x_A))] <- "#9A8822"
Batch[which(names(Batch) %in% rownames(LI_CAP_3h_2x_A))] <- "#9A8822"
Batch[which(names(Batch) %in% rownames(LU_HEK_30m_1x_A))] <- "#9A8822"
Batch[which(names(Batch) %in% rownames(LU_HEK_30m_2x_A))] <- "#9A8822"
Batch[which(names(Batch) %in% rownames(LU_HEK_3h_2x_A))] <- "#9A8822"
Batch[which(names(Batch) %in% rownames(LU_CAP_30m_1x_A))] <- "#9A8822"
Batch[which(names(Batch) %in% rownames(LU_CAP_30m_2x_A))] <- "#9A8822"
Batch[which(names(Batch) %in% rownames(LU_CAP_3h_2x_A))] <- "#9A8822"
Batch[which(names(Batch) %in% rownames(S_HEK_30m_1x_B))] <- "#74A089"
Batch[which(names(Batch) %in% rownames(S_HEK_30m_2x_B))] <- "#74A089"
Batch[which(names(Batch) %in% rownames(S_HEK_3h_2x_B))] <- "#74A089"
Batch[which(names(Batch) %in% rownames(S_CAP_30m_1x_B))] <- "#74A089"
Batch[which(names(Batch) %in% rownames(S_CAP_30m_2x_B))] <- "#74A089"
Batch[which(names(Batch) %in% rownames(S_CAP_3h_2x_B))] <- "#74A089"
Batch[which(names(Batch) %in% rownames(LI_HEK_30m_1x_B))] <- "#74A089"
Batch[which(names(Batch) %in% rownames(LI_HEK_30m_2x_B))] <- "#74A089"
Batch[which(names(Batch) %in% rownames(LI_HEK_3h_2x_B))] <- "#74A089"
Batch[which(names(Batch) %in% rownames(LI_CAP_30m_1x_B))] <- "#74A089"
Batch[which(names(Batch) %in% rownames(LI_CAP_30m_2x_B))] <- "#74A089"
Batch[which(names(Batch) %in% rownames(LI_CAP_3h_2x_B))] <- "#74A089"
Batch[which(names(Batch) %in% rownames(LU_HEK_30m_1x_B))] <- "#74A089"
Batch[which(names(Batch) %in% rownames(LU_HEK_30m_2x_B))] <- "#74A089"
Batch[which(names(Batch) %in% rownames(LU_HEK_3h_2x_B))] <- "#74A089"
Batch[which(names(Batch) %in% rownames(LU_CAP_30m_1x_B))] <- "#74A089"
Batch[which(names(Batch) %in% rownames(LU_CAP_30m_2x_B))] <- "#74A089"
Batch[which(names(Batch) %in% rownames(LU_CAP_3h_2x_B))] <- "#74A089"

Batch.dark <- col2rgb(Batch)
Batch.dark <- rgb2hsv(Batch.dark)
Batch.dark[3,] <- Batch.dark[3,]*0.5
Batch.dark <- hsv(Batch.dark[1,], Batch.dark[2,], Batch.dark[3,], 1)
Batch.dark[which(is.na(Batch))] <- NA

#############################################################################

CellFeat <- as.data.frame(matrix(NA, nrow(data), 20))
CellFeat[,1] <- as.character(Organ)
CellFeat[,2] <- as.character(Organ.dark)
CellFeat[,3] <- as.character(CellType)
CellFeat[,4] <- as.character(CellType.dark)
CellFeat[,5] <- as.character(TimePoint)
CellFeat[,6] <- as.character(TimePoint.dark)
CellFeat[,7] <- as.character(Dose)
CellFeat[,8] <- as.character(Dose.dark)
CellFeat[,9] <- as.numeric(NoDetGenes)
CellFeat[,10] <- as.character(DetGenes)
CellFeat[,11] <- as.character(DetGenes.dark)
CellFeat[,12] <- as.numeric(ReadNumbers)
CellFeat[,13] <- as.character(ReadCount)
CellFeat[,14] <- as.character(ReadCount.dark)
CellFeat[,15] <- as.character(clusters)
CellFeat[,16] <- as.character(clusters.dark)
CellFeat[,17] <- as.character(RandomForest)
CellFeat[,18] <- as.character(RandomForest.dark)
CellFeat[,19] <- as.character(Batch)
CellFeat[,20] <- as.character(Batch.dark)
rownames(CellFeat) <- rownames(data)

colnames(CellFeat) <- c("Organ", "Organ.dark", "CellType", "CellType.dark", "TimePoint", "Timepoint.dark","Dose", "Dose.dark",
                          "NoDetGenes", "DetGenes", "DetGenes.dark","ReadNumbers", "ReadCount", "ReadCount.dark",
                          "clusters", "clusters.dark", "RandomForest", "RandomForest.dark", "Batch", "Batch.dark")

save(CellFeat, file = "QCCells_Features.RData")
