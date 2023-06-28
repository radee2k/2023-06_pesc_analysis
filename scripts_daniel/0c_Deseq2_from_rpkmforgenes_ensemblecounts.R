setwd("~/Desktop/bioinformatics/R/CD63GFP/WD/")
output <- "../WD/"
args<-commandArgs(TRUE)
#source("https://bioconductor.org/biocLite.R")
#biocLite('DESeq2')
library(DESeq2)

#Data<-read.csv("../WD/BulksCounts.csv", header = T, sep = ";")
#Data<-Data[-(which(duplicated(Data[,1]))),]
#Data<-t(Data)
#colnames(Data)<-Data[1,]
#Data<-Data[-1,]
#row.names<-rownames(Data)
#Data <- apply(Data, 2, function(x) as.numeric(as.character(x)))
#rownames(Data)<-row.names

Data<-read.csv("../WD/ENCODELiverSpleenTPMs.csv", header = T, sep = ";")
Data<-t(Data)
colnames(Data)<-Data[1,]
Data<-Data[-1,]
rownames(Data) <- (c("Liver1","Liver2,","Spleen1","Spleen2"))

ControlSample<-"Liver"
TestSample<-"Spleen"

# subset to consider on the samples you want to.
Vector<-vector(mode = "logical", length = nrow(Data))
for (i in 1:nrow(Data)){
  if(grepl(paste(ControlSample, "*", sep=""), rownames(Data)[i]))
    Vector[i]<-TRUE
  if(grepl(paste(TestSample, "*", sep=""), rownames(Data)[i]))
    Vector[i]<-TRUE
}
Data<-subset(Data, Vector)

Groups<-as.matrix(paste(rownames(Data)))
Groups[which(grepl(paste(ControlSample, "*", sep=""), Groups[,1])),1]<-as.character(ControlSample)
Groups[which(grepl(paste(TestSample, "*", sep=""), Groups[,1])),1]<-as.character(TestSample)

# make the data structure
countData<-as.matrix(t(Data))
colData<-cbind(as.data.frame(colnames(countData)),Groups)
colnames(colData)<-c("Sample","Group")
countData[c(1:4),] <- apply(countData[c(1:4),],2, function(x) as.numeric(as.character(x)))
dds<-DESeqDataSetFromMatrix(countData, colData, design=~Group)
colData(dds)$Group <- factor(colData(dds)$Group, levels=c(ControlSample, TestSample))

# Run Deseq2
dds <- DESeq(dds)
res <- results(dds)

# plot the results
png("DESeq2.png", width=500, height=500, units='px')
plotMA(dds,ylim=c(-4,4),main=sprintf("%s_vs_%s", ControlSample, TestSample))
dev.off()

# select only the significant genes
save(res, file="DESeq2.RData")
res.sig.unique<-subset(res, res[,5] < 0.05)
save(res.sig.unique, file="DESeq2sig.RData")
write.csv(res.sig.unique, "DESeq-Results_%s_vs_%s_siggenes.csv")
