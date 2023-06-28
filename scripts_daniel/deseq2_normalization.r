#!/bin/R

deseq2_normalizaton <- function(countdata, control, case)
{

	library(DESeq2)
	################# main calculation #################

	# setup the sample condition in order to perform wald test

	sampleCondition<-c(replicate(control,"treated"),replicate(case,"untreated"))
	sampleCondition<-as.factor(sampleCondition)


	coldata_u <- data.frame(sampleCondition)
	rownames(coldata_u) <- colnames(countdata)

	# Loading the samples
	dds <- DESeqDataSetFromMatrix(countData = countdata , colData = coldata_u, design = ~ sampleCondition)

	# perform the calculation
	dds <- estimateSizeFactors(dds)

    # size factor : dds$sizeFactor

	# baseMean is average of normalized counts

	# The shrinkage is greater for the log2 fold change estimates from genes with low counts and high dispersion
	normCounts <- counts(dds, normalized=TRUE)

	return (normCounts)
}

