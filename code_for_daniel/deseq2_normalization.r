#!/bin/R








deseq2_normalizaton <- function(express_matrix, control, case)
{

	library(DESeq2)
	################# main calculation #################

	# setup the sample condition in order to perform wald test

	sampleCondition<-c(replicate(control,"treated"),replicate(case,"untreated"))


	coldata_u <- data.frame(sampleCondition)
	rownames(coldata_u) <- colnames(express_matrix)

	# Loading the samples
	dds <- DESeqDataSetFromMatrix(countData = express_matrix , colData = coldata_u, design = ~ sampleCondition)

	# perform the calculation
	dds <- estimateSizeFactors(dds)

    # size factor : dds$sizeFactor

	# baseMean is average of normalized counts

	# The shrinkage is greater for the log2 fold change estimates from genes with low counts and high dispersion
	normCounts <- counts(dds, normalized=TRUE)

	return (normCounts)
}

