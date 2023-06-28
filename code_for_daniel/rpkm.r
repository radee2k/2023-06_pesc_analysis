#!/bin/R




 rpkm <- function(genes_length_filename, input_matrix, ...)
 {

		library(edgeR) 

		gene_length <- read.csv(genes_length_filename, header = T, sep = "\t")

		d <- DGEList(counts = input_matrix)

		#cpm(d)
		#cpm(d,log=TRUE)

		d$genes$Length <- gene_length[,1]

		rpkm_express <- rpkm(d)

 		return (rpkm_express)
 
 }