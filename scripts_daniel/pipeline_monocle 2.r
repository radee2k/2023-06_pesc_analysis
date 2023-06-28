#!/bin/R



library(Rtsne)
library(org.Mm.eg.db)
library(DESeq2)
library(edgeR)
library(gplots)
library(plotrix)
library(monocle)



source("~/R_lib/scrna_lib.r")

x <- read.csv("CD63GFPcounts.csv", header = T)
express_matrix <- x[,2:ncol(x)]
rownames(express_matrix) <- x[,1]

keys <- as.character(rownames(express_matrix))
annotations <- select(org.Mm.eg.db, keys,  c("SYMBOL","GENENAME"), keytype="ENSEMBL")


gene_length <- data.frame(read.csv("genes_length.txt", header = T, sep = "\t")[1:46609,1])
match_index <- match(annotations[,1], keys)

gene_length <- gene_length[match_index,]
#RPKM <- RPKM[match_index,]
#rownames(RPKM) <- annotations[,2]
new_annotation <- data.frame(annotations, gene_length)

cell_titles <- read.csv("gene_expression_lung_titles_formatted.txt")
cell_titles_2 <- unique(cell_titles[,1])

####### calculate the gene length #########

express_matrix_tabular <- read.csv("gene_expression_lung.txt", header = T, sep	= "\t")

match_index <- match(rownames(express_matrix_tabular), new_annotation[,2])
express_matrix_tabular <- express_matrix_tabular[which(!is.na(match_index)),]
write.table(express_matrix_tabular, "gene_expression_lung_formatted.txt", col.names = T, row.names = T, quote = F, sep = "\t")
gene_length <- gene_length[match_index[which(!is.na(match_index))]]
write.table(gene_length, "gene_expression_lung_formatted_gene_length.txt", col.names = T, row.names = F, quote = F, sep = "\t")










library(monocle)
cell_titles_2 <- c("stromal_cell", "epithelial_cell_of_lung", "leukocyte", "monocyte", "T_cell", "classical_monocyte", "ciliated_columnar_cell_of_tracheobronchial_tree", "natural_killer_cell")
for (i in 1: length(cell_titles_2))
{
    
    group_1 <- which(cell_titles[,1] == cell_titles_2[i]) ### HFD
    
    group_2 <- which(cell_titles[,1] != cell_titles_2[i])### ND
    
    
    cluster_express_matrix <- express_matrix_tabular[,c(group_1, group_2)]
    
    ################# monocle pipeline #################
    
    ######################## cal rpkm ########################
    
        
    
    d <- DGEList(counts = cluster_express_matrix)
    #cpm(d)
    #cpm(d,log=TRUE)
    
    d$genes$Length <- gene_length
    
    rpkm_express <- rpkm(d)
    
    
    rpkm_express_symbol <- rpkm_express
    
     ################################# calculate Differential expression #######################################
    
    
    group <- c(replicate(length(group_1),"HFD"), replicate(length(group_2),"ND"))
    annot <- data.frame(group) ######## 1 vs 7 big
    
    rownames(annot) <- colnames(rpkm_express_symbol)
        
    colnames(annot) <- "samples"
 
    pd <- new("AnnotatedDataFrame", data = annot)
    
    
    fd_data <- data.frame(rownames(rpkm_express_symbol)) ###### remove the genes with NA symbol
    rownames(fd_data) <- rownames(rpkm_express_symbol)
    
    colnames(fd_data) <- "gene_short_name"
    
        
    fd <- new("AnnotatedDataFrame", data = fd_data)
    
    HSMM <- newCellDataSet(as.matrix(rpkm_express_symbol), phenoData = pd, featureData = fd)
    HSMM <- detectGenes(HSMM, min_expr = 0.1)
        
    HSMM <- estimateSizeFactors(HSMM)
    
    diff_ND_vs_HFD <- differentialGeneTest(HSMM, fullModelFormulaStr="~samples") #########
    
    match_index_2 <- match(rownames(diff_ND_vs_HFD), rownames(rpkm_express_symbol))


	log2fc <- apply(rpkm_express_symbol, 1, function(b) {log2(mean(b[1:length(group_1)])) - log2(mean(b[(length(group_1)+1):(length(group_1)+length(group_2))]))})######## 1 vs 7 big
    
    HFD_num_qualifed <- apply(rpkm_express_symbol, 1, function(b) {length(b[which(b[1:length(group_1)] > 0.1)])})
    ND_num_qualifed <- apply(rpkm_express_symbol, 1, function(b) {length(b[which(b[(length(group_1)+1):(length(group_1)+length(group_2))] > 0.1)])})
    
    result_matrix <- data.frame(diff_ND_vs_HFD,HFD_num_qualifed[match_index_2], ND_num_qualifed[match_index_2], as.matrix(log2fc)[match_index_2,], rpkm_express_symbol[match_index_2,]) #############
    
    ################################# converting gene id #######################################
    
    # annotate the genes
    #keys <- as.character(rownames(result_matrix))
    #annotations <- select(org.Mm.eg.db, keys,  c("GENENAME"), keytype="SYMBOL")
    
    #match_index <- match(annotations[,1], keys)
    #output <- data.frame(annotations, result_matrix[match_index,])
    output <- result_matrix
    filename = paste(cell_titles_2[i], "_vs_rest.txt", sep="")
    write.table(output, filename, row.names = F, col.names=T, quote = F,sep="\t")  #############
    system(paste("sed 's/\\./,/g' ", cell_titles_2[i], "_vs_rest.txt > ", cell_titles_2[i], "_vs_rest_comma.txt", sep = ""))
    
    
    filename = paste(cell_titles_2[i], "_vs_rest_short.txt", sep="")
    write.table(output[,1:9], filename, row.names = F, col.names=T, quote = F,sep="\t")  #############
    system(paste("sed 's/\\./,/g' ", cell_titles_2[i], "_vs_rest_short.txt > ", cell_titles_2[i], "_vs_rest_short_comma.txt", sep = ""))
        
    filename = paste(cell_titles_2[i], "_vs_rest.txt", sep="")
    x <- output
    
    num_hfd = max(x[,7])
    num_nd = max(x[,8])
    









    total_num = length(colnames(x))
    

#if((num_hfd + num_nd + 10) == total_num)
#   {
        
        output_2 <- data.frame( x[,1:9], rowMeans(x[,10:(10 + num_hfd - 1)]), rowMeans(x[,(10 + num_hfd): total_num]))
        colnames(output_2)[c(ncol(output_2)-2, ncol(output_2)-1, ncol(output_2))] <- c("Log2FC","AVERAGE RPKM HFD", "AVERAGE RPKM ND")
        
        write.table(output_2, paste(gsub('.{4}$', '', filename), "_short_with_rpkm.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t" )
        system(    paste("sed 's/\\./,/g' ", paste(gsub('.{4}$', '', filename), "_short_with_rpkm.txt", sep = ""), ">",  gsub ('.{4}$', '_comma.txt', paste(gsub('.{4}$', '', filename), "_short_with_rpkm.txt", sep = "")), sep = ""))
        
        #   }else
        #{
        
        #print ("please check the number of cells MANUALLY!!!")
        
        #}


}














