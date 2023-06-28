





deseq2_normalizaton <- function(express_matrix, control, wt)
{
    
    library(DESeq2)
    ################# main calculation #################
    
    # setup the sample condition in order to perform wald test
    
    sampleCondition<-c(replicate(control,"treated"),replicate(wt,"untreated"))
    
    
    coldata_u <- data.frame(sampleCondition)
    rownames(coldata_u) <- colnames(express_matrix)
    
    # Loading the samples
    dds <- DESeqDataSetFromMatrix(countData = express_matrix , colData = coldata_u, design = ~ sampleCondition)
    
    # perform the calculation
    dds <- estimateSizeFactors(dds)
    
    # baseMean is average of normalized counts
    
    # The shrinkage is greater for the log2 fold change estimates from genes with low counts and high dispersion
    normCounts <- counts(dds, normalized=TRUE)
    
    return (normCounts)
}









detect.genes <- function(x, cut=1) {
    length(which(x > cut))
}









convert.to.color <- function(x,colscale,col.range = NULL){
    x.range <- range(na.omit(x))
    by=0.1
    if (is.null(col.range)){
        col.range <- seq(x.range[1],x.range[2],by=by)
    }
    col.def <- colscale(length(col.range))
    col.idx <- round((x-x.range[1])/by)+1
    col.idx[col.idx > length(col.range)] <- length(col.range)
    cols <- col.def[col.idx]
    return(list(cols=cols,col.def=col.def,col.range=col.range))
}



run.pca<-function(data,seln=0,log.vals=TRUE,samples.col=TRUE,center=TRUE){
	if (!samples.col){ data<-t(data) }

	# remove zero read genes
	z<-which(rowSums(data)==0)
	if (length(z>0)){
	    data<-data[-z,]
	}

	if (log.vals) { data <- log2(data+0.1) }
	
	# select top varied genes
	if (seln>0){
	    sD<-apply(data,1,sd)
	    o<-order(sD,decreasing=T)
	    data<-data[o[1:seln],]
	}


        myPca <- prcomp(t(data),center=center,scale.=FALSE)
        vars<- myPca$sdev^2
        vars<- vars/sum(vars)
        pcnames<-mat.or.vec(1,length(vars))
        for (i in 1:length(vars)) {
            pcnames[i]<-sprintf("PC%d\n%.5f",i,vars[i])
        }

	myPca$pc.names<-as.vector(pcnames)
        return(myPca)
}


pca.plot <- function (data, seln=0, selpc=c(1,2),cex=0.7, log.vals=TRUE,center=TRUE,samples.col=TRUE, ...){


	 if (class(data) != "prcomp"){
	    data<-run.pca(data,seln=seln,log.vals=log.vals,samples.col=samples.col,center=center)
	 }
	 
	 tmpPca <- as.data.frame(data$x[,selpc])
	 colnames(tmpPca)<-data$pc.names[selpc]
    	 plot(tmpPca, ...)
         invisible(data)
}

mgsub <- function(pattern, replacement, x, ...) {
    if (length(pattern)!=length(replacement)) {
        stop("pattern and replacement do not have the same length.")
    }
    result <- x
    for (i in 1:length(pattern)) {
        result <- gsub(pattern[i], replacement[i], result, ...)

    }
    result
}


colplot.clust <- function( hclust, lab=substr(hclust$labels,1,5), lab.col=rep(1,length(hclust$labels)), hang=0.1,symbol="",...){
 if (length(symbol)>1) {
    lab<-paste(symbol,lab,sep="  ")
 }
 y <- rep(hclust$height,2)
 x <- as.numeric(hclust$merge)
 y <- y[which(x<0)]
 x <- x[which(x<0)]
 x <- abs(x)
 y <- y[order(x)]
 x <- x[order(x)]
 plot( hclust, labels=FALSE, hang=hang, ... )
 text( x=x, y=y[hclust$order]-(max(hclust$height)*hang), labels=lab[hclust$order],cex=.2, col=lab.col[hclust$order], srt=90, adj=c(1,0.5), xpd=NA, ... )
 return (lab.col[hclust$order])
}


top_var <- function (data, seln){
    
    if (seln>0){
        sD<-apply(data,1,sd)
        o<-order(sD,decreasing=T)
        data<-data[o[1:seln],]
    }
    return (data)
}
