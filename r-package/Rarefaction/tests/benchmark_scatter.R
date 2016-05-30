rm(list = ls())
require("rarefaction")
require('vegan')
library(reshape2)
require("Matrix")
require(ggplot2)





randommatrix <- function(x){
	input <- as.matrix(round(rsparsematrix(100,700, 0.1, rand.x=runif) * x, 0))
	if(is.null(colnames(input))){
		colnames(input) <- paste("col ", seq(1:ncol(input)), sep="")
	}
	if(is.null(rownames(input))){
		rownames(input) <- paste("row ", seq(1:nrow(input)), sep="")
	}
	return(input)
}


rmatrs <- lapply(c( 10000),randommatrix) 





run.vegan <- function(data, samplesize=NULL){
  if(is.null(samplesize)){
    samplesize <- min(rowSums(data))
  }
	v.ret			<- vegan::rrarefy(x = data, sample = samplesize)
	return(c(v.ret))
}

run.rare <- function(data, samplesize=NULL){
	if(is.null(samplesize)){
	  samplesize <- min(rowSums(data))
	}
  
	data <- t(data)
	res 		<- rarefaction::rare(input = data, output = "", 
					rareDepth = samplesize, repeats = 1, 
					NoOfMatrices = 1, returnObject = T, verbose = F)
	return(c(res$raremat[[1]]))
}


x <- c(sapply(rmatrs, run.rare))
y <- c(sapply(rmatrs, run.vegan))



df <- data.frame(rare=x, vegan=y)
m <- ggplot(df, aes(rare, vegan)) +
	xlab("rarefaction values by 'rarefaction'") + 
	ylab("rarefaction values by 'vegan'") +
	ggtitle("Comparison of rarefaction results") +
	scale_x_log10(breaks = c(10,100,1000,10000))+
	scale_y_log10(breaks = c(10,100,1000,10000))+
  coord_fixed()+
	geom_point(shape=1)
	
ggsave(filename = "/Users/saary/projekt/Rarefaction/tests/scatter.pdf", plot = m, device = "pdf", width = 7.5, height = 7.5)
