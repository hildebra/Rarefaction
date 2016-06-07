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


#rmatrs <- lapply(c( 10000),randommatrix) 

path <- "/Users/saary/testData/testdataSingle.csv"
table     <- read.table(file = path, header = TRUE, row.names = 1)
table.m   <- as.matrix(table)

run.vegan <- function(data, samplesize=NULL){
  if(is.null(samplesize)){
    samplesize <- min(colSums(data))
  }
  data2 <- t(data)
	v.ret			<- vegan::rrarefy(x = data2, sample = samplesize)
	return(c(v.ret))
}

run.rare <- function(data, samplesize=NULL){
	if(is.null(samplesize)){
	  samplesize <- min(colSums(data))
	}
  
	
	res 		<- rarefaction::rare(input = data, 
					depth = samplesize, repeats = 1, 
					NoOfMatrices = 1, returnObject = T, verbose = F, margin=2)
	return(c(t(res$raremat[[1]])))
}


#x <- c(sapply(rmatrs, run.rare))
#y <- c(sapply(rmatrs, run.vegan))

x <- run.rare(table.m, 100000)
y <- run.vegan(table.m, 100000)

df <- data.frame(rare=x, vegan=y)
correaltion <- cor(x,y)
dfNull <- df[rowSums(df)>1,]

dfS <- dfNull[sample(1:nrow(dfNull), 1000,replace=FALSE),]
dfS$rare <- as.numeric(dfS$rare)
dfS$vegan <- as.numeric(dfS$vegan)
rownames(dfS) <- NULL
#dfS <- log(dfS, base = 10)

m <- ggplot(dfS, aes(rare, vegan)) +
	xlab("rarefaction values by 'rarefaction'") + 
	ylab("rarefaction values by 'vegan'") +
	ggtitle("Comparison of rarefaction results") +
	geom_point(shape=1) +
  geom_smooth(method = "lm", formula = "y~x") + 
  scale_x_log10()+  scale_y_log10()

ggsave(filename = "/Users/saary/projekt/Rarefaction/tests/scatter_real_ggplot.png", plot = m, device = "png", width = 7.5, height = 7.5)


# simple R plot, works
png(file = "/Users/saary/projekt/Rarefaction/tests/scatter_real_R.png")
reg1 <- lm(rare~vegan, dfS)
par(cex=.8)
plot(dfS$rare, dfS$vegan, log ="xy", main="Comparison of rarefaction results", xlab="rarefaction values by 'rarefaction'", ylab="rarefaction values by 'vegan'")
abline(reg1)
text(1500,200, paste("cor(x,y) = ", round(correaltion,3), sep=""))
dev.off()


# qplot ist auch doof
n <- qplot(rare, vegan, data = dfNull, log="xy") +
  geom_smooth(method="lm")



