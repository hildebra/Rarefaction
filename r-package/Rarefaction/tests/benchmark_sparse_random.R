require("rarefaction")
require('vegan')
require("Matrix")
require("reshape2")
require("ggplot2")


comp <- function(x){
  
  data<- as.matrix(round(rsparsematrix(nrow = 100,ncol = 75, density = 0.3, rand.x = runif)*100, digits = 0))
  
  rownames(data) <- rep("bpb", nrow(data))
  
  
  samplesize      <- min(rowSums(data))
  result.vegan    <- vegan::rrarefy(x = data, sample = samplesize)
  result.rarefy   <- rarefaction::rare(input = t(data), output = "", rareDepth = samplesize, repeats = 3, NoOfMatrices = 1, 
                                       returnObject = T, verbose = F)
  
  x <- cor(c(result.vegan),as.numeric(unlist(result.rarefy$countsDF)))
  return(x)
}

x <- seq(1,20)

y <- sapply(x, comp)

df <- as.data.frame(x=x, y=y)
ggplot(df, aes(x,y)) +
  geom_point() +
  
