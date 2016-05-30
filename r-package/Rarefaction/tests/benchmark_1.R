rm(list = ls())
require("rarefaction")
require('vegan')
library(reshape2)


require(ggplot2)

perform.analysis <- function(x){
  # rarefy matrices
  path            <- "/Users/saary/projekt/testData/OTU.txt"
  output          <- ""
  data            <- read.table(file = path, header = TRUE, row.names = 1)
  
  veganData       <- t(data)
  samplesize      <- 1000
  result.vegan    <- vegan::rrarefy(x = veganData, sample = samplesize)
  result.rarefy   <- rarefaction::rare(input = path, output = "", rareDepth = samplesize, repeats = 10, NoOfMatrices = 1, returnObject = T, verbose = F)
  
  df.vegan        <- as.data.frame(t(result.vegan))
  df.rarefy       <- as.data.frame(result.rarefy$countsDF)
  
  # remove columns which are ommited by ou package
  df.vegan        <- df.vegan[,!(!(colnames(df.vegan) %in% colnames(df.rarefy)))]
  
  return(df.vegan, - df.rarefy)
}



perform.analysis.cor <- function(x){
  # rarefy matrices
  path            <- "/Users/saary/projekt/testData/OTU.txt"
  output          <- ""
  data            <- read.table(file = path, header = TRUE, row.names = 1)
  
  veganData       <- t(data)
  samplesize      <- 1000
  result.vegan    <- vegan::rrarefy(x = veganData, sample = samplesize)
  result.rarefy   <- rarefaction::rare(input = path, output = "", rareDepth = samplesize, repeats = 1, NoOfMatrices = 1, returnObject = T, verbose = F)
  
  df.vegan        <- as.data.frame(t(result.vegan))
  df.rarefy       <- as.data.frame(result.rarefy$countsDF)
  
  # remove columns which are ommited by ou package
  df.vegan        <- df.vegan[,!(!(colnames(df.vegan) %in% colnames(df.rarefy)))]
  
  return(cor(df.vegan, df.rarefy))
}


tests <- rep(NA, 5)


#comparison.cor <- lapply(tests, perform.analysis )
#comparison.m <- matrix(unlist(comparison, F, F), ncol = length(comparison[[1]]), byrow = TRUE)
#a <- data.frame(diff=apply(comparison.m, 2, mean))
#m <- ggplot(a,aes(a$diff)) +
#  geom_histogram(binwidth = 0.1) +
#  xlab("Difference to vegan") +
#  ylab("Count")
#ggsave(filename = "/Users/saary/projekt/Rarefaction/tests/accuaracy.pdf", plot = m, device = "pdf", width = 10, height = 10)





#comparison <- lapply(tests, perform.analysis)
comparison.cor  <- lapply(tests, perform.analysis.cor)
c.array         <- array(unlist(comparison.cor), c(70,70, length(tests)))
c.means         <- rowMeans( c.array , dims = 2 )
c.means.m       <- melt(c.means)
m <- ggplot(data = c.means.m, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
ggsave(filename = "/Users/saary/projekt/Rarefaction/tests/heatmap.pdf", plot = m, device = "pdf", width = 10, height = 10)



