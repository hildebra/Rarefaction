rm(list = ls())
require("rarefaction")
require('vegan')

library(microbenchmark)


# benchmark the speed
speed.vegan <- function(data, samplesize){
  result.vegan    <- vegan::rrarefy(x = veganData, sample = samplesize)
  return(result.vegan)
}
speed.rare <- function(data, samplesize){
  result.rarefy   <- rarefaction::rare(input = data, output = "", rareDepth = samplesize, repeats = 10, NoOfMatrices = 1, returnObject = T, verbose = F)
}

path            <- "/Users/saary/projekt/testData/OTU.txt"
output          <- ""
data            <- read.table(file = path, header = TRUE, row.names = 1)
data.rare       <- as.matrix(data)
veganData       <- t(data)
samplesize      <- 1000

microbenchmark(
  #verbose<-rare(path, output, T, 10),
  mute<-speed.vegan(veganData, samplesize),
  mute<-speed.rare(data.rare, samplesize ),
  times=10
)
