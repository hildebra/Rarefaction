rm(list = ls())
require("rarefaction")


path    <- "/Users/saary/projekt/testData/9x9.csv"
output  <- "/Users/saary/Rtest_fri"

print("Normal")
output1 <- rare(path, output, verbose = FALSE, rareDepth = 2, 
                repeats = 10, returnObject = T, NoOfMatrices =1, transpose=F)
output1$raremat

print("Normal")
output1 <- rare(path, output, verbose = FALSE, rareDepth = 2, 
                repeats = 10, returnObject = T, NoOfMatrices =1, transpose=T)
output1$raremat