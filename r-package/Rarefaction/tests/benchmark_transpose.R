rm(list = ls())
require("rarefaction")
require("vegan")


path    <- "/Users/saary/projekt/testData/9x9.csv"
output  <- "/Users/saary/Rtest_fri"

d <- 20

print("Normal")
output1 <- rare(path, output, verbose = FALSE, rareDepth = d, 
                repeats = 10, returnObject = T, NoOfMatrices =1, margin=1)
output1$raremat


table     <- as.matrix(read.table(file = path, sep = "\t", header = T, row.names = 1))
output1 <-  vegan::rrarefy(table, 20)
output1

print("Transpose")
output2 <- rare(path, output, verbose = FALSE, rareDepth = d, 
                repeats = 5, returnObject = T, NoOfMatrices =1, margin=2)
output2$raremat
