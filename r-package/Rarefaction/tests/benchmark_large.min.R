path <- '/home/saary/data/testdata.csv'


print("Rarefaction with path, without returning the values")
require("rarefaction")
ptm <- proc.time()
result.rarefy   <- rarefaction::rare(input = path,
                                     output = "",
                                     rareDepth = 1000,
                                     repeats = 10,
                                     NoOfMatrices = 1,
                                     returnObject = F,
                                     verbose = F,
									 margin = 1)
proc.time() - ptm
str(result.rarefy)
rm(list = ls())
