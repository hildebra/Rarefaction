path <- '/home/saary/data/testdata.csv'


print("Rarefaction with path, without matrix transformation")
require("rarefaction")
ptm <- proc.time()
result.rarefy   <- rarefaction::rare(input = path,
                                     output = "",
                                     rareDepth = 1000,
                                     repeats = 10,
                                     NoOfMatrices = 1,
                                     returnObject = F,
                                     verbose = T,
									 margin = 1)
proc.time() - ptm
str(result.rarefy)
rm(list = ls())


path <- '/home/saary/data/testdata.csv'


print("Rarefaction with path, with tranforming")
require("rarefaction")
ptm <- proc.time()
result.rarefy   <- rarefaction::rare(input = path,
                                     output = "",
                                     rareDepth = 1000,
                                     repeats = 10,
                                     NoOfMatrices = 1,
                                     returnObject = F,
                                     verbose = T,
									 margin = 1)
proc.time() - ptm
str(result.rarefy)
rm(list = ls())
