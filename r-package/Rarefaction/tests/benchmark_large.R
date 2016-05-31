
path <- '/home/saary/rareFalk/testData/OTU.txt'


print("Rarefaction with path, without returning the values")
require("rarefaction")
ptm <- proc.time()
result.rarefy   <- rarefaction::rare(input = path,
                                     output = "",
                                     rareDepth = 1000,
                                     repeats = 10,
                                     NoOfMatrices = 1,
                                     returnObject = F,
                                     verbose = F)
proc.time() - ptm
str(result.rarefy)
rm(list = ls())




print("Rarefaction with path, while returning the values")
path <- '/home/saary/rareFalk/testData/OTU.txt'

require("rarefaction")
ptm <- proc.time()
result.rarefy   <- rarefaction::rare(input = path,
                                     output = "",
                                     rareDepth = 1000,
                                     repeats = 10,
                                     NoOfMatrices = 1,
                                     returnObject = T,
                                     verbose = F)
proc.time() - ptm
#str(result.rarefy)
rm(list = ls())


print("vegan")
path <- '/home/saary/rareFalk/testData/OTU.txt'
require("vegan")
ptm <- proc.time()
data            <- read.table(file = path, header = TRUE, row.names = 1)
samplesize      <- min(rowSums(data))
result.vegan    <- vegan::rrarefy(x = data, sample = samplesize)
proc.time() - ptm
rm(list = ls())
