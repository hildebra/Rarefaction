
path <- '/Users/saary/projekt/testData/mh14/mh.14.gene.profile.screened.hg19.on.1000RefGeneCat.1-2.fastx.allbest.l30.p95.base.only.unique.raw.gene'


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
path <- '/Users/saary/projekt/testData/mh14/mh.14.gene.profile.screened.hg19.on.1000RefGeneCat.1-2.fastx.allbest.l30.p95.base.only.unique.raw.gene'

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
path <- '/Users/saary/projekt/testData/mh14/mh.14.gene.profile.screened.hg19.on.1000RefGeneCat.1-2.fastx.allbest.l30.p95.base.only.unique.raw.gene'
require("vegan")
ptm <- proc.time()
data            <- read.table(file = path, header = TRUE, row.names = 1)
samplesize      <- min(rowSums(data))
result.vegan    <- vegan::rrarefy(x = data, sample = samplesize)
proc.time() - ptm
rm(list = ls())






