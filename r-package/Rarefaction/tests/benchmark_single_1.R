path <- '/home/saary/data/testdataSingle.csv'
# benchmark the ram usage of rarefaction using vagrind
cat("rarefaction")
require("rarefaction")
result.rarefy   <- rarefaction::rare(input = path,
                                     output = "",
                                     rareDepth = 1000,
                                     repeats = 10,
                                     NoOfMatrices = 1,
                                     returnObject = T,
                                     verbose = F,
									 margin = 1)
str(result.rarefy$raremat)
