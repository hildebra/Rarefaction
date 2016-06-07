path <- '/home/saary/data/testdataSingle.csv'
path <- '/Users/saary/testData/testdataSingle.csv'
# benchmark the ram usage of rarefaction using vagrind


# benchmark the ram usage of rarefaction using vagrind
cat("rarefaction")
require("rarefaction")

args = commandArgs(trailingOnly=TRUE)

#path <- args[2]

result.rarefy   <- rarefaction::rare(input = args[2],
                                     depth = as.numeric(args[1]),
                                     repeats = 1,
                                     NoOfMatrices = 1,
                                     returnObject = T,
                                     verbose = T,
                                 margin = 2)
str(result.rarefy$raremat)
