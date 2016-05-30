# type test
require("rarefaction")

# rarefy matrices
path      <- "/Users/saary/projekt/testData/OTU.txt"
output    <- "/Users/saary/Rtest_fri"
data      <- read.table(file = path, header = TRUE, row.names = 1)
data.mat   <- data.matrix(data)
samplesize <- 10

print("DF")
result.rarefy<- rarefaction::rare(input = data, output = "", repeats = samplesize, NoOfMatrices = 1, returnObject = T, verbose = F)
str(result.rarefy$countsDF)
print("Matrix")
result.rarefy<- rarefaction::rare(input = data.mat, output = "", repeats = samplesize, NoOfMatrices = 1, returnObject = T, verbose = F)
str(result.rarefy)
