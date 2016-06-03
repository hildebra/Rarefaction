require("rarefaction")

# generate semi sparse example data
data            <- matrix(sample(x = c(rep(0, 1500),rep(1:10, 500),1:1000),size = 120, replace = T), 10)
# find the column with the lowest aboundance
samplesize      <- min(colSums(data))
# rarefy the dataset, so each column contains the same number of samples
data.rarefied   <- rare(input = data, rareDepth = samplesize, returnObject = T, NoOfMatrices = 1)


