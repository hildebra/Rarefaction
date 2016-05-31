path <- '/home/saary/data/testdataSingle.csv'
# benchmark the ram usage of rarefaction using vagrind
cat("vegan")
require("vegan")

# load data in matrix
data     		<- t(as.matrix(read.table(file = path, sep = "\t", header = T, row.names = 1)))
samplesize 		<- min(rowSums(data))
result.rarefy   <- rrarefy(x = data, sample = samplesize)
str(result.rarefy)
