path <- '/Users/saary/testData/OTU.txt'
#path <- '/Users/saary/testData/testdataSingle.csv'

# benchmark the ram usage of rarefaction using vagrind
cat("vegan")
args = commandArgs(trailingOnly=TRUE)
cat(args[1])

require("vegan")


repeats     <- args[3]
path        <- args[2]
samplesize  <- args[1]

# load data in matrix
data     		<- t(read.table(file = path, sep = "\t", header = T, row.names = 1))

i <- 0
while(i < repeats){
  result.rarefy   <- vegan::rrarefy(x = data, sample = samplesize)
  shannonn        <- diversity(x = data, index = "shannon")
  simpson         <- diversity(x = data, index = "simpson")
  invsimpson      <- diversity(x = data, index = "invsimpson")
  richness        <- rarefy(x = data, sample = samplesize)
  i <- i + 1
}


cat("vegan, done")
