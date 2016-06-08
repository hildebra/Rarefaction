rm(list = ls())


require("rarefaction")


path    <- "/Users/saary/projekt/testData/OTU.txt"
output  <- "/Users/saary/Rtest_fri"



print("Test: Rarefaction with path")
output1 <- rare(path,  verbose = T, depth = 1000,
                repeats = 5, returnObject = T, NoOfMatrices = 3)
str(output1)


print("Test: Rarefaction with matrix")
print(" Read table")
table     <- read.table(file = path, header = TRUE, row.names = 1)
table.m   <- data.matrix(table)
print("call rare")
output2 <- rare(table.m, verbose = FALSE, depth = 1000, repeats = 10, returnObject = T)




print("done")
