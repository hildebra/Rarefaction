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
#str(output2)




#print("Test different orientations:")
#m.norm  <- rare(table.m, output, verbose = FALSE, returnObject = T, repeats = 10, NoOfMatrices = 3)
#m.t     <- rare(t(table.m), output, verbose = FALSE, returnObject = T, repeats = 10, NoOfMatrices = 3)
#str(m.norm$countsDF)




#print("Test: wrong path")
#path <- "/Users/saary/notamatrix"
#rare(path, output, T)



print("done")
