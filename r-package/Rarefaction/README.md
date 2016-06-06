# Rarefaction - R module
This R module is a wrapper for the rarefaction program.

## Dependencies
This package depends on the R package RCPP, which enables using Cpp-Code easily in R
packages.

## How to install
To build and install this package from source, just run these commands in a directory
above the dir containing the source files.

```
R CMD build Rarefaction
R CMD install Rarefaction_*.tar.gz
```

## Running the package
```R
require("rarefaction")
# generate semi sparse example data
data            <- matrix(sample(x = c(rep(0, 1500),rep(1:10, 500),1:1000),size = 120, replace = T), 10)
# find the column with the lowest aboundance
samplesize      <- min(colSums(data))
# rarefy the dataset, so each column contains the same number of samples
data.rarefied   <- rare(input = data, rareDepth = samplesize, returnObject = T, NoOfMatrices = 1)

path 		<- "/path/to/a/file.csv"
data.rarefied   <- rare(input = path, rareDepth = 1000)

```

## State of the project
The project is currently developed in the development branch.
No stable version is currently provided.
