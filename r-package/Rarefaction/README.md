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
R CMD check Rarefaction_*.tar.gz
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
data.rarefied   <- rare(input = data, rareDepth = samplesize, NoOfMatrices = 1)

path 		<- "/path/to/a/file.csv"
data.rarefied   <- rare(input = path, rareDepth = 1000)

```
More documentation if provided inside of the R package. Please look into `man/`.

## Development
To develop this package it is recommended to instal the packages `Rcpp` and  `testthat` for R.
Clone the git repo. The main software `rare` is located in `Rarefaction/Rare`. The files used in the r-package are link via symbolic links. This might not work on your system, so please bear in mind, that changes on those files should always committed so that this structure is preserved, even if your system might not support this.

Unit tests are performed using testthat in the `test/testthat/` location. Read up on unit tests before adding new one.
