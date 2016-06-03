\name{rarefaction-package}
\alias{rarefaction-package}
\alias{rarefaction}
\docType{package}
\title{
\packageTitle{rarefaction}
}
\description{
\packageDescription{rarefaction}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{rarefaction}
\packageIndices{rarefaction}
This package might be used to rarefy data. Providing a table containing columns or rows with
different sums it returns an object containing one or multiple rarefied version of the table.

Also entropy

Quick by use of Rcpp and Cpp11

}
\author{
\packageAuthor{rarefaction}

Maintainer: \packageMaintainer{rarefaction}
}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
require("rarefaction")

# generate semi sparse example data
data            <- matrix(sample(x = c(rep(0, 1500),rep(1:10, 500),1:1000),size = 120, replace = T), 10)
# find the column with the lowest aboundance
samplesize      <- min(colSums(data))
# rarefy the dataset, so each column contains the same number of samples
data.rarefied   <- rare(input = data, rareDepth = samplesize, returnObject = T, NoOfMatrices = 1)

}
