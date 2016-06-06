\name{rare}
\alias{rare}


\title{Rarefy tables}

\description{ Rarefy datasets in R or from a path. }

\usage{
rare(x, sample, .. options with defaults MARGIN = 2)
}

\arguments{
  \item{input}{Either a R object like a numeric matrix or a dataframe or a path to a \code{.csv}-file.}
  \item{output}{Should we not remove this option?}
  \item{repeats}{Number of times to compute diversity measures.}
  \item{depth}{Number of elements per row/column to rarefy to. The so called rarefaction depth or samplesize.}
  \item{NoOfMatrices}{Number of rarefied matrices which are returned to R.}
  \item{verbose}{If extra output should be printed to std::out or not to see progress of rarefaction.}
  \item{returnObject}{Should be removed and TRUE by default?}
  \item{margin}{In which direction should the dataset be rarefied.
  Default is to rarefy all columns. Can be set to rows by
  \code{margin=1} to get the same number of species in each row instead.}

}
\details{

  Function \code{rare}  does stuff

}

\value{
  Returns a list containing ....
}

\references{
  Our paper here

}

\seealso{Nowhere? or vegan}

\author{Falk Hildebrand, Paul Saary}

\examples{
require("rarefaction")

# generate semi sparse example data
data            <- matrix(sample(x = c(rep(0, 1500),rep(1:10, 500),1:1000),size = 120, replace = T), 10)
# find the column with the lowest aboundance
samplesize      <- min(colSums(data))
# rarefy the dataset, so each column contains the same number of samples
data.rarefied   <- rare(input = data, rareDepth = samplesize, returnObject = T, NoOfMatrices = 1)

path 		<- "/path/to/a/file.csv"
data.rarefied   <- rare(input = path, rareDepth = 1000)


}
