\name{rare}
\alias{rare}


\title{Rarefy tables}

\description{ Rarefy datasets in R or from a path. }

\usage{
rare(x, sample, .. options with defaults MARGIN = 2)
}

\arguments{
  \item{input}{Either a R object like a numeric matrix or a path to a \code{.csv}-file. (\code{required})}
  \item{repeats}{Number of times to compute diversity measures. (\code{default: 10})}
  \item{depth}{Number of elements per row/column to rarefy to. The so called rarefaction depth or samplesize. Can also be a vector of ints. (\code{default: 1000})}
  \item{NoOfMatrices}{Number of rarefied matrices which are returned to R. Set to zero to only measure diversity. (\code{default: 1})}
  \item{margin}{In which direction the dataset will be rarefied.
  Default is to rarefy all columns (margin=2). Can be set to rows by
  \code{margin=1} to get the same number of species in each row instead. (\code{default: 2 (columns)})}
  \item{verbose}{If extra output should be printed to std::out or not to see progress of rarefaction. (\code{default: TRUE)}}

}
\details{

  Function \code{rare} takes a datset and calcualtes the diversity measures, namely the shaonnon diversity, richness,  simpson index, the inverse simpson index, chao1 and eve.

}

\value{
The function \code{rare} returns an list containing a list of matrices with the length specified in the options. All matrices are rarefied to the desired depth.
Moreover it always returns the diversity measures for all samples (depending on the margin setting) with the number of repeats specified in repeats.
A median for all measurements is also returned.
}

\references{
  Our paper here

}

\seealso{Nowhere? or vegan}

\author{Falk Hildebrand, Paul Saary}

\examples{
require("rarefaction")

# generate semi sparse example data
data            <- matrix(sample(x = c(rep(0, 1500),rep(1:10, 500),1:1000),size = 120, replace = TRUE), 10)
# find the column with the lowest aboundance
samplesize      <- min(colSums(data))
# rarefy the dataset, so each column contains the same number of samples
data.rarefied   <- rare(input = data, depth = samplesize, NoOfMatrices = 1, verbose=FALSE)

#path 		<- "/path/to/a/file.csv"
#data.rarefied   <- rare(input = path, depth = 1000)


}
