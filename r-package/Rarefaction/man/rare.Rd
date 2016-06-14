\name{rare}
\alias{rare}


\title{Rarefy tables}

\description{ Rarefy datasets in R or from a path. }

\usage{
rare(input, repeats=10, depth = 1000, ReturnMatrix  = 0, margin = 2, verbose=TRUE, threads=1 )
}

\arguments{
  \item{input}{This can be either a numeric matrix or a path to a text file in tab-delimited format on the harddrive. The later option is for very big matrices, to avoid unnecessary memory consumption in R. (required)}
  \item{repeats}{Number of times to compute diversity measures. (\code{default: 10})}
  \item{depth}{Number of elements per row/column to rarefy to. The so called rarefaction depth or samplesize. Can also be a vector of ints. (\code{default: 1000})}
  \item{ReturnMatrix}{Number of rarefied matrices which are returned to R. Set to zero to only measure diversity. (\code{default: 1})}
  \item{margin}{Indicates which margin in the matrix represents the Samples and Species. Default is to rarefy assuming columns represent single samples (margin=2). If margin=1, rows are assumed to be samples. (default: 2 (columns))}
  \item{verbose}{If extra output should be printed to std::out or not to see progress of rarefaction. (\code{default: TRUE)}}

}
\details{
Function \code{rare} takes a datset and calcualtes the diversity measures, namely the shannon diversity, richness, simpson index, the inverse simpson index, chao1 and eve.
The returned diversity indexes are:
-	Shannon diversity – this is.. blab la, formulat
-	Richness – bla bla

}

\value{The function \code{rare} returns an object of class 'rarefaction',  containing the following objects: $divvs with sublists of the different diversity indexes with N repeats (N being controlled by parameter "repeats"), $raremat are the returned  different rarefied matrices as specified by parameter ReturnMatrix. A median for all diversity measurements is also returned in $div.median.
}

\references{
  Our paper here

}

\seealso{Use \code{\link{plot.rarefaction}} for how to plot your results.}

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
