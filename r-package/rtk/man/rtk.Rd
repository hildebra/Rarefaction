\name{rtk}
\alias{rtk}


\title{Rarefy tables}

\description{ Rarefy datasets in R or from a path. }

\usage{
rtk(input, repeats = 10, depth = 0, ReturnMatrix = 0, margin = 2,
    verbose = FALSE, threads = 1, tmpdir = NULL )
}

\arguments{
  \item{input}{This can be either a numeric matrix or a path to a text file in tab-delimited format on the locally available storage. The later option is for very big matrices, to avoid unnecessary memory consumption in R.}
  \item{repeats}{Number of times to compute diversity measures. (\code{default: 10})}
  \item{depth}{Number of elements per row/column to rarefy to. The so called rarefaction depth or samplesize. Can also be a vector of ints. (\code{default: 1000})}
  \item{ReturnMatrix}{Number of rarefied matrices which are returned to R. Set to zero to only measure diversity. (\code{default: 1})}
  \item{margin}{Indicates which margin in the matrix represents the Samples and Species. Default is to rarefy assuming columns represent single samples (margin=2). If margin=1, rows are assumed to be samples. (default: 2 (columns))}
  \item{verbose}{If extra output should be printed to std::out or not to see progress of rarefaction. (\code{default: TRUE)}}
  \item{threads}{Number of threads to use during rarefaction}
  \item{tmpdir}{Location to store temporary files}

}
\details{
Function \code{rare} takes a dataset and calcualtes the diversity measures, namely the shannon diversity, richness, simpson index, the inverse simpson index, chao1 and evenness.

If wished for the function can also return one or multiple rarefied matrices rarefied to one or multiple depths. Those can then also be used to create collectorcurves (see \code{\link{collectors.curve}}).
}

\value{The function \code{rare} returns an object of class 'rarefaction',  containing the  objects \code{divvs}, \code{raremat}, \code{skipped}, \code{div.median} and \code{depths}. If more than one depth was computed the elements 1-4 are inside a list themself and can be acessed by the index of the desired depth.

The object \code{divvs} contains a list of diversity measures for each sample provieded.

\code{raremat} is one or multiple rarefied matrices. Samples with not enough counts are removed, thus not all \code{raremat}-matrices for different depths might be of the same size. If and which sampels where excluded is denoted in the element \code{skipped} using the names of the respective samples.

\code{depths} just contains the input variable and might be usefull for further analysis of the results.

It is possible to plot the results of the rarefaction, depending on the parameters passed to \code{rare}. See \code{\link{plot.rtk}} for examples.


}

\references{
To be filled someday
}

\seealso{\code{\link{plot.rtk}}, \code{\link{collectors.curve}}}

\author{\packageAuthor{rtk}}

\examples{
require("rtk")
# generate semi sparse example data
data            <- matrix(sample(x = c(rep(0, 1500),rep(1:10, 500),1:1000),
                          size = 120, replace = TRUE), 10)
# find the column with the lowest aboundance
samplesize      <- min(colSums(data))
# rarefy the dataset, so each column contains the same number of samples
data.rarefied   <- rtk(input = data, depth = samplesize, ReturnMatrix = 1)

richness   <- get.diversity(data.rarefied, div = "richness")
eveness    <- get.diversity(data.rarefied, div = "eveness")

}
