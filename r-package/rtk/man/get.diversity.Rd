\name{get.diversity}
\alias{get.diversity}
\alias{get.mean.diversity}
\alias{get.median.diversity}



\title{get.diversity}

\description{ Collectorscurves visualize the richness gained by picking more samples.}

\usage{
	get.diversity(obj, div = "richness", multi = FALSE)
	get.mean.diversity(obj, div = "richness")
	get.median.diversity(obj, div = "richness")
}

\arguments{
  \item{obj}{Object of type rtk}
  \item{div}{diversity measure as string e.g "richness"}
	\item{multi}{Argument set to true if called recursivly and class should not be checked. Should not be set in normal use case.}
}
\details{
This set of functions allows fast and easy access to calculated diversity measures by rtk.
It returns a matrix, when rarefaction was only performed to one
depth and a list of matrices or vectors if rarefaction was done for multiple depths.
}

\references{
Saary, Paul, et al. "RTK: efficient rarefaction analysis of large datasets." Bioinformatics (2017): btx206.	
}

\seealso{Use \code{\link{rt}} before calling this function.}

\author{Falk Hildebrand, Paul Saary}

\examples{
require("rtk")
# Collectors Curve dataset should be broad and contain many samples (columns)
data            <- matrix(sample(x = c(rep(0, 15000), rep(1:10, 100)),
                          size = 10000, replace = TRUE), ncol = 80)
data.r 			<- rtk(data, depth = min(colSums(data)))
get.diversity(data.r)
get.median.diversity(data.r)
get.mean.diversity(data.r)
}
