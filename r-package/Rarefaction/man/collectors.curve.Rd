\name{collectors.curve}
\alias{collectors.curve}



\title{collectors.curve}

\description{ Collectorscurves visualize the richness gained by picking more samples.}

\usage{
	collectors.curve(x, col = 0, times = 10, bin = 5 ...)
}

\arguments{
  \item{x}{Input a rarefaction object with one matrix and one depth or dataframe/matrix}
  \item{col}{color of the boxplots}
  \item{times}{Number of times the sampeling of samples should be perfomed}
  \item{bin}{Number of samples to be added each step. Usefull to adjust for a quick glance.}
  \item{...}{Options passed to plot or boxplot }

}
\details{
The function \code{collectors.curve} can visualize the richness a dataset has, if sampels are picked at random. It can handle rareafaction results as well as normal dataframes.

}

\references{
  Our paper here

}

\seealso{Use \code{\link{plot.rarefaction}} for how to plot your results.}

\author{Falk Hildebrand, Paul Saary}

\examples{
require("rarefaction")
# Collectors Curve dataset should be broad and contain many samples (columns)
data            <- matrix(sample(x = c(rep(0, 15000),rep(1:10, 100)),size = 10000, replace = TRUE), ncol = 80)
data.r 			<- rare(data, ReturnMatrix = 1, depth = min(colSums(data)))
# collectors curve on dataframe/matrix
collectors.curve(data, xlab = "No. of samples", ylab = "richness")
# same with rarefaction results (one matrix recommended)
collectors.curve(data.r, xlab = "No. of samples (rarefied data)", ylab = "richness")

}
