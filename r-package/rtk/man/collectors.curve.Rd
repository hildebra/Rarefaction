\name{collectors.curve}
\alias{collectors.curve}



\title{collectors.curve}

\description{ Collectorscurves visualize the richness gained by picking more samples.}

\usage{
	collectors.curve(x, y = NULL, col = 1, times = 10, bin = 3, add = FALSE, 
	                 ylim = NULL, xlim = NULL, doPlot = TRUE, rareD = NULL, 
	                 cls = NULL, pch = 20, col2 = NULL, accumOrder = NULL, ...)
}

\arguments{
  \item{x}{Input a rarefaction object with one matrix and one depth or dataframe/matrix or the output of collectors.curve itself}
  \item{y}{secondary input matrix for comparative plots}
  \item{col}{fill color of the boxplots (set to c(0) for no color)}
  \item{times}{Number of times the sampeling of samples should be perfomed}
  \item{bin}{Number of samples to be added each step. Usefull to adjust for a quick glance.}
  \item{add}{add the plot to an existing plot?}
  \item{ylim}{Limits for Y-scale}
  \item{xlim}{Limits for X-scale}
  \item{doPlot}{should this function plot the collectors curve, or just return an object that can be plotted later with this function?}
  \item{rareD}{Depth to which rarefy the dataset using rtk}
  \item{cls}{vector describing the class of each input sample}
  \item{pch}{Plotting symbols}
  \item{col2}{Color for the border of the boxplot, defaults to col}
  \item{accumOrder}{accumulate successively within each class, given by cls in the order given in this vector. All classes in cls must be represented in this vector.}
  \item{...}{Options passed to plot or boxplot }

}
\details{
The function \code{collectors.curve} can visualize the richness a dataset has, if sampels are picked at random. It can handle rareafaction results as well as normal dataframes.

}

\references{
Saary, Paul, et al. "RTK: efficient rarefaction analysis of large datasets." Bioinformatics (2017): btx206.
}

\seealso{Use \code{\link{plot.rtk}} for how to plot your results.}

\author{Falk Hildebrand, Paul Saary}

\examples{
require("rtk")
# Collectors Curve dataset should be broad and contain many samples (columns)
data       <- matrix(sample(x = c(rep(0, 15000),rep(1:10, 100)),
                     size = 10000, replace = TRUE), ncol = 80)
data.r     <- rtk(data, ReturnMatrix = 1, depth = min(colSums(data)))
# collectors curve on dataframe/matrix
collectors.curve(data, xlab = "No. of samples", ylab = "richness")
# same with rarefaction results (one matrix recommended)
collectors.curve(data.r, xlab = "No. of samples (rarefied data)", ylab = "richness")

# if you want to have an accumulated order, t compare various studies to one another:
cls          <- rep_len(c("a","b","c","d"), ncol(data))  # study origin of each sample
accumOrder   <- c("b","a","d","c")      # define the order, for the plot
colors       <- c(1,2,3,4)
names(colors) <- accumOrder # names used for legend
collectors.curve(data, xlab = "No. of samples",
                 ylab = "richness", col = colors, bin = 1,cls = cls, 
                 accumOrder = accumOrder)

}
