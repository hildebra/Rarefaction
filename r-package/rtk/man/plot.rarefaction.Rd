\name{plot.rarefaction}
\alias{rarefactioncurve}
\alias{plot}
\alias{plot.rarefaction}


\title{Plot rarfeaction results}

\description{ Rarefy datasets in R or from a path. }

\usage{
plot.rarefaction(x, div = c("richness"),  groups = NA, col = NULL, lty = 1, pch = NA, fit = "arrhenius", legend = TRUE, legend.pos = "topleft", log.dim = "", boxplot = FALSE, ...)
plot(...)
}

\arguments{
  \item{x}{a rare result object}
  \item{div}{Diversity measure to plot. Can be any of \code{c('richness', 'shannon', 'simpson', 'invsimpson', 'chao1', 'eve')}}
  \item{groups}{If grouping is desired a vector of factors corresponting to the input samples}
  \item{col}{Colors used for plotting. Can be a vector of any length which will be recycled if it is to small. By default a rainbow is used.}
  \item{lty}{Linetypes used for plotting. Can be a vector of any length which will be recycled if it is to small.}
  \item{pch}{Symbols used for plotting. Can be a vector of any length which will be recycled if it is to small.}
  \item{fit}{Fit the rarefaction curve. Possible values: \code{c("arrhenius", "michaelis-menten", "logis")}}
  \item{legend}{Logical indicating if a legend should be created or not}
  \item{legend.pos}{Position of the said legend}
  \item{log.dim}{Character vector indicating which scale log log transform for plotting rarefaction curves.}
  \item{...}{Other plotting input will be passed to \code{plot} or \code{boxplot} repectivly}

}
\details{
To create plots from the rarefaction results you can easily just call a plot on the resulting elements. This will either produce a rarefaction curve, if mor than one depth was rarefied to, or a boxplot for a single depth. Grouping of samples is possible by simply passing a vetor of the length of the samples to the option \code{groups}.

Rarefaction curves can be fittet to either the arrhenius-equation, the michaelis-menten (\link{SSmicmen}) equation or the logis function \link{SSlogis}. To disable fitting \code{fit} must be set to \code{FALSE}.
}



\references{
  Our paper here

}


\seealso{\code{\link{rtk}}, \code{\link{collectors.curve}}}

\author{Falk Hildebrand, Paul Saary}

\examples{
require("rarefaction")
# generate semi sparse example data
data            <- matrix(sample(x = c(rep(0, 1500),rep(1:10, 500),1:1000),size = 120, replace = TRUE), 40)
# find the column with the lowest aboundance
samplesize      <- min(colSums(data))
# rarefy the dataset, so each column contains the same number of samples
d1  <- rtk(input = data, depth = samplesize)
# rarefy to different depths between 1 and samplesize
d2  <- rtk(input = data, depth = round(seq(1, samplesize, length.out = 200)))

# just the richness of all three samples as boxplot
plot(d1, div = "richness")
#rarefaction curve for each sample with fit
plot(d2, div = "eveness", fit = "arrhenius", pch = c(1,2,3))
# Rarefaction curve with boxplot, sampels pooled together (grouped)
plot(d2, div = "richness", fit = FALSE, boxplot = TRUE, col = 1, groups = rep(1, ncol(data)))

}
