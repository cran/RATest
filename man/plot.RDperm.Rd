% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.RDperm.R
\name{plot.RDperm}
\alias{plot.RDperm}
\title{Plot RDperm}
\usage{
\method{plot}{RDperm}(x, w, plot.class = "both", ...)
}
\arguments{
\item{x}{Object of class "RDperm"}

\item{w}{Character. Name of variable to be plotted}

\item{plot.class}{Character. Can be: "both" for a histogram and cdf plot, "hist" for a histogram or "cdf" for only the cdf plot}

\item{...}{Additional ggplot2 controls}
}
\description{
Plots a histogram and empirical cdf
}
\examples{
\dontrun{
permtest<-RDperm(W=c("demshareprev","demwinprev"),z="difdemshare",data=lee2008)
plot(permtest,w="demshareprev")
}
}
\references{
Canay, I and Kamat V, (2017) Approximate Permutation Tests and Induced Order Statistics in the Regression Discontinuity Design. \url{http://faculty.wcas.northwestern.edu/~iac879/wp/RDDPermutations.pdf}
}
\author{
Maurcio Olivares Gonzalez

Ignacio Sarmiento Barbieri
}
\keyword{permutation}
\keyword{rdperm}