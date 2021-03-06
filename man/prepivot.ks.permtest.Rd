% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepivot.ks.permtest.R
\name{prepivot.ks.permtest}
\alias{prepivot.ks.permtest}
\title{Permutation Test for the two-sample goodness-of-fit problem under covariate-adaptive randomization}
\usage{
prepivot.ks.permtest(Y1, Y0, alpha, B, n.perm)
}
\arguments{
\item{Y1}{Numeric. A vector containing the response variable of the treatment group.}

\item{Y0}{Numeric. A vector containing the response variable of the control group.}

\item{alpha}{Numeric. Nominal level for the test. The default is 0.05.}

\item{B}{Numeric. Number of weighted bootstrap samples.}

\item{n.perm}{Numeric. Number of permutations needed for the stochastic approximation of the p-values. The default is n.perm=999.}
}
\value{
An object of class "prepivot.ks.permtest" containing at least the following components:

 \item{n_populations}{Number of grups.}
 \item{N}{Sample Size.}
 \item{T.obs}{Observed test statistic.}
 \item{cv}{Critical Value. This value is used in the general construction of a randomization test.}
 \item{pvalue}{P-value.}
 \item{rejectrule}{Rule. Binary decision for randomization test, where 1 means "to reject"}
 \item{T.perm}{Vector. Test statistic recalculated for all permutations used in the stochastic approximation.}
 \item{n.perm}{Number of permutations.}
 \item{B}{Bayesian bootstrap samples.}
 \item{sample_sizes}{Groups size.}
}
\description{
A permutation test of the two-sample goodness-of-fit hypothesis when the randomization scheme is covariate-adaptive. The permutation test considered here is based on prepivoting the Kolmogorov-Smirnov test statistic following Beran (1987,1988), and adapted by Olivares (2020). Current version includes the following randomization schemes: simple randomization, Efron's biased-coin design, Wei's biased-coin design, and stratified block randomization. This implementation uses a Bayesian bootstrap approximation for prepivoting.
}
\examples{
\dontrun{
Y0 <- rnorm(100, 1, 1)
Y1 <- rbeta(100,2,2)
Tx = sample(100) <= 0.5*(100)
# Observed Outcome 
Y = ifelse( Tx, Y1, Y0 )
dta <- data.frame(Y = Y, A = as.numeric(Tx))
pKS.GoF<-prepivot.ks.permtest(dta$Y[dta$A==1],dta$Y[dta$A==0],alpha=0.05,B=1000,n.perm = 999)
summary(pKS.GoF)
}
}
\references{
Beran, R. (1987). Prepivoting to reduce level error of confidence sets. Biometrika, 74(3): 457–468.
Beran, R. (1988). Prepivoting test statistics: a bootstrap view of asymptotic refinements. Journal of the American Statistical Association, 83(403):687–697.
Olivares, M. (2020). Asymptotically Robust Permutation Test under Covariate-Adaptive Randomization. Working Paper.
}
\author{
Maurcio Olivares
}
\keyword{Bayesian}
\keyword{bootstrap}
\keyword{goodness-of-fit}
\keyword{permutation}
\keyword{prepivoting}
\keyword{test}
