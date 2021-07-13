
#' @title  Permutation Test for the two-sample goodness-of-fit problem under covariate-adaptive randomization
#'
#' @description A permutation test of the two-sample goodness-of-fit hypothesis when the randomization scheme is covariate-adaptive. The permutation test considered here is based on prepivoting the Kolmogorov-Smirnov test statistic following Beran (1987,1988), and adapted by Olivares (2020). Current version includes the following randomization schemes: simple randomization, Efron's biased-coin design, Wei's biased-coin design, and stratified block randomization. This implementation uses a Bayesian bootstrap approximation for prepivoting.
#' 
#' @param Y1 Numeric. A vector containing the response variable of the treatment group.
#' @param Y0 Numeric. A vector containing the response variable of the control group.  
#' @param alpha Numeric. Nominal level for the test. The default is 0.05.
#' @param B Numeric. Number of weighted bootstrap samples.
#' @param n.perm Numeric. Number of permutations needed for the stochastic approximation of the p-values. The default is n.perm=999.

#' @return An object of class "prepivot.ks.permtest" containing at least the following components:
#' 
#'  \item{n_populations}{Number of grups.}
#'  \item{N}{Sample Size.}
#'  \item{T.obs}{Observed test statistic.}
#'  \item{cv}{Critical Value. This value is used in the general construction of a randomization test.}
#'  \item{pvalue}{P-value.}
#'  \item{rejectrule}{Rule. Binary decision for randomization test, where 1 means "to reject"}
#'  \item{T.perm}{Vector. Test statistic recalculated for all permutations used in the stochastic approximation.}
#'  \item{n.perm}{Number of permutations.}
#'  \item{B}{Bayesian bootstrap samples.}
#'  \item{sample_sizes}{Groups size.}
#'  
#' @author Maurcio Olivares
#' @references  
#' Beran, R. (1987). Prepivoting to reduce level error of confidence sets. Biometrika, 74(3): 457–468.
#' Beran, R. (1988). Prepivoting test statistics: a bootstrap view of asymptotic refinements. Journal of the American Statistical Association, 83(403):687–697.
#' Olivares, M. (2020). Asymptotically Robust Permutation Test under Covariate-Adaptive Randomization. Working Paper.
#' 
#' @keywords permutation test goodness-of-fit prepivoting Bayesian bootstrap
#' @include group.action.R
#' @include randomization.test.R
#' @importFrom stats rexp runif
#' @importFrom compiler cmpfun
#' @examples
#'\dontrun{
#' Y0 <- rnorm(100, 1, 1)
#' Y1 <- rbeta(100,2,2)
#' Tx = sample(100) <= 0.5*(100)
#' # Observed Outcome 
#' Y = ifelse( Tx, Y1, Y0 )
#' dta <- data.frame(Y = Y, A = as.numeric(Tx))
#' pKS.GoF<-prepivot.ks.permtest(dta$Y[dta$A==1],dta$Y[dta$A==0],alpha=0.05,B=1000,n.perm = 999)
#' summary(pKS.GoF)
#' }
#' @export


prepivot.ks.permtest <- function(Y1,Y0, alpha, B, n.perm){
  
  if (anyNA(c(Y1,Y0))) stop("NAs in first or second argument")
  if (!is.numeric(c(Y1,Y0)))  stop("Arguments must be numeric")
  
  # First group's sample size
  m <- length(Y1)
  # Pooled sample
  Z <- c(Y1,Y0)
  N <- length(Z)
  # Calculate permutations
  Sn <- group.action(Z,n.perm,"permutations")
  # Matrix containing data and all permutations of it. First column is the identity permutation
  Z.pi <- cbind(Z,Sn)
  # Observed KS statistic (including the one from actual data)
  K.obs <- apply(Z.pi, 2, function(Z) SKS.stat(Z,m))
  # Pre-pivoted statistics
  bootstrap.samples <-apply(Z.pi, 2, function(z) weighted.bootstrap(z,B))
  w.SKS <- apply(rbind(Z.pi,bootstrap.samples),2, function(W) boot.SKS.stat( matrix(W,ncol=(B+1)), m))
  Tmn   <- apply(rbind(K.obs,w.SKS),2, function(x) mean((x[-1]<=x[1])))
  # Permutation Test
  pt.out  <- unname(randomization.test(Tmn[1],Tmn[-1],alpha))
  # P-value
  pvalue <- mean(Tmn>=Tmn[1])
  
  object_perm<-list() #Generates an empty list to collect all the required info for summary
  object_perm$N<-N
  object_perm$T.obs<-Tmn[1]
  object_perm$cv <- pt.out[2]
  object_perm$pvalue<-pvalue
  object_perm$rejectrule<-pt.out[1]
  object_perm$T.perm<- Tmn[-1]
  object_perm$n_perm<- n.perm
  object_perm$B<- B
  object_perm$sample_sizes<-c(m,N-m)
  
  
  class(object_perm)<-"prepivot.ks.permtest"
  return(object_perm)
}





"SKS.stat" <- function(Z, m){
  # ---------------------------------------------------------------------------------- #
  # Description: calculates the classical KS statistic
  # ---------------------------------------------------------------------------------- #
  # INPUTS: - Y: vector of pooled data from the treatment and control group, respectively
  #         - m: the sample size of the treatment group
  # ---------------------------------------------------------------------------------- #
  # RETURNS:- ks.stat: the numerical calculation of the KS statistic
  # ---------------------------------------------------------------------------------- #  
  N <- length(Z)
  Y1 = Z[1:m]
  Y0 = Z[(m+1):length(Z)]
  
  unique.points = c(Y1, Y0)
  Fn1 = ecdf(Y1)
  Fn0 = ecdf(Y0)
  difference = Fn1(unique.points) - Fn0(unique.points)
  
  ks.stat <- sqrt((m*(N-m))/N)*max(abs(difference))
  return(ks.stat)
}


"pre.boot.SKS.stat" <- function(bootstrap.data,m){
  # ---------------------------------------------------------------------------------- #
  # Description: calculates the KS statistic based on the weighted bootstrap process.
  # ---------------------------------------------------------------------------------- #
  # INPUTS: - bootstrap.data: this is a matrix containing the bootstrap samples. It 
  #           is assumed that the first column is the observed data.
  #         - m: number of observations from the treatment group.
  # ---------------------------------------------------------------------------------- #
  # RETURNS:- w.ks.stat: the numerical calculation of the KS statistic based on the
  #                      bootstrap empirical process.
  # ---------------------------------------------------------------------------------- #
  N <- nrow(bootstrap.data)
  B <- ncol(bootstrap.data)
  # Define the observed data
  Y1 = bootstrap.data[1:m,1]
  Y0 = bootstrap.data[(m+1):N,1]
  # Evaluations of the empirical CDFs are based on unique points from the observed data
  unique.points = c(Y1, Y0)
  # Define the empirical CDFs for observed and bootstrap samples
  Fn1   <-  ecdf(Y1)
  Fn0   <-  ecdf(Y0)
  w.Fn1 <- apply(bootstrap.data[,-1], 2, function(z) w.Fn1=ecdf(z[1:m]))
  w.Fn0 <- apply(bootstrap.data[,-1], 2, function(z) w.Fn0=ecdf(z[(m+1):N]))
  # Evaluate the empirical CDFs at unique points
  Fn1.u    <- Fn1(unique.points)
  Fn0.u    <- Fn0(unique.points)
  w.Fn1.u  <- lapply(1:(B-1), function(b) w.Fn1[[b]](unique.points))
  w.Fn0.u  <- lapply(1:(B-1), function(b) w.Fn0[[b]](unique.points))
  # weighted bootstrap Empirical Process
  boot.process <- mapply(function(X,Y) sqrt((m*(N-m))/N)*( (X-Y)-(Fn1.u - Fn0.u) ),w.Fn1.u, w.Fn0.u)
  # Calculation of the two-sample Kolmogorov-Smirnov statistic
  w.ks.stat <- apply(boot.process,2, function(K) max(abs(K)))
  # Output
  return(w.ks.stat)
}

"boot.SKS.stat" <- compiler::cmpfun(pre.boot.SKS.stat)



"pre.weighted.bootstrap" <- function(Z,B=1000){
  # ---------------------------------------------------------------------------------- #
  # Description: Generate bootstrap samples. 
  # ---------------------------------------------------------------------------------- #
  # INPUTS: - Z: pooled data. This function only considers the outcome of interest as
  #              data input, but it can easily be generalized to be a data.frame. See
  #              for example the documentation of bayes_boot in the "bayesboot" package
  #         - B: number of bootstrap samples
  # ---------------------------------------------------------------------------------- #
  # RETURNS:- w.bootstrap.sample: a M \times B matrix whose columns are the boostrap samples.
  # ---------------------------------------------------------------------------------- # 
  
  # Sample size
  N <- length(Z)
  # Draw from a uniform Dirichlet dist. with alpha set to rep(1, n_dim).
  # Using the facts that Dirichlet draws arises from gamma distribution 
  # and that rgamma(n, 1) is equivalent to rexp(n, 1)
  dirichlet.weights  <- matrix( rexp(N* B, 1) , ncol = N, byrow = TRUE)
  dirichlet.weights  <- dirichlet.weights/rowSums(dirichlet.weights)
  w.bootstrap.sample <- apply(dirichlet.weights, 1, function(w) sample(Z,N,replace = TRUE,prob = w))
  return(w.bootstrap.sample)
}

"weighted.bootstrap" <- compiler::cmpfun(pre.weighted.bootstrap)





