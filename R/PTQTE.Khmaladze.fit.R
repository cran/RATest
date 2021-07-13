#' @title  Quantile-Based Permutation Test with an Estimated Nuisance Parameter
#'
#' @description A permutation test for testing whether the quantile treatment effects are constant across quantiles. The permutation test considered here is based on the Khmaladze transformation of the quantile process (Koenler and Xiao (2002)), and adapted by Chung and Olivares (2021).
#' 
#' @param Y Numeric. Vector of responses.
#' @param Z Numeric. Treatment indicator. Z=1 if the unit is in the treatment group, and Z=0 if the unit is in the control group.
#' @param taus quantiles at which the process is to be evaluated, if any of the taus lie outside (0,1) then the full process is computed for all distinct solutions.
#' @param alpha Significance level.
#' @param n.perm Numeric. Number of permutations needed for the stochastic approximation of the p-values. The default is n.perm=999.

#' @return An object of class "PTQTE.Khmaladze" containing at least the following components:
#' 
#'  \item{n_populations}{Number of grups.}
#'  \item{N}{Sample Size.}
#'  \item{KS.obs}{Observed two-sample Kolmogorov-Smirnov test statistic based on the quantile process.}
#'  \item{shift}{The estimated nuisance parameter.}
#'  \item{rej.rule}{Binary decision for the permutation test, where 1 means rejection.}
#'  \item{pvalue}{P-value.}
#'  \item{KS.perm}{Vector. Test statistic recalculated for all permutations used in the stochastic approximation.}
#'  \item{n_perm}{Number of permutations.}
#'  \item{sample_sizes}{Groups size.}
#'  
#' @author Maurcio Olivares
#' @references  
#' Khmaladze, E. (1981). Martingale Approach in the Theory of Goodness-of-fit Tests. Theory of Probability and its Application, 26: 240–257.
#' Koenker, R. and Xiao, Z. (2002) Inference on the Quantile Regression Process. Econometrica, 70(4): 1583-1612.
#' Chung, E. and Olivares, M. (2021). Comment on "Can Variation in Subgroups' Average Treatment Effects Explain Treatment Effect Heterogeneity? Evidence from a Social Experiment."
#' 
#' @keywords Permutation Test Khmaladze Transformation Quantile Process
#' @include group.action.R
#' @include randomization.test.R
#' @import quantreg
#' @importFrom stats runif
#' @examples
#'\dontrun{
#' dta <- data.frame(Y=rnorm(100),Z=sample(c(0,1), 100, replace = TRUE))
#' pt.QTE<-PTQTE.Khmaladze.fit(dta$Y,dta$Z,taus=seq(.1,.9,by=0.05),alpha=0.05,n.perm = 499)
#' summary(pt.QTE)
#' }
#' @export



PTQTE.Khmaladze.fit <- function(Y,Z,taus=seq(.1,.9,by=0.05), alpha=0.05, n.perm=999){
  
  # (Pooled) Sample size and subgroup sample size
  N<-length(Y);
  lengths <- c(sum(Z),N-sum(Z));
  
  # ---------------- #
  #   Permutations   #
  
  # Generate random permutations of the data
  Z.pi  <- group.action(Z,n.perm,"permutations")
  # Stack all permuted indeces and data. The first column is
  # the original/observed data, second is the group indicator, 
  # rest are permutations of indeces.
  mf<-cbind(Y,Z.pi); 

  
  # ---------------------------- #
  #   Calculates the statistic   #
  # ---------------------------- #
  
  te.hat <- mean(Y[Z==1])-mean(Y[Z==0]);
  Y.star <- Y - Z*te.hat;
  
  
  # Calculation of the test statistic for all permutations. 
  stat <- apply( mf[,-1], 2, function(X) KhmaladzeTest(Y.star~X,taus=taus)$Tn )
  
  # Observed Test statistic
  KS.obs<-stat[1] 
  #Test statistic for the permutated samples
  KS.perm<-stat[-1] 
  
  # Critical Value
  rej.rule <- randomization.test(KS.obs,KS.perm,alpha)
  # P-value
  pval <- mean(ifelse(KS.perm>=KS.obs,1,0))
  
  #Generates an empty list to collect all the required info for summary
  object_perm<-list() 
  object_perm$N<-N
  object_perm$KS.obs<-KS.obs
  object_perm$shift <- te.hat
  object_perm$rej.rule <- rej.rule[1]
  object_perm$pvalue<-pval
  object_perm$T.perm<- KS.perm
  object_perm$n_perm<- n.perm
  object_perm$alpha<- alpha
  object_perm$sample_sizes<-c(lengths)
  
  
  class(object_perm)<-"PTQTE.Khmaladze.fit"
  
  return(object_perm)
  
}