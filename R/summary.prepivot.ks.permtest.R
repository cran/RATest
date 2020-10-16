#' Summarizing Two-sample Goodness-of-fit Permutation Test under Covariate-adaptive Randomization
#' 
#' \code{summary} method for class \code{"prepivot.ks.permtest"}
#' 
#' @method summary prepivot.ks.permtest
#' @param object an object of class \code{"prepivot.ks.permtest"}, the result of calling \code{\link{prepivot.ks.permtest}}
#' @param digits number of digits to display
#' @param ... unused
#' @return \code{summary.prepivot.ks.permtest} returns an object of \link{class} "\code{summary.prepivot.ks.permtest}" which has the following components
#'  \item{results}{Matrix with the Testing Problem, Sample Sizes, Number of Permutations, Number of Bootstrap samples, Observed test Statistic, Critical value and P-value.}
#' @author Maurcio Olivares
#' @export



summary.prepivot.ks.permtest<-function(object, ..., digits=max(3, getOption("digits") - 3)){
  
  cat("\n")
  cat("***********************************************\n")
  cat("**     Two-sample Goodness-of-fit testing    **\n")   
  cat("**   under Covariate-adaptive Randomization  **\n")
  cat("***********************************************\n")
  cat("\n")
  cat("* ---------------------------------------------------------- *\n")
  cat("* Testing Problem: testing equality of two distributions     *\n")
  cat("* under covariate-adaptive randomization using a permutation *\n")
  cat("* test based on a prepivoted Kolmogorov-Smirnov statistic    *\n")
  cat("\n")
  cat(paste("Total Number of Observations (Pooled Sample): ", object$N ,sep=""))
  cat("\n")
  z <- as.data.frame(cbind(c("Treatment", "Control"),object$sample_sizes))   
  colnames(z)<-c("Group","Obs")
  print(z,row.names=FALSE)
  cat("\n")
  cat(paste("Number of Permutations: ",object$n_perm,sep=""))
  cat("\n")
  cat(paste("Number of Bootstrap samples: ",object$B,sep=""))
  cat("\n\n")
  cat("* ----------------*\n")
  cat("*   H0: F1 = F0   *\n")
  cat("* ----------------*\n\n")
  cat(paste("Test Statistic: ", round(object$T.obs,3),sep=""))
  cat("\n")
  cat(paste("Critical Value: "),round(object$cv,3), sep="")
  cat("\n")
  cat(paste("P-Value: ", round(object$pvalue,3) ,sep=""))
  cat("\n")
  
}