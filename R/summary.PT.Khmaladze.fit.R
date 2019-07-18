#' Summarizing Permutation Test for Heterogeneous Treatment Effects with Estimated Nuisance Parameter
#' 
#' \code{summary} method for class \code{"PT.Khmaladze.fit"}
#' 
#' @method summary PT.Khmaladze.fit
#' @param object an object of class \code{"PT.Khmaladze.fit"}, usually a result of a call to \code{\link{PT.Khmaladze.fit}}
#' @param digits number of digits to display
#' @param ... unused
#' @return \code{summary.PT.Khmaladze.fit} returns an object of \link{class} "\code{summary.PT.Khmaladze.fit}" which has the following components
#'  \item{results}{Matrix with the Testing Problem, Sample Sizes, Number of Permutations, ATE, Test Statistic, Critical value and P-value.}
#' @author Maurcio Olivares Gonzalez
#' @author Ignacio Sarmiento Barbieri
#' @export



summary.PT.Khmaladze.fit<-function(object, ..., digits=max(3, getOption("digits") - 3)){
  
  cat("\n")
  cat("**************************************************************************\n")
  cat("**   Testing Heterogeneous Treatment Effect with a nuisance parameter:  **\n")   
  cat("**                   Asymptotically Robust Permutation Test             **\n")
  cat("**************************************************************************\n")
  cat("\n")
  cat("* Testing Problem: testing goodness of fit in the presence *\n")
  cat("* of a nuisance parameter using a permutation test that is *\n")
  cat("* asymptotically valid under fairly weak assumptions.      *\n")
  cat("\n")
  cat(paste("Total Number of Observations (Pooled Sample): ", object$N ,sep=""))
  cat("\n")
  z <- as.data.frame(cbind(c("Treatment", "Control"),object$sample_sizes))   
  colnames(z)<-c("Group","Obs")
  print(z,row.names=FALSE)
  cat("\n")
  cat(paste("Number of Permutations: ",object$n_perm,sep=""))
  cat("\n\n")
  cat("* --------------------------------------*\n")
  cat(paste("*   H0: F1(y+", "\u03B4",") = F0(y), for some ","\u03B4  *\n"))
  cat("* --------------------------------------*\n\n")
  cat(paste("Estimated Nuisance Parameter (ATE): ", round(object$shift,3), sep=""))
  cat("\n")
  cat(paste("Test Statistic: ", round(object$T.obs,3),sep=""))
  cat("\n")
  cat(paste("Critical Value: "),round(object$cv,3), sep="")
  cat("\n")
  cat(paste("P-Value: ", round(object$pvalue,3) ,sep=""))
  cat("\n")
  
}


