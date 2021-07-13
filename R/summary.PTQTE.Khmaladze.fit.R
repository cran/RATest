#' Summarizing Quantile-Based Permutation Test with an Estimated Nuisance Parameter
#' 
#' \code{summary} method for class \code{"PTQTE.Khmaladze.fit"}
#' 
#' @method summary PTQTE.Khmaladze.fit
#' @param object an object of class \code{"PTQTE.Khmaladze.fit"}, the result of calling \code{\link{PTQTE.Khmaladze.fit}}
#' @param digits number of digits to display
#' @param ... unused
#' @return \code{PTQTE.Khmaladze.fit} returns an object of \link{class} "\code{PTQTE.Khmaladze.fit}" which has the following components
#'  \item{results}{Matrix with the Testing Problem, Sample Sizes, Number of Permutations, Observed test Statistic, Binary Rule and Significance Level.}
#' @author Maurcio Olivares
#' @export



summary.PTQTE.Khmaladze.fit<-function(object, ..., digits=max(3, getOption("digits") - 3)){
  
  cat("\n")
  cat("**********************************************\n")
  cat("**      Quantile-Based Permutation Test     **\n")   
  cat("**   with an Estimated Nuisance Parameter   **\n")
  cat("***********************************************\n")
  cat("\n")
  cat("* --------------------------------------------- *\n")
  cat("* Testing Problem: testing whether the quantile *\n")
  cat("* treatment effects with an estimated nuisance  *\n")
  cat("* parameter are constant across quantiles       *\n")
  cat("\n")
  cat(paste("Total Number of Observations (Pooled Sample): ", object$N ,sep=""))
  cat("\n")
  z <- as.data.frame(cbind(c("Treatment", "Control"),object$sample_sizes))   
  colnames(z)<-c("Group","Obs")
  print(z,row.names=FALSE)
  cat("\n")
  cat(paste("Number of Permutations: ",object$n_perm,sep=""))
  cat("\n")
  cat("* --------------------------------------------*\n")
  cat("*   H0: constant quantile treatment effects   *\n")
  cat("* --------------------------------------------*\n\n")
  cat(paste("Test Statistic: ", round(object$KS.obs,3),sep=""))
  cat("\n\n")
  if(object$rej.rule==1){
    description<-"Reject H0"
  }else if(object$rej.rule==0){
    description<-"Do not Reject H0"
  }
  cat(paste("Conclusion: ", description ,sep=""))
  cat("\n")
  cat(paste("Significance Level: ", round(object$alpha,3) ,sep=""))
  cat("\n")
  
}