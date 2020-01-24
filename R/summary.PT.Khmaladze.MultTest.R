#' Summarizing Permutation Test for Within-grpup Treatment Effect Heterogeneity in the presence of an Estimated Nuisance Parameter
#' 
#' \code{summary} method for class \code{"PT.Khmaladze.fit"}
#' 
#' @method summary PT.Khmaladze.MultTest
#' @param object an object of class \code{"PT.Khmaladze.MultTest"}, usually a result of a call to \code{\link{PT.Khmaladze.MultTest}}
#' @param digits number of digits to display
#' @param ... unused
#' @return \code{summary.PT.Khmaladze.MultTest} returns an object of \link{class} "\code{summary.PT.Khmaladze.MultTest}" which has the following components
#'  \item{results}{Matrix with the Testing Problem, Number of Permutations for the test and the multiple testing procedure, number of subgroups, (raw) p-values, adjusted p-values, Test Statistic.}
#' @author Maurcio Olivares
#' @export

summary.PT.Khmaladze.MultTest<-function(object, ..., digits=max(3, getOption("digits") - 3)){
  
  cat("\n")
  cat("**************************************************************************\n")
  cat("**   Testing Heterogeneous Treatment Effect with a nuisance parameter:  **\n")   
  cat("**               Within-grpup Treatment Effect Heterogeneity            **\n")
  cat("**************************************************************************\n")
  cat("\n")
  cat("* Testing Problem: testing the joint null hypothesis that the   *\n")
  cat("* treatment outcome distribution is the same as the control     *\n")
  cat("* outcome distribution shifted by the group-specific, constant  *\n")
  cat("* average treatment effect for al mutually exclusive subgroups  *\n")  
  cat("\n")
  cat(paste("Number of Subgroups: ", object$n.subgroups,sep=""))
  cat("\n")
  cat(paste("Number of Permutations for the Permutation Test: ",object$n.perm,sep=""))
  cat("\n")
  cat(paste("Multiple Testing Procedure: ", object$description,sep=""))
  cat("\n")
  if(object$description %in% c("maxT","minP")){
    cat(paste("Number of Permutations for the Westfall-Young procedure: ",object$B,sep=""))
    cat("\n")
  }
  cat("\n")
  cat("* ----------------------------------------------------------------- *\n")
  decision <- ifelse(sum((object$adj.pvalues<object$alpha))>0,"Reject the joint null hypothesis H0","Do not reject the joint null hypothesis H0")
  cat(paste("Decision: ", decision," at level ","\u0251","=",object$alpha ,sep=""))
  cat("\n")
}


