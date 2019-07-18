#' @title  General Construction of Permutation Tests: Group Actions
#'
#' @description Calculates the pre-specified actions on data. Consider data \eqn{Z}{Z} taking values in a sample space \eqn{\Omega}. Let \deqn{\mathbf{G}} be a finite group of transformations from \eqn{\Omega}  onto itself, with \deqn{M=\vert \mathbf{G}\vert}. This function applies \deqn{gZ} as \eqn{g}{g} varies in \deqn{\bf{G}}. If \eqn{Z}{Z} is a vector of size\eqn{N}{N} and the actions \eqn{g}{g} are permutations, \deqn{M=N!}. If the actions \eqn{g}{g} are sign changes, then \deqn{M=\{1,-1\}^{N}}.                
#' 
#' @param Z Numeric. A vector of size \eqn{N}{N} to which the group action will act on. In the two-sample testing problem, \eqn{Z}{Z} is the pooled sample.
#' @param M Numeric. Number of actions to be performed. This is the number of transformations used in the stochastic approximation to the test. This is due to the fact that in some cases \deqn{M=\vert \mathbf{G}\vert} is too large, which makes the application of the actions computationally expensive.
#' @param type Character. The action to be performed. It represents \deqn{gx}, the action the action  of \deqn{g\in\mathbf{G}} on \deqn{x\in\Omega}. It can be either permutations or sign changes.
#' @return Numeric. A matrix of size \deqn{N\times M} where \eqn{N}{N} is the size of input \eqn{Z}{Z} and \eqn{M}{M} is the number of actions to be performed on \eqn{Z}{Z}.
#'
#' @author Maurcio Olivares Gonzalez
#' @author Ignacio Sarmiento Barbieri
#' @references
#' Lehmann, Erich L. and Romano, Joseph P (2005) Testing statistical hypotheses.Springer Science & Business Media.
#' @keywords permutation test rdperm group action
#' @export


group.action <- function(Z,M,type="permutations"){
  if(!(type%in%c("permutations","sign.changes"))){
    print("Must specify group action. Options include 'permutations' or sign.changes'' ")
    stop()
  }
  
  if(type=="permutations"){
    # Generate random permutations of data
    sample.indexes = lapply(1:M, function(x) sample(1:length(Z)) )
    S_action_list<-lapply(sample.indexes, function(x,db){db[x]},Z)
    
  } else if(type=="sign.changes"){
    # Generate random permutations of {1,-1} (sign changes)
    sample.indexes = lapply(1:M, function(x) sample(c(1,-1),length(Z),replace = TRUE) )
    S_action_list<-lapply(sample.indexes, function(x,y){y*x},Z)
  }
  Sn<-do.call(cbind,S_action_list)
  return(Sn)
}