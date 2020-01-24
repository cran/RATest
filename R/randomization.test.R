#' @title  General Construction of Randomization Tests
#'
#' @description Calculates the randomization test. Further discussion can be found in chapter 15 of Lehmann and Romano (2005, p 633). Consider data \eqn{X}{X} taking values in a sample space \eqn{\Omega}. Let \eqn{\mathbf{G}} be a finite group of transformations from \eqn{\Omega}  onto itself, with \eqn{M=\vert \mathbf{G}\vert}. Let \eqn{T(X)} be a real-valued test statistic such that large values provide evidence against the null hypothesis. Denote by \deqn{T^{(1)}(X)\le T^{(2)}(X)\le\dots\le T^{(M)}(X)} the ordered values of \eqn{\{T(gX)\,:\,g\in\mathbf{G}\}}. Let \eqn{k=M-\lfloor M\alpha\rfloor} and define \deqn{M^{+}(x)} and \deqn{M^{0}(x)} be the number of values \eqn{T^{(j)}(X)}, \eqn{j=1,\dots,M}, which are greater than \eqn{T^{(k)}(X)} and equal to \eqn{T^{(k)}(X)} respectively. Set \deqn{a(X)=\frac{\alpha M-M^{+}(X)}{M^{0}(X)}}. The randomization test is given by \deqn{\phi(X)=1} if \eqn{T(x)> T^{(k)}(X)}, \deqn{\phi(X)=0} if \eqn{T(X)< T^{(k)}(X)}, and \deqn{\phi(X)=a(X)} if \eqn{T(X)= T^{(k)}(X)}.  
#' 
#' @param Tn Numeric. A scalar representing the observed test statistic \eqn{T(X)}.
#' @param Tng Numeric. A vector containing \eqn{\{T(gX)\,:\,g\in\mathbf{G}\}}.
#' @param alpha Numeric. Nominal level for the test. The default is 0.05.
#' @return Numeric. A scalar \eqn{\phi(X)\in\{0,1\}}. The test rejects the null hypothesis if \eqn{\phi(X)=1}, and does not reject otherwise.
#'
#' @author Maurcio Olivares Gonzalez
#' @author Ignacio Sarmiento Barbieri
#' @references
#' Lehmann, Erich L. and Romano, Joseph P (2005) Testing statistical hypotheses.Springer Science & Business Media.
#' @keywords permutation test rdperm randomization test
#' @export

randomization.test <- function(Tn,Tng,alpha=0.05){
  x <- c(Tn,Tng)
  M <- length(x)
  y <- sort(x)
  k <- M - floor(M*alpha)
  cv = y[k]
  M.plus <- sum(y>cv)
  M.zero <- sum(y==cv)
  a <-  (M*alpha - M.plus)/M.zero
  phi <- (Tn>cv) + (Tn==cv)*(runif(1)<=a)
  return(phi)
}






