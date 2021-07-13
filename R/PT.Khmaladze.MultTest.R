
#' @title Martingale transformed Permutation Test: Multiple Testing procedures.
#'
#' @description This function applies the martingale transformed Permutation test (Chung and Olivares (2020)) to test whether there exists within-group treatment effect heterogeneity. 
#' The method jointly tests the null hypotheses that treatment effects are constant within mutually exclusive subgroups while allowing them to be different across subgroups. More formally, 
#' assume the mutually exclusive subgroups are formed from observed covariates, and are taken as given. Denote \eqn{\mathcal{J}} the total number of such subgroups. 
#' Let \eqn{F_0^{j}(y)} and \eqn{F_1^{j}(y)} be the CDFs of the control and treatment group for subgroup \eqn{1\le j\le \mathcal{J}}. The null hypothesis of interest is given by the joint hypothesis
#' \deqn{\mathbf{H}_{0}: F_1^{j}(y + \delta_{j}) = F_0^{j}(y)}
#' for all mutually exclusive \eqn{j\in\{1,\dots,\mathcal{J}\}}, for some \eqn{\delta_j}. We are treating \eqn{\mathbf{H}_0} as a multiple testing problem in which every individual hypothesis \eqn{j\in\{1,\dots,\mathcal{J}\}}, given by
#' \deqn{H_{0,j}: F_1^{j}(y + \delta_{j}) = F_0^{j}(y)}
#' for some \eqn{\delta_j} specifies whether the treatment effect is heterogeneous for a particular subgroup.
#'
#' To achieve control of the family-wise error rate, the function considers several multiple testing procedures, such as Bonferroni, maxT and minP (Westfall and Young (1993)), and Holm (1979).
#' For further details, see Chung and Olivares (2020). 
#'                                                                                        
#'
#' @param data List. Data are presented in the form of a list, where each sublist contains the treatment and control group observations for a specific subgroup.
#' @param procedure multiple testing procedure. Several options are available, including maxT and minP (Westfall and Young (1993)), Bonferroni adjustment, and Holm (1979) procedure. The default is Bonferroni.
#' @param alpha Significance level.
#' @param n.perm Numeric. Number of permutations needed for the stochastic approximation of the p-values. See Remark 4 in Chung and Olivares (2020). The default is n.perm=499.	
#' @param B Numeric. Number of permutations needed for the stochastic approximation in the Westfall-Young procedures. See Remark 11 in Chung and Olivares (2020). The default is B=499.
#' @param na.action a function to filter missing data. This is applied to the model.frame . The default is na.omit, which deletes observations that contain one or more missing values.

#' @return An object of class "PT.Khmaladze.MultTest" is a list containing at least the following components:


#'  \item{description}{Type of multiple testing adjustment. It can be Westfall-Young's maxT, minP, Holm or Bonferroni.}
#'  \item{n.subgroups}{Number of subgrups for a specific covariate.}
#'  \item{T.obs}{Vector. Observed test statistic for each subgroup.}
#'  \item{pvalues}{Vector. P-value for each individual test.}
#'  \item{adj.pvalue}{Vector. Adjusted p-values according to the user-chosen multiple testing procedure.}
#'  \item{n.perm}{Number of permutations.}
#'  \item{B}{Number of permutations used in the Westfall-Young procedure.}
#'  \item{sample.sizes}{Subgroup sample sizes.}
#'  \item{alpha}{Significance level.}

#' @author Maurcio Olivares
#' @references 
#' Chung, E. and Olivares, M. (2021). Permutation Test for Heterogeneous Treatment Effects with a Nuisance Parameter. Forthcoming in Journal of Econometrics.
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics, pages 65-70.
#' Westfall, P.H. and Young, S.S. (1993). Resampling-based multiple testing: Examples and methods for p-value adjustment, Volume 279, John & Wiley Sons.
#' @keywords Khmaladze Permutation Test Multiple Testing Westfall-Young
#' @include PT.Khmaladze.fit.R
#' @import quantreg
#' @importFrom stats p.adjust
#' @examples
#'\dontrun{
#' subgroup1 <- list()
#' subgroup1$Y0 <- rnorm(11)
#' subgrpup1$Y1 <- rnorm(8,1,1) 
#' subgroup2 <- list()
#' subgroup2$Y0 <- rnorm(9)
#' subgroup2$Y1 <- rnorm(7,1,2)
#' data <- list(subgroup1,subgroup2)
#' res.minP <- PT.Khmaladze.MultTest(data,"minP",n.perm=100,B=100)
#' summary(res.minP)
#' adjusted.p.values <- res.minP$adj.pvalues
#' adjusted.p.values
#' }
#' @export



PT.Khmaladze.MultTest <- function(data,procedure="maxT",alpha=0.05, n.perm=499,B=499, na.action){
  if(!(procedure%in%c("maxT","minP","Holm","Bonferroni"))){
    print("Must specify multiple testing procedure. Options include 'maxT', 'minP', 'Holm', and 'Bonferroni'")
    stop()
  }
  # Calculate the individual permutation tests for every subgroup (sublist).
  individual.Test <- lapply(data, function(s) PT.Khmaladze.fit(s[[1]],s[[2]],alpha=0.05,n.perm=n.perm))
  T.obs <- unlist(lapply(individual.Test, function(s) s$T.obs))
  ind.pvalues <- unlist(lapply(individual.Test, function(s) s$pvalue)) 
  # Composition of the data set
  n.subgroups <- length(data)
  sample.sizes <- unlist(lapply(individual.Test, function(s) s$N))

  # --------------------------------------------------------------------------------- #
  #                           Multiple Testing Procedures                             #
  # --------------------------------------------------------------------------------- # 
  
  # ---------------------------- # 
  #      Westfall-Young maxT     #
  # ---------------------------- # 
  if(procedure=="maxT"){
    T.ord <- sort(T.obs,decreasing = TRUE,index.return=TRUE)$x #max T
    H0.ord<- sort(T.obs,decreasing = TRUE,index.return=TRUE)$ix
    # Stack observations from each subgroup j into a vector of size m_j+n_j
    Z.data <- lapply(data,function(d) unname(unlist(d)))
    # Calculate m_j and n_j for every subgroup j 
    subgroup.samples <- lapply(data, function(d) c(length(d[[1]]),length(d[[2]])))
    # Check you have more observations in every subgroup than permutations B
    B <- ifelse(B>factorial(min(sample.sizes)),factorial(min(sample.sizes)),B)
    # Generate index permutations for each subgroup
    Sn <-  lapply(1:n.subgroups,function(i) lapply(1:B, function(j) sample(1:sample.sizes[i])))
    Pi <- lapply(Sn, function(s) do.call(cbind,s))
    # Generate permutations of the data. Every element "j" of the list is a matrix whose
    # columns represent the permuted samples of subgroup "j"
    Z.pi <- mapply(function(x,y,z) matrix(x[y],nrow = z),Z.data,Pi,sample.sizes)
    # Calculate the martingale transformed KS statistic for each permuted sample and each subgroup.
    # The resulting matrix is of dimension B x j
    KS.pi <- mapply(function(x,y) apply(x,2,function(z) PT.Khmaladze.fit(z[1:y[1]],z[(y[1]+1):(y[1]+y[2])],alpha = alpha, n.perm = n.perm)$T.obs ),Z.pi,subgroup.samples)
    # Define Khat: this is step 1.ii) in Algorithm 2 in the main manuscript.
    KS.pi.hat <- t( apply(KS.pi,1,function(x) sort(x,decreasing = TRUE)))
    # Define H_j as in Step 2 in manuscript.
    H <- unlist(lapply(1:n.subgroups, function(s) sum(T.ord[s]<=KS.pi.hat[,s])))/B
    # Adjusted p-values
    adj.pvalues <- rep(NA,n.subgroups)
    adj.pvalues[1] <- H[1]
    for (r in 2:n.subgroups) {
      adj.pvalues[r] <- max(adj.pvalues[r-1],H[r])
    }
    adj.pvalues <- adj.pvalues[order(H0.ord)]
    # ---------------------------- # 
    #      Westfall-Young minP     #
    # ---------------------------- # 
  } else if(procedure=="minP"){
    pval.ord <- sort(ind.pvalues,decreasing = FALSE,index.return=TRUE)$x #max T
    H0.ord<- sort(ind.pvalues,decreasing = FALSE,index.return=TRUE)$ix
    # Stack observations from each subgroup j into a vector of size m_j+n_j
    Z.data <- lapply(data,function(d) unname(unlist(d)))
    # Calculate m_j and n_j for every subgroup j 
    subgroup.samples <- lapply(data, function(d) c(length(d[[1]]),length(d[[2]])))
    # Check you have more observations in every subgroup than permutations B
    B <- ifelse(B>factorial(min(sample.sizes)),factorial(min(sample.sizes)),B)
    # Generate index permutations for each subgroup
    Sn <-  lapply(1:n.subgroups,function(i) lapply(1:B, function(j) sample(1:sample.sizes[i])))
    Pi <- lapply(Sn, function(s) do.call(cbind,s))
    # Generate permutations of the data. Every element "j" of the list is a matrix whose
    # columns represent the permuted samples of subgroup "j"
    Z.pi <- mapply(function(x,y,z) matrix(x[y],nrow = z),Z.data,Pi,sample.sizes)
    # Calculate the pvalue from the permutation test based on martingale transformed KS statistic 
    # for each permuted sample and each subgroup. The resulting matrix is of dimension B x j
    pval.pi <- mapply(function(x,y) apply(x,2,function(z) PT.Khmaladze.fit(z[1:y[1]],z[(y[1]+1):(y[1]+y[2])],alpha = alpha, n.perm = n.perm)$pvalue ),Z.pi,subgroup.samples)
    # Define p.hat: this is step 1.ii) in Algoritm 1 of the main manuscript.
    pval.pi.hat <- t( apply(pval.pi,1,sort))
    # Define L_j as in Step 2 in manuscript.
    L <- unlist(lapply(1:n.subgroups, function(s) sum(pval.ord[s]>=pval.pi.hat[,s])))/B
    # Adjusted p-values
    adj.pvalues <- rep(NA,n.subgroups)
    adj.pvalues[1] <- L[1]
    for (r in 2:n.subgroups) {
      adj.pvalues[r] <- max(adj.pvalues[r-1],L[r])
    }
    adj.pvalues <- adj.pvalues[order(H0.ord)]
    
    
  } else if(procedure=="Holm"){
    adj.pvalues <- p.adjust(ind.pvalues,method = "holm")
  } else if(procedure=="Bonferroni"){
    adj.pvalues <- ifelse(ind.pvalues*n.subgroups>1,1.0,ind.pvalues*n.subgroups) 
  } 
  
  
  object_perm<-list() #Generates an empty list to collect all the required info for summary
  object_perm$description<-procedure
  object_perm$n.subgroups<-n.subgroups
  object_perm$T.obs<-T.obs
  object_perm$pvalues<-ind.pvalues
  object_perm$adj.pvalues<-adj.pvalues
  object_perm$n.perm<- n.perm
  if(procedure=="maxT"|procedure=="minP") object_perm$B<-B
  object_perm$sample.sizes<-sample.sizes
  object_perm$alpha <- alpha

  class(object_perm)<-"PT.Khmaladze.MultTest"
  return(object_perm)
  
}














