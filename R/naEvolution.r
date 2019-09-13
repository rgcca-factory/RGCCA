#'Analysis of the comparison of different NA methods on RGCCA for increasing percent of missing data in each block
#'
#' @param A A list of complete blocks
#' @param listNAdataset if TRUE, no RGCCA on complete data is run
#' @param C,tau,scale,scheme,sameBlockWeight parameters of RGCCA
#' @param listMethods
#' @param nDatasets Number of simulated dataset
#' @return \item{resultComparison} A list of dataset indicators containg a list of rgcca with a list of criterion. 
#' @return \item{crit} Convergence criterion : abs(1-obj_k/obj_{k-1})
#' @return \item{obj} Vector containing the mean square error between the predict values and the original non missing values at each iteration
#' @title Evolution of fairness of rgcca with increasing missing values
#' @examples 

naEvolution=function(A,prctNA=c(0.1,0.2,0.3),listMethods=c("mean"))
{
  resultComparison=list()
  for(prct in prctNA)
  {
    resultComparison[[as.character(prct)]]=list()
    resultComparison[[as.character(prct)]]=whichNAmethod(A,listMethods=listMethods,patternNA=rep(prct,length(A)))
  }
  return(resultComparison)
}