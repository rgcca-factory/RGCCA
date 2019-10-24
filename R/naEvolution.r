#'Analysis of the comparison of different NA methods on RGCCA for increasing percent of missing data in each block
#' @param A A list of complete blocks
#' @param listNAdataset if TRUE, no RGCCA on complete data is run
#' @param C,tau,scale,scheme,sameBlockWeight parameters of RGCCA
#' @param listMethods Vector containing the names of the methods to use : "complete","nipals","em","knn1","knn10"....
#' @param nDatasets Number of simulated dataset
#' @return \item{resultComparison}{A list of dataset indicators containg a list of rgcca with a list of criterion. }
#' @return \item{crit}{Convergence criterion : abs(1-obj_k/obj_{k-1})}
#' @return \item{obj}{Vector containing the mean square error between the predict values and the original non missing values at each iteration}
#' @title Evolution of fairness of rgcca with increasing missing values
#' @examples 
#' data(Russett)
#' library(FactoMineR)
#' library(parallel)
#' X_agric =as.matrix(Russett[,c("gini","farm","rent")])
#' X_ind = as.matrix(Russett[,c("gnpr","labo")])
#' X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
#' A = list(agri=X_agric, ind=X_ind, polit=X_polit)
#' #ponctual
#' listResults=naEvolution(A=A,listMethods=c("complete","nipals","pca"),
#' prctNA=c(0.05,0.1,0.15,0.2,0.25,0.3,0.4),typeNA="ponc",ncomp=rep(1,3),
#' sameBlockWeight=FALSE)
#' plotEvol(listResults,output="a",barType = "stderr",ylim=c(0,0.2))
#' @export naEvolution
naEvolution=function(A,prctNA=c(0.1,0.2,0.3),listMethods=c("mean"),typeNA="block",ncomp=rep(1,length(A)),sameBlockWeight=TRUE,scale=TRUE,nDatasets=20,tol=1e-6,verbose=FALSE,scheme="centroid")
{
  resultComparison=list()
  for(prct in prctNA)
  { print(paste("pourcent=",prct))
    resultComparison[[as.character(prct)]]=list()
    resultComparison[[as.character(prct)]]=whichNAmethod(A=A,listMethods=listMethods,patternNA=rep(prct,length(A)),typeNA=typeNA,ncomp=ncomp,sameBlockWeight=sameBlockWeight,scale=scale,nDatasets=nDatasets,tol=tol,verbose=verbose,scheme=scheme)
  }
  return(resultComparison)
}