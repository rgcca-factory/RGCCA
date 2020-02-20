#'Analysis of the comparison of different NA methods on RGCCA for increasing percent of missing data in each block
#' @inheritParams whichNAmethod
#' @param prctNA if number, percentage of missing data required for each block. If a vector with the same size as blocks, the percentage can be adapted for each block. If a list of values, this percent is calculated per variable
#' @param seed NULL by default (no reproducibility). A number representing the seed (for reproducibility)
#' @return \item{resultComparison}{A list of dataset indicators containg a list of rgcca with a list of criterion. }
#' @return \item{crit}{Convergence criterion : abs(1-obj_k/obj_{k-1})}
#' @return \item{obj}{Vector containing the mean square error between the predict values and the original non missing values at each iteration}
#' @title Evolution of fairness of rgcca with increasing missing values
#' @examples 
#' data(Russett)
#' library(parallel)
#' X_agric =as.matrix(Russett[,c("gini","farm","rent")])
#' X_ind = as.matrix(Russett[,c("gnpr","labo")])
#' X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
#' A = list(agri=X_agric, ind=X_ind, polit=X_polit)
#' #ponctual
#' listResults=naEvolution(blocks=A,listMethods=c("complete","nipals","mean"),
#' prctNA=c(0.05,0.1,0.15,0.2,0.25,0.3,0.4),typeNA="ponc",ncomp=rep(1,3),
#' sameBlockWeight=FALSE)
#' plotEvol(listResults,output="a",barType = "stderr",ylim=c(0,0.2))
#' @export naEvolution
naEvolution=function(blocks,prctNA=c(0.1,0.2,0.3),listMethods=c("mean"),typeNA="block",ncomp=rep(1,length(blocks)),sameBlockWeight=TRUE,scale=TRUE,nDatasets=20,tol=1e-6,verbose=FALSE,scheme="centroid",seed=NULL,connection=matrix(1,length(blocks),length(blocks))-diag(length(blocks)),tau=rep(1,length(blocks)))
{
     if(any(prctNA>1)){stop("prctNA should be a vector of proportion of missing data (between 0 and 1)")}
    match.arg(typeNA,c("block","ponc","byVar","rand"))
    check_ncomp(ncomp,blocks)
    check_boolean("sameBlockWeight",sameBlockWeight)
    check_boolean("scale",scale)
    check_integer("nDatasets",nDatasets)
    check_boolean("verbose",verbose)
    check_integer("tol",tol,float=TRUE,min=0)
    choices <- c("horst", "factorial", "centroid")
    if (!scheme %in% (choices) && !is.function(scheme))
        stop(paste0(scheme, " must be one of ", paste(choices, collapse = ", "), "' or a function."))
    check_connection(connection,blocks)
    check_tau(tau,blocks)
    if(!is.null(seed)){check_integer("seed",seed)}
    
    resultComparison=list()
    i=0
    for(prct in prctNA)
    {
       if(verbose)
            print(paste("pourcent=",prct))
        
        resultComparison[[as.character(prct)]]=list()
        resultComparison[[as.character(prct)]]=whichNAmethod(blocks=blocks,connection=connection,tau=tau,listMethods=listMethods,patternNA=rep(prct,length(blocks)),typeNA=typeNA,ncomp=ncomp,sameBlockWeight=sameBlockWeight,scale=scale,nDatasets=nDatasets,tol=tol,verbose=verbose,scheme=scheme,seed=seed+i)
        i=i+10
    }
    return(resultComparison)
}