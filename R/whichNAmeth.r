#' Analysis of the comparison of different NA methods on RGCCA
#' @param listNAdataset if TRUE, no RGCCA on complete data is run
#' @inheritParams rgcca
#' @param listMethods vector containing a list of methods ("mean","complete","nipals"...)
#' @param nDatasets Number of simulated dataset
#' @param patternNA number of missing data required
#' @param typeNA structure of missing data required ("ponc" or "block")
#' @param seed if filled (by a number), the randomness is reproducible.
#' @param typeRGCCA type of RGCCA ("sgcca","rgcca"...).
#' @return \item{A}{A list of dataset indicators containg a list of rgcca with a  list of criterion. }
#' @return \item{crit}{Convergence criterion : abs(1-obj_k/obj_{k-1})}
#' @return \item{obj}{Vector containing the mean square error between the predicted values and the original non missing values at each iteration}
#' @title comparison of two RGCCA results
#' @examples 
#' set.seed(42);X1=matrix(rnorm(350),70,5);X2=matrix(rnorm(280),70,4);
#' colnames(X1)=paste("blocks",1:5);colnames(X2)=paste("B",1:4);
#' rownames(X1)=rownames(X2)=paste("S",1:70)
#' A=list(X1,X2);
#' res=whichNAmethod(A,listMethods=c("nipals","mean"),patternNA=rep(0.1,2))
#' @export
#' @importFrom parallel mclapply

whichNAmethod=function(blocks,listNAdataset=NULL,connection=matrix(1,length(blocks),length(blocks))-diag(length(blocks)), tau=rep(1,length(blocks)),
                       listMethods,nDatasets=20,patternNA=NULL,typeNA="block",ncomp=rep(2,length(blocks)),sameBlockWeight=TRUE,scale=TRUE,tol=1e-6,
                       verbose=FALSE,scheme="centroid",seed=NULL,typeRGCCA="rgcca",sparsity=NULL)
{
  check_connection(connection,blocks)
  check_tau(tau,blocks)
  check_ncomp(tau,blocks)
  check_integer("nDatasets",nDatasets)
  check_boolean("sameBlockWeight",sameBlockWeight)
  check_boolean("scale",scale)
  check_boolean("verbose",verbose)
  check_integer("tol",tol,float=TRUE,min=0)
  choices <- c("horst", "factorial", "centroid")
  if (!scheme %in% (choices) && !is.function(scheme))
      stop(paste0(scheme, " must be one of ", paste(choices, collapse = ", "), "' or a function."))
  
#  if(length(seed)!=0){check_integer("seed",seed)}
  match.arg(typeNA,c("block","ponc","rand","byVar"))
  match.arg(typeRGCCA,c("rgcca","sgcca"))
  if(is.null(patternNA)){patternNA=determine_patternNA(blocks,graph=FALSE)$pctNAbyBlock}
  if(is.vector(patternNA)){if(length(patternNA)!=length(blocks)){stop("patternNA should have the same size as length(blocks)")}}
  referenceDataset=intersection(blocks)

  # Getting list of datasets stemming from referenceDataset with the same pattern of missing values
  if(is.null(listNAdataset))
  {
    if(verbose){
        print("creation of datasets with NA...")}
    listNAdataset=lapply(1:nDatasets,function(i)
        {
            if(!is.null(seed)&&length(seed)==nDatasets)
                {   createNA(blocks=referenceDataset,typeNA=typeNA,pNA=patternNA,nAllRespondants=10,output="list",seed=seed[i])}
            else
                {   createNA(blocks=referenceDataset,typeNA=typeNA,pNA=patternNA,nAllRespondants=10,output="list",seed=NULL)}
        }
    )
  }

  if(typeRGCCA=="rgcca")
  {
      referenceRgcca=rgccad(referenceDataset,C=connection,tau=tau,ncomp=ncomp,verbose=verbose,sameBlockWeight=sameBlockWeight,scale=scale,tol=tol,scheme=scheme)
      if(verbose)    {  print("comparisons of RGCCA with the different methods...(this could take some time)")}
  }
  if(typeRGCCA=="sgcca")
  {
      referenceRgcca=sgcca(referenceDataset,C=connection,sparsity=sparsity,ncomp=ncomp,verbose=verbose,sameBlockWeight=sameBlockWeight,scale=scale,tol=tol,scheme=scheme)
      if(verbose)    {  print("comparisons of SGCCA with the different methods...(this could take some time)")}
  }
  resultComparison=NULL
  resultComparison=mclapply(1:nDatasets,function(i)
  {  
    selectCompletePatient=listNAdataset[[i]]$subjectKept
    indicators=NULL
    for(method in listMethods)
    {  
        if(typeRGCCA=="rgcca")
        {
        methodRgcca=rgccaNa(blocks=listNAdataset[[i]]$dat,connection=connection,tau=tau,method=method,ncomp=ncomp,sameBlockWeight=sameBlockWeight,scale=scale,tol=tol,verbose=FALSE,scheme=scheme)
        }
        if(typeRGCCA=="sgcca")
        {
            methodRgcca=sgccaNa(blocks=listNAdataset[[i]]$dat,connection=connection,sparsity=sparsity,method=method,ncomp=ncomp,sameBlockWeight=sameBlockWeight,scale=scale,tol=tol,verbose=FALSE,scheme=scheme)
        }
        indicators[[method]]=comparison(rgcca1=referenceRgcca,rgcca2=methodRgcca$rgcca,selectPatient=selectCompletePatient,indNA=methodRgcca$indNA)
    }
    return(indicators)
  })
  class(resultComparison)<- "whichNAmethod"
  return(resultComparison)
}
