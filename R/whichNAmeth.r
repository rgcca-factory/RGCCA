#' Analysis of the comparison of different NA methods on RGCCA
#' @param A A list of complete blocks
#' @param listNAdataset if TRUE, no RGCCA on complete data is run
#' @param C,tau,scale,scheme,sameBlockWeight,ncomp,tol,verbose parameters of RGCCA
#' @param listMethods vector containing a list of methods ("mean","complete","nipals"...)
#' @param nDatasets Number of simulated dataset
#' @param patternNA number of missing data required
#' @param typeNA structure of missing data required ("ponc" or "block")
#' @param seed if filled (by a number), the randomness is reproducible.
#' @return \item{A}{A list of dataset indicators containg a list of rgcca with a  list of criterion. }
#' @return \item{crit}{Convergence criterion : abs(1-obj_k/obj_{k-1})}
#' @return \item{obj}{Vector containing the mean square error between the predicted values and the original non missing values at each iteration}
#' @title comparison of two RGCCA results
#' @examples 
#' set.seed(42);X1=matrix(rnorm(350),70,5);X2=matrix(rnorm(280),70,4);
#' colnames(X1)=paste("A",1:5);colnames(X2)=paste("B",1:4);
#' rownames(X1)=rownames(X2)=paste("S",1:70)
#' A=list(X1,X2);
#' res=whichNAmethod(A,listMethods=c("nipals","mean"),patternNA=rep(0.1,2))
#' @export
#' @importFrom parallel mclapply

whichNAmethod=function(A,listNAdataset=NULL,C=matrix(1,length(A),length(A))-diag(length(A)), tau=rep(1,length(A)),listMethods,nDatasets=20,patternNA=NULL,typeNA="block",ncomp=rep(2,length(A)),sameBlockWeight=TRUE,scale=TRUE,tol=1e-6,verbose=FALSE,scheme="centroid",seed=NULL)
{
  if(is.null(patternNA)){patternNA=determine_patternNA(A,graph=FALSE)$pctNAbyBlock}
  if(is.vector(patternNA)){if(length(patternNA)!=length(A)){stop("patternNA should have the same size as length(A)")}}
  referenceDataset=intersection(A)

  # Getting list of datasets stemming from referenceDataset with the same pattern of missing values
  if(is.null(listNAdataset))
  {
    print("creation of datasets with NA...")
    listNAdataset=lapply(1:nDatasets,function(i)
        {
            if(!is.null(seed)&&length(seed)==nDatasets)
                {   createNA(A=referenceDataset,option=typeNA,pNA=patternNA,nAllRespondants=10,output="list",seed=seed[i])}
            else
                {   createNA(A=referenceDataset,option=typeNA,pNA=patternNA,nAllRespondants=10,output="list",seed=NULL)}
        }
    )
  }

  referenceRgcca=rgcca(referenceDataset,C=C,tau=tau,ncomp=ncomp,verbose=verbose,sameBlockWeight=sameBlockWeight,scale=scale,tol=tol,scheme=scheme)
  print("comparisons of RGCCA with the different methods...(this could take some time)")
  resultComparison=NULL
  resultComparison=mclapply(1:nDatasets,function(i)
  {  
    selectCompletePatient=listNAdataset[[i]]$subjectKept
    indicators=NULL
    for(method in listMethods)
    {  
        methodRgcca=rgccaNa(A=listNAdataset[[i]]$dat,C=C,tau=tau,method=method,ncomp=ncomp,sameBlockWeight=sameBlockWeight,scale=scale,tol=tol,verbose=verbose,scheme=scheme)
       indicators[[method]]=comparison(rgcca1=referenceRgcca,rgcca2=methodRgcca$rgcca,selectPatient=selectCompletePatient,indNA=methodRgcca$indNA)
    }
    return(indicators)
  })
  return(resultComparison)
}
