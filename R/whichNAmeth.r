#' whichNAmethod
#' 
#' Analysis of the comparison of different methods to deal with missing data in RGCCA or SGCCA.
#' @inheritParams rgcca
#' @param listMethods vector containing a list of methods ("mean","complete","nipals"...)
#' \itemize{
#' \item{\code{"mean"}}{ corresponds to an imputation by the colmeans}
#' \item{\code{"complete"}}{ corresponds to run RGCCA only on the complete subjects (subjects with missing data are removed)}
#' \item{\code{"nipals"}}{ corresponds to run RGCCA on all available data (NIPALS algorithm)}
#' \item{\code{"em"}}{ corresponds to impute the data with EM-type algorithms}
#' \item{\code{"sem"}}{ corresponds to impute the data with EM-type algorithms with superblock approach}
#' \item{\code{"knn1"}}{ corresponds to impute the data with the 1-Nearest Neighbor. 1 can be replace by another number (such as knn3) to impute with the 3-Nearest Neighbors.}}
#' @param nDatasets Number of simulated datasets
#' @param patternNA pattern of missing values required. If NULL, result of \link{get_patternNA},else, vector containing the percent of missing data for each block if typeNA="ponc" or "block. If type NA="byVar",list of the same size than blocks. Each element of this list should be a vector whose length is the number of variable in the block, containing for each variable the proportion of missing values.
#' @param typeNA structure of missing data required ("ponc" or "block" or "byVar")
#' @param seed if filled (by a number), the randomness is reproducible.
#' @param typeRGCCA type of analysis to be run ("sgcca"or rgcca"...).
#' @return a whichNAmethod object: a list of length nDataset containing. Each element of this list corresponding to each simulated dataset is
#' a list whose names are the chosen missing methods. Each element of such a list is also a list containing
#' \itemize{ 
#' \item{a}{: A vector of size J (number of blocks) corresponding to the norm of difference between the original weight (first axis) and the simulated one for each block }
#'  \item{rv}{: A vector of size J (number of blocks) corresponding to the rv coefficient between the original individual map (first axis) and the simulated one for each block }
#'  \item{bm}{: A vector of size J (number of blocks) corresponding to the percent of different top ten variables correlated with first axis and the simulated one for each block}
#' \item{rvComplete}{: A vector of size J (number of blocks) corresponding to the rv coefficient between the original individual map (first axis) and the simulated one for each block ONLY on the complete individuals}
#' \item{rmse}{: Available only if imputation. A vector of size J (number of blocks) corresponding to the Root Mean Square Error between the original scores and the simulated ones for each block}
#' }
#' @examples 
#' set.seed(42);X1=matrix(rnorm(350),70,5);X2=matrix(rnorm(280),70,4);
#' colnames(X1)=paste("blocks",1:5);colnames(X2)=paste("B",1:4);
#' rownames(X1)=rownames(X2)=paste("S",1:70)
#' A=list(X1,X2);
#' res=whichNAmethod(A,listMethods=c("nipals","mean"),patternNA=rep(0.1,2))
#' @export
#' @importFrom parallel mclapply
#' @seealso \link{plot.whichNAmethod}, \link{naEvolution}
whichNAmethod=function(blocks,listMethods=c("complete","nipals"),typeNA="block",nDatasets=20,patternNA=NULL,connection=matrix(1,length(blocks),length(blocks))-diag(length(blocks)), tau=rep(1,length(blocks)),
                       ncomp=rep(2,length(blocks)),scale_block=TRUE,scale=TRUE,tol=1e-6,
                       verbose=FALSE,scheme="centroid",seed=NULL,typeRGCCA="rgcca",sparsity=NULL)
{
  check_connection(connection,blocks)
  check_tau(tau,blocks)
  check_ncomp(tau,blocks)
  check_integer("nDatasets",nDatasets)
  check_boolean("scale_block",scale_block)
  check_boolean("scale",scale)
  check_boolean("verbose",verbose)
  check_integer("tol",tol,float=TRUE,min=0)
  choices <- c("horst", "factorial", "centroid")
  if (!scheme %in% (choices) && !is.function(scheme))
      stop_rgcca(paste0(scheme, " must be one of ", paste(choices, collapse = ", "), "' or a function."))
  
#  if(length(seed)!=0){check_integer("seed",seed)}
  match.arg(typeNA,c("block","ponc","rand","byVar"))
  match.arg(typeRGCCA,c("rgcca","sgcca"))
  if(is.null(patternNA)){patternNA=get_patternNA(blocks)$pctNAbyBlock}
  if(is.vector(patternNA)){if(length(patternNA)!=length(blocks)){stop_rgcca("patternNA should have the same size as length(blocks)")}}
  referenceDataset=intersection_list(blocks)

  # Getting list of datasets stemming from referenceDataset with the same pattern of missing values
#  if(is.null(listNAdataset))
#  {
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
#  }

  if(typeRGCCA=="rgcca")
  {
      referenceRgcca=rgccad(referenceDataset,C=connection,tau=tau,ncomp=ncomp,verbose=verbose,scale_block=scale_block,scale=scale,tol=tol,scheme=scheme)
      if(verbose)    {  print("comparisons of RGCCA with the different methods...(this could take some time)")}
  }
  if(typeRGCCA=="sgcca")
  {
      referenceRgcca=sgcca(referenceDataset,C=connection,sparsity=sparsity,ncomp=ncomp,verbose=verbose,scale_block=scale_block,scale=scale,tol=tol,scheme=scheme)
      if(verbose)    {  print("comparisons of SGCCA with the different methods...(this could take some time)")}
    
    }
  resultComparison=NULL
  resultComparison=parallel::mclapply(1:nDatasets,function(i)
  {  
    selectCompletePatient=listNAdataset[[i]]$subjectKept
    indicators=NULL
    for(method in listMethods)
    {  
        if(typeRGCCA=="rgcca")
        {
        methodRgcca=rgccaNa(blocks=listNAdataset[[i]]$dat,connection=connection,tau=tau,method=method,ncomp=ncomp,scale_block=scale_block,scale=scale,tol=tol,verbose=FALSE,scheme=scheme)
        }
        if(typeRGCCA=="sgcca")
        {
            methodRgcca=sgccaNa(blocks=listNAdataset[[i]]$dat,connection=connection,sparsity=sparsity,method=method,ncomp=ncomp,scale_block=scale_block,scale=scale,tol=tol,verbose=FALSE,scheme=scheme)
        }
        
        indicators[[method]]=comparison(rgcca1=referenceRgcca,rgcca2=methodRgcca$rgcca,selectPatient=selectCompletePatient,indNA=methodRgcca$indNA)
    }
    return(indicators)
  })
  class(resultComparison)<- "whichNAmethod"
  return(resultComparison)
}
