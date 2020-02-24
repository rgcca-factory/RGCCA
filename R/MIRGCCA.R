#'Multiple imputation for RGCCA
#'
#' This method allows multiple imputation for RGCCA with several options.
#' @param ni number of imputations
#' @param k Integer representing the number of neighbors or "auto" or "all"
#' @inheritParams rgcca
#' @param naxis number of component to select in the superblock for the estimation of missing data for "em" option
#' @param superblock TRUE if the A list should be returned in the output, FALSE ifelse
#' @param klim TRUE if the A list should be returned in the output, FALSE ifelse
#' @param output TRUE if the A list should be returned in the output, FALSE ifelse
#' @param option "knn" for k Nearest Neigbors or "em" for Expectation Maximization
#' @return \item{rgcca0}{RGCCA results for the reference dataset}
#' @return \item{data}{list of imputed data obtained}
#' @return \item{rgccaList}{list of RGCCA obtained}
#' @title MIRGCCA: Multiple imputation for RGCCA
#' @examples 
#' set.seed(42);X1=matrix(rnorm(500),100,5);
#' set.seed(22);X2=matrix(rnorm(400),100,4);
#' set.seed(2);X3=matrix(rnorm(700),100,7);
#' rownames(X1)=rownames(X2)=rownames(X3)=paste("S",1:100,sep="")
#' colnames(X1)=paste("A",1:5)
#' colnames(X2)=paste("B",1:4)
#' colnames(X3)=paste("C",1:7)
#' X1[1,]=NA
#' X2[7,1]=NA
#'  X2[5,1]=NA
#'  X3[3,1:2]=NA
#'  A=list(X1,X2,X3)
#' res=MIRGCCA(A,k=3,ni=5,scale=TRUE,sameBlockWeight=TRUE,tau=rep(0,3))
#' @export
MIRGCCA=function(blocks,option="knn",superblock=TRUE,k=5,ni=5,scale=TRUE,sameBlockWeight=TRUE,tau=rep(1:length(A)),klim=NULL,output="mean",scheme="centroid",tol=1e-8,connection=NULL,ncomp=rep(2,length(A)),naxis=1)
{
    A=blocks
    C=connection
    match.arg(option,c("knn","em"))
    check_boolean("superblock",superblock)
    check_boolean("sameBlockWeight",sameBlockWeight)
    check_tau(tau,A)
    check_integer("tol",tol,float=TRUE,min=0)
    check_integer("naxis",naxis)
    check_boolean("scale",scale)
    check_ncomp(ncomp,A)
    check_integer("ni",ni)
    choices <- c("horst", "factorial", "centroid")
    if (!scheme %in% (choices) && !is.function(scheme))
        stop(paste0(scheme, " must be one of ", paste(choices, collapse = ", "), "' or a function."))
    
     if(option=="knn")
  {
    dataTest0=imputeNN(A=A,output=output,k=k,klim=klim)
    if(!is.null(dataTest0))
    {
      rgcca0=rgccad(dataTest0,ncomp=rep(2,length(A)),scale=scale,sameBlockWeight=sameBlockWeight,tau=tau,verbose=FALSE,scheme=scheme,tol=tol)
      #plotRGCCA2(rgcca0,indnames=TRUE,varnames=TRUE)
      dataTest=resRgcca2=resprocrustes=list()
      for(i in 1:ni)
      {
        dataTest[[i]]=imputeNN(A=A,output="random",k=k,klim=klim)
        resRgcca2[[i]]=rgccad(dataTest[[i]],ncomp=rep(2,length(dataTest[[i]])),scale=scale,sameBlockWeight=sameBlockWeight,tau=tau,verbose=FALSE,scheme=scheme,tol=tol)
      }
      return(list(rgcca0=rgcca0,data=dataTest,rgccaList=resRgcca2))
    }
    else{stop("not enough neighbors with complete data (<5)")}
  }
  if(option=="em")
  {

     
       dataTest=resRgcca2=list()
      resImpute=imputeEM(A=A,tau=tau,C=C,scheme=scheme,ncomp=ncomp,superblock=superblock,naxis = naxis)
      dataTest0=resImpute$A
      rgcca0=rgccad(dataTest0,ncomp=rep(2,length(A)),scale=scale,sameBlockWeight=sameBlockWeight,tau=tau,verbose=FALSE,scheme=scheme,tol=tol)
      
      for(i in 1:ni)
      {
        print(i)
         dataTest[[i]]=addNoise(resImpute)
        resRgcca2[[i]]=rgccad(dataTest[[i]],ncomp=rep(2,length(dataTest[[i]])),scale=scale,sameBlockWeight=sameBlockWeight,tau=tau,verbose=FALSE,scheme=scheme,tol=tol)
      }
  }
    obj=list(rgcca0=rgcca0,data=dataTest,rgccaList=resRgcca2)
    class(obj) <- "mirgcca"
  return(obj)
}
