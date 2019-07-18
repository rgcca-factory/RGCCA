#'Analysis of the comparison of different NA methods on RGCCA
#'
#' @param refData A list of complete blocks
#' @param noIntersect if TRUE, no RGCCA on complete data is run
#' @param C,tau,scale,scheme,sameBlockWeight parameters of RGCCA
#' @
#' @return \item{A} A list containing: a: the correlation between axes, rv: the rv coefficient, bm the biomarkers
#' @return \item{crit} Convergence criterion : abs(1-obj_k/obj_{k-1})
#' @return \item{obj} Vector containing the mean square error between the predict values and the original non missing values at each iteration
#' @title comparison of two RGCCA results
#' @examples 
#'  data();...

analysis=function(refData,noIntersect=TRUE,C=NULL,tau=NULL,scale=TRUE,nAxe=2,scheme="centroid",sameBlockWeight=TRUE,wd=getwd(),blocknames=NULL,nbTestFiles=20)
{
  nBlock=length(refData)
  if(is.null(C)) {C=matrix(1,nBlock,nBlock);diag(C)=0;}
  if(is.null(tau)){tau=rep(1,nBlock)}
  listFinale=list()
  library(parallel)
  refRgcca=rgcca(refData[1:nBlock],C=C,ncomp=rep(nAxe,nBlock),verbose=FALSE,tau=tau,scale=scale,scheme=scheme,sameBlockWeight=sameBlockWeight,returnA=TRUE)
  
  if (Sys.info()["sysname"] == "Windows"){
    listFinale= lapply(1:nbTestFiles, FUN=function(j) {
      print(paste("dataset",j))
      setwd(paste0(wd,"/",j))
      testData=readDataset(bloc = blocknames)
      listRgcca=list()
      #	lAxes=lapply(testData,function(m){resPca=PCA(m,graph=FALSE);eigenValues=resPca$eig[,1];nbAxes=critereCoude(eigenValues,graph=TRUE);return(nbAxes)})	
      #	nbAxesRgcca=as.vector(unlist(lAxes))		
      # listRgcca[["MI-kNN2"]]=MIRGCCA(testData,k=2,niter=5,scale=scale,sameBlockWeight=TRUE,tau,output="weightedMean",scheme=scheme)$rgcca0
      listRgcca[["MI-kNNAll"]]=MIRGCCA(testData,k= "all",ni=5,scale=scale,sameBlockWeight=sameBlockWeight,tau=tau,output="mean",scheme=scheme,returnA=TRUE,tol=1e-8)$rgcca0
      testDataSB=imputeSB(testData,ncomp=rep(nAxe,nBlock),scale=scale,sameBlockWeight=sameBlockWeight,tau=tau,tol=1e-8,ni=10)
      listRgcca[["EM"]]=rgcca(testDataSB$A,C=C,ncomp=rep(nAxe,nBlock),verbose=FALSE,scale=scale,sameBlockWeight=sameBlockWeight,tau=tau,scheme=scheme,returnA=TRUE,tol=1e-8)
      listRgcca[["Nipals"]]=rgcca(testData,C=C,ncomp=rep(nAxe,nBlock),verbose=FALSE,na.rm=TRUE,scale=scale,sameBlockWeight=sameBlockWeight,tau=tau,scheme=scheme,returnA=TRUE,tol=1e-8)
      if(!noIntersect){
        completeData=intersection(testData)
        listRgcca[["Complete"]]=rgcca(completeData,C=C,ncomp=rep(nAxe,nBlock),verbose=FALSE,scale=scale,sameBlockWeight=sameBlockWeight,tau=tau,returnA=TRUE,tol=1e-8)
      }
      graphics.off()
      if(noIntersect){selectPatient=NULL}else{selectPatient=rownames(listRgcca[["Complete"]][["Y"]][[1]])}
      res=lapply(listRgcca, function(x)      {
        comp=comparison(rgcca1=refRgcca,rgcca2=x,selectPatient=selectPatient);return(comp)
      })
      setwd("./..")
      names(res)=names(listRgcca)
      return(res)
    })
  } else {
    listFinale= mlcapply(1:nbTestFiles, mc.cores=7,
                          FUN=function(j)
    {
      print("dataset")
      print(j)
      setwd(as.character(j))
      
      testData=readDataset(blocknames)
      listRgcca=list()
      #	lAxes=lapply(testData,function(m){resPca=PCA(m,graph=FALSE);eigenValues=resPca$eig[,1];nbAxes=critereCoude(eigenValues,graph=TRUE);return(nbAxes)})	
      #	nbAxesRgcca=as.vector(unlist(lAxes))		
      # listRgcca[["MI-kNN2"]]=MIRGCCA(testData,k=2,niter=5,scale=scale,sameBlockWeight=TRUE,tau,output="weightedMean",scheme=scheme)$rgcca0
      listRgcca[["MI-kNNAll"]]=MIRGCCA(testData,k= ,ni=5,scale=scale,sameBlockWeight=sameBlockWeight,tau=tau,output="weightedMean",scheme=scheme,returnA=TRUE,tol=1e-8)$rgcca0
      testDataSB=imputeSB(testData,ncomp=rep(nAxe,nBlock),scale=scale,sameBlockWeight=sameBlockWeight,tau=tau,tol=1e-8,ni=10)
      listRgcca[["EM"]]=rgcca(testDataSB$A,C=C,ncomp=rep(nAxe,nBlock),verbose=FALSE,scale=scale,sameBlockWeight=sameBlockWeight,tau=tau,scheme=scheme,returnA=TRUE,tol=1e-8)
      listRgcca[["Nipals"]]=rgcca(testData,C=C,ncomp=rep(nAxe,nBlock),verbose=FALSE,na.rm=TRUE,scale=scale,sameBlockWeight=sameBlockWeight,tau=tau,scheme=scheme,returnA=TRUE,tol=1e-8)
      if(!noIntersect)
      {
        completeData=intersection(testData)
        listRgcca[["Complete"]]=rgcca(completeData,C=C,ncomp=rep(nAxe,nBlock),verbose=FALSE,scale=scale,sameBlockWeight=sameBlockWeight,tau=tau,returnA=TRUE,tol=1e-8)
      }
      
      graphics.off()
      if(noIntersect){selectPatient=NULL}else{selectPatient=rownames(listRgcca[["Complete"]][["Y"]][[1]])}
      res=lapply(listRgcca, function(x)
      {
        comp=comparison(rgcca1=refRgcca,rgcca2=x,selectPatient=selectPatient);return(comp)
      }
      )
      setwd("./..")
      names(res)=names(listRgcca)
      return(res)
    })
  }
  return(listFinale)
}
