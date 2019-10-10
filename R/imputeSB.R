#'Impute with superblock method
#'
#'This method is used for the implementation of EM algorithm for missing data
#'
#' @param A A list of J blocks
#' @param tau A vector of tau values with the same length as A
#' @param tol The stopping value for convergence.
#' @param graph if graph = TRUE, 
#' @param ncomp vector containing the number of components per block in RGCCA
#' @param naxis number of component to select in the superblock for the estimation of missing data
#' @param  scale  If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
#' @param sameBlockWeight A logical value indicating if the different blocks should have the same weight in the analysis (default, sameBlockWeight=TRUE)
#' @param bias A logical value indicating if variance should be biased or not
#' @return \item{A} A list of blocks imputed
#' @return \item{crit} Convergence criterion : abs(1-obj_k/obj_{k-1})
#' @return \item{obj} Vector containing the mean square error between the predict values and the original non missing values at each iteration
#' @title imputeSB: impute with superblock method
#' @examples 
#'  data();...
imputeSB=function(A,tau,ni=50,tol=1e-8,graph=FALSE,ncomp=NULL,naxis=1,scale=TRUE,sameBlockWeight=TRUE,bias=TRUE)
{
    #listWithoutNA
  nvar=sapply(A,NCOL)
  indNA=which(is.na(do.call(cbind,A)),arr.ind=TRUE)
   # to get the superblock X1NA, with missing values imputed by the colmeans
  Alist=lapply(A,function(X){
      if(is.list(X)){X=as.matrix(X)}
      if(is.matrix(X))
      {
        m=apply(X, 2, function(x){return(mean(x, na.rm = TRUE))} )
        indNA=which(is.na(X),arr.ind = TRUE)
        X[indNA] = m[indNA[, 2]]
        return(X) 
      }
      if(is.vector(X))
      {
        m=mean(X,na.rm=TRUE)
        X[is.na(X)]=m
        return(X)
      }
      }
    )
  X1NA = Reduce(cbind,Alist)
  # center and normalized ???
  D=matrix(1,dim(X1NA)[1],dim(X1NA)[2])
  X2NA = scale2(X1NA, scale=scale,bias = TRUE)
   if(sameBlockWeight)
   {
      group=unlist(lapply(A,"NCOL"))
      debutBlock=c(1,1+cumsum(group)[1:(length(group)-1)])
      finBlock=cumsum(group)
      for(u in 1:length(finBlock))
      {
         var_group=sum(apply(X2NA[,debutBlock[u]:finBlock[u]],2,"cov2"))
         D[,debutBlock[u]:finBlock[u]]=1/sqrt(var_group)
      }
   } 
  # initialization
  i=1
  diff=objective=old=criterion=list()
  criterionFinal=objectiveFinal=list()
  continue=TRUE
  
  J=length(A)
  C2=matrix(1,J+1,J+1);C2[1:J,1:J]=0; C2[J+1,J+1]=0
  tau2=c(tau,0) # mode A for the blocks and mode B for superblock (design for CPCA2 )
  if(is.null(ncomp)){ncomp2=c(rep(1, J), naxis)}else{ncomp2=c(ncomp,naxis)}
  critRGCCA=c()
  old[[1]]=-Inf
  
  while (continue)
  { 
    diff[[i]]=objective[[i]]=old[[i+1]]=criterion[[i]]=list()
    
    # building of a list with superblock
 #  ASB = c(Alist, list(X1NA))
 #  names(ASB)=c(names(Alist),"Superblock")
 #  fit.rgcca = rgcca(A=ASB, tau = tau2,C=C2,
 #                            ncomp = ncomp2,
 #                            scheme = "factorial",
 #                            scale = scale, init = "svd",
 #                            verbose = FALSE, tol = tol,sameBlockWeight=sameBlockWeight)
    fit.rgcca = rgcca(A=Alist, tau = tau,C="superblock",
                ncomp = ncomp,
                        scheme = "factorial",
                           scale = scale, init = "svd",
                          verbose = FALSE, tol = tol,sameBlockWeight=sameBlockWeight)
   # si on veut le critere comme somme des deux composantes
  # critRGCCA=c(critRGCCA,fit.rgcca$crit[[1]][length(fit.rgcca$crit[[1]])]+fit.rgcca$crit[[2]][length(fit.rgcca$crit[[2]])])
  
   # Criteria for one component only
   nIterCrit=length(fit.rgcca$crit[[1]])
   critRGCCA=c(critRGCCA,fit.rgcca$crit[[1]][nIterCrit])
   
     y = fit.rgcca$Y[[J+1]]
     # if at least 2 components are required
    if(is.matrix(y)&&dim(y)[2]>1){
      y=y[,1:naxis]
      a=NULL
      for(k in 1:naxis)
      { # regression of the relevant block on the y
        a=cbind(a,apply(D*scale2(X1NA, bias = TRUE,scale=scale), 2, function(x) lm(x~y[,k])$coefficients[2]))
      }
    } # if only 1 component is required
    else{a=apply(D*scale2(X1NA, bias = TRUE,scale=scale), 2, function(x) lm(x~y)$coefficients[2])}
    
     X2NA = scale2(X1NA, scale=scale,bias = TRUE)
    # mean and standard deviation of each variables are saved
    moy = matrix(attr(X2NA, "scaled:center"),
                 nrow = NROW(X2NA), ncol=NCOL(X2NA), byrow = TRUE)
    stdev = matrix(attr(X2NA, "scaled:scale"),
                   nrow = NROW(X2NA),ncol=NCOL(X2NA), byrow = TRUE)
    
# Reconstruction ! ! !
    Xhat = (y%*%t(a))*stdev*(1/D) + moy
    
    # 
    X1NA[indNA] = Xhat[indNA]
    Alist=list()
    for(j in 1:length(nvar))
    {
      if(j==1){sel=1:nvar[1]}else{debut=sum(nvar[1:(j-1)])+1;fin=debut+(nvar[j]-1);sel=debut:fin}
      Alist[[j]]=as.matrix(X1NA[,sel])
      colnames( Alist[[j]])=colnames(A[[j]])
    }
    names(Alist)=names(A)
    
    diff[[i]] <-X1NA - Xhat
    diff[[i]][indNA] <- 0
    objective[[i]]<-sqrt(sum(diff[[i]]^2)/(dim(X1NA)[1]*dim(X1NA)[2]-sum(is.na(X1NA))))
    criterion[[i]] <- abs(1 - objective[[i]]/old[[i]])
    old[[i+1]] <- objective[[i]]

    if (!is.nan(criterion[[i]])) 
    {
      if (criterion[[i]] < tol && (i > 2) ) 
      {
        continue <- FALSE
      }
      if(objective[[i]]<tol)
      {
        continue=FALSE
        if ( verbose) 
        {
          cat("The algorithm converged to a stationary point after", i - 1, "iterations \n")
        }
     }
    }
    if (i > ni)
    {
      continue <- FALSE
      warning(paste("The RGCCA imputation algorithm did not converge after ", ni, " iterations"))
    }
    i<- i + 1
  }	
  
   if(graph)
   {   vec=rep(NA,i-1)
   for(k in 1:(i-1)){vec[k]=criterion[[k]]}
  # x11()
   png("objective.png")
   plot(vec[-1], xlab = "iteration", ylab = "objective",pch=16)
   dev.off()
  
   png("criterion.png")
   plot(critRGCCA, xlab = "iteration", ylab = "criterion",pch=16)
   dev.off()
  }

  return(list(A=Alist,crit=unlist(critRGCCA),stab=unlist(criterion)))
}
