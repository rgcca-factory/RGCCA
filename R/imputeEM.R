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
imputeEM=function(A,C,tau,ni=50,tol=1e-8,graph=FALSE,ncomp=NULL,naxis=1,scale=TRUE,sameBlockWeight=TRUE,bias=TRUE,superblock=FALSE)
{
    #listWithoutNA
  nvar=sapply(A,NCOL)
  indNA=which(is.na(do.call(cbind,A)),arr.ind=TRUE)
  
   # to get the superblock concatenedBlocks", with missing values imputed by the colmeans
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
  concatenedBlocks = Reduce(cbind,Alist)
  # center and normalized ???
  D=matrix(1,dim(concatenedBlocks)[1],dim(concatenedBlocks)[2])
  scaledConcatenedBlocks = scale2(concatenedBlocks, scale=scale,bias = TRUE)

    group=unlist(lapply(A,"NCOL"))
    debutBlock=c(1,1+cumsum(group)[1:(length(group)-1)])
    finBlock=cumsum(group)
    if(sameBlockWeight)
    {
      for(u in 1:length(finBlock))
      {
         var_group=sum(apply( scaledConcatenedBlocks[,debutBlock[u]:finBlock[u]],2,"cov2"))
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
    fit.rgcca = rgcca(A=Alist, tau = tau,C=C,
                ncomp = ncomp,
                        scheme = "factorial",
                           scale = scale, init = "svd",
                          verbose = FALSE, tol = tol,sameBlockWeight=sameBlockWeight)
   # si on veut le critere comme somme des deux composantes
   
   # Criteria for one component only
   nIterCrit=length(fit.rgcca$crit[[1]])
   if(ncomp[1]==1)
   {
     print("fit")
     critRGCCA=c(critRGCCA,fit.rgcca$crit[[1]][1])
     
   }
   if(ncomp[1]==2)
   {
     critRGCCA=c(critRGCCA,fit.rgcca$crit[[1]][1]+fit.rgcca$crit[[2]][1])
     
   }
    
   
   # Reconstruction to do
   # if(superblock)
   # {
   #   w = fit.rgcca$astar[[J+1]]
   #   if(is.matrix(w)&&dim(y)[2]>1){
   #     y=y[,1:naxis]
   #     gamma=NULL
   #     for(k in 1:naxis)
   #     { # regression of the relevant block on the y
   #       gamma=cbind(a,apply(D*scaledConcatenedBlocks, 2, function(x) lm(x~w[,k])$coefficients[2]))
   #     }
   #   } # if only 1 component is required
   #   else{a=apply(D*scaledConcatenedBlocks, 2, function(x) lm(x~w)$coefficients[2])}
   #   
   # }
   if(!superblock)
   {
     Xhat=list()
     for(j in 1:J)
     {
       w = fit.rgcca$astar[[j]]
       if(is.matrix(w)&&dim(w)[2]>1){
         if(naxis>1)
         {
           w=w[,1:naxis]
           gamma=cbind(gamma,apply(D[,debutBlock[j]:finBlock[j]]*scaledConcatenedBlocks[,debutBlock[j]:finBlock[j]], 1, function(x) lm(x~w[,1:naxis])$coefficients[2]))
         }
         else
         {
           w=w[,1:naxis]
           gamma=cbind(gamma,apply(D[,debutBlock[j]:finBlock[j]]*scaledConcatenedBlocks[,debutBlock[j]:finBlock[j]], 1, function(x) lm(x~w)$coefficients[2]))
         }
              } # if only 1 component is required
       else
       {
         gamma=apply(scaledConcatenedBlocks[,debutBlock[j]:finBlock[j]], 1, function(x) lm(x~w)$coefficients[2])
        # res=apply(D[,debutBlock[j]:finBlock[j]]*scale2(fit.rgcca$A[[j]], bias = TRUE,scale=scale), 1, function(x) lm(x~w)$residuals)
        # sigma=var(res)  
         
                  # mean and standard deviation of each variables are saved
     
       }
       
       moy = matrix(attr(scaledConcatenedBlocks, "scaled:center")[debutBlock[j]:finBlock[j]],
                    nrow = NROW(scaledConcatenedBlocks), ncol=length(debutBlock[j]:finBlock[j]), byrow = TRUE)
       stdev = matrix(attr(scaledConcatenedBlocks, "scaled:scale")[debutBlock[j]:finBlock[j]],
                      nrow = NROW(scaledConcatenedBlocks),ncol=length(debutBlock[j]:finBlock[j]), byrow = TRUE)
       Xhat[[j]] =gamma%*%t(w)*stdev*(1/D[,debutBlock[j]:finBlock[j]]) + moy
     }
     
     concatenedXhat= Reduce(cbind,Xhat)
   }                          
     # if at least 2 components are required
      
   concatenedBlocks[indNA] = concatenedXhat[indNA]
   # Faut il reinitialiser les scaledConcatenedBlocks dans la boucle ? 
   scaledConcatenedBlocks = scale2(concatenedBlocks, scale=scale,bias = TRUE)
   
   Alist=list()
   for(j in 1:length(nvar))
   {
     if(j==1){sel=1:nvar[1]}else{debut=sum(nvar[1:(j-1)])+1;fin=debut+(nvar[j]-1);sel=debut:fin}
     Alist[[j]]=as.matrix(scaledConcatenedBlocks[,sel])
     colnames( Alist[[j]])=colnames(A[[j]])
   }
   names(Alist)=names(A)
   
      
    diff[[i]] <-concatenedBlocks - concatenedXhat
    diff[[i]][indNA] <- 0
    objective[[i]]<-sqrt(sum(diff[[i]]^2)/(dim(concatenedBlocks)[1]*dim(concatenedBlocks)[2]-sum(is.na(concatenedBlocks))))
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
    print(objective[[i]])
    print(criterion[[i]])
    i<- i + 1
    print(i)
  }	
  
   if(graph)
   { 
     vec=rep(NA,i-1)
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
