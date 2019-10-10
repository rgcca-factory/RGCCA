#'Impute with superblock method
#'
#'This method is used for the implementation of EM algorithm for missing data
#'
#' @param A A list of J blocks
#' @param A A matrix with J dimensions (to be changed for superblock ?)
#' @param tau A vector of tau values with the same length as A
#' @param ni An integer for the maximal number of iterations before convergence
#' @param tol The stopping value for convergence.
#' @param graph if graph = TRUE, graphics are plotted/saved
#' @param ncomp vector containing the number of components per block in RGCCA
#' @param naxis number of component to select for the estimation of missing data
#' @param  scale  If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
#' @param sameBlockWeight A logical value indicating if the different blocks should have the same weight in the analysis (default, sameBlockWeight=TRUE)
#' @param scheme scheme chosene for RGCCA (is not useful when superblock=TRUE)
#'  @param bias A logical value indicating if variance should be biased or not
#' @param superblock Boolean, if TRUE, the missing values are estimated with the superblock
#' @return \item{A} A list of blocks imputed
#' @return \item{stab} Convergence criterion : abs(1-obj_k/obj_{k-1})
#' @return \item{obj} Vector containing the mean square error between the predict values and the original non missing values at each iteration
#'@return \item{crit} RGCCA criterion
#'@return \item{moy} Estimated mean obtained for each variable  (required for the addNoise function)
#'@return \item{stdev} Estimated standard deviations of each variable obtained  (required for the addNoise function)
#'@return \item{sigma}Estimated standard deviations for the noise obtained  (required for the addNoise function)
#'@return \item{indNA} Position of missing values

#' @title imputeSB: impute with superblock method
#' @examples
#'  data();...
imputeEM = function(A,
                    C,
                    tau,
                    ni = 50,
                    tol = 1e-8,
                    graph = FALSE,
                    ncomp = NULL,
                    naxis = 1,
                    scale = TRUE,
                    sameBlockWeight = TRUE,
                    scheme = "centroid",
                    bias = TRUE,
                    superblock = FALSE)
{
  #listWithoutNA
  nvar = sapply(A, NCOL)
  nsuj = dim(A[[1]])[1]
  J = length(A)
  # Getting missing values pattern in the superblock
  indNA = which(is.na(do.call(cbind, A)), arr.ind = TRUE)
  indNA2 = lapply(A, function(x) {
    return(which(is.na(x), arr.ind = TRUE))
  })
  # Imputing Alist by the mean
  Alist = lapply(A, function(X) {
    if (is.list(X)) {
      X = as.matrix(X)
    }
    if (is.matrix(X))
    {
      m = apply(X, 2, function(x) {
        return(mean(x, na.rm = TRUE))
      })
      indNA = which(is.na(X), arr.ind = TRUE)
      X[indNA] = m[indNA[, 2]]
      return(X)
    }
    if (is.vector(X))
    {
      m = mean(X, na.rm = TRUE)
      X[is.na(X)] = m
      return(X)
    }
  })
  # Calculating starts and ends of each block in the superbblock
  group = unlist(lapply(A, "NCOL"))
  debutBlock = c(1, 1 + cumsum(group)[1:(length(group) - 1)])
  finBlock = cumsum(group)
  
  # to get the superblock concatenedBlocks", with missing values imputed by the colmeans
  concatenedBlocks = Reduce(cbind, Alist)
  
  
  # Normalize the blocks
  scaledConcatenedBlocks = scale2(concatenedBlocks, scale = scale, bias = TRUE)
  
  # calculating the D matrix with the weight by number of variables in each block
  # D=matrix(1,dim(concatened
 # Blocks)[1],dim(concatenedBlocks)[2])
# if(sameBlockWeight)
# {
#   for(u in 1:length(finBlock))
#   {
#     var_group=sum(apply( scaledConcatenedBlocks[,debutBlock[u]:finBlock[u]],2,"cov2"))
#     D[,debutBlock[u]:finBlock[u]]=1/sqrt(var_group)
#   }
# }
#
# initialization
i=1

diff=objective=old=criterion=list()
criterionFinal=objectiveFinal=list()
continue=TRUE
if(is.null(ncomp)){ncomp2=c(rep(1, J), naxis)}else{ncomp2=c(ncomp,naxis)}
critRGCCA=c()
old[[1]]=-Inf

if(superblock)
{
  C2=matrix(1,J+1,J+1);C2[1:J,1:J]=0; C2[J+1,J+1]=0
  tau2=c(tau,0) # mode A for the blocks and mode B for superblock (design for CPCA2 )
  
}
while (continue)
{
  diff[[i]]=objective[[i]]=old[[i+1]]=criterion[[i]]=list()
  
  # building of a list with superblock
  if(superblock)
  {
    # creation of the superblock
    ASB = c(Alist, list(concatenedBlocks))
    names(ASB)=c(names(Alist),"Superblock")
    fit.rgcca = rgcca(A=ASB, tau = tau2,C=C2,
                      ncomp = ncomp2,
                      scheme = "factorial",
                      scale = scale, init = "svd",
                      verbose = FALSE, tol = tol,sameBlockWeight=sameBlockWeight)
    
  }
  if(!superblock)
  {
    fit.rgcca = rgcca(A=Alist, tau = tau,C=C,
                      ncomp = ncomp,
                      scheme = scheme,
                      scale = scale, init = "svd",
                      verbose = FALSE, tol = tol,sameBlockWeight=sameBlockWeight)
  }
  
  
  # si on veut le critere comme somme des deux composantes
  
  # Criteria for one component only
  if(naxis==1)
  {
    critRGCCA=c(critRGCCA,fit.rgcca$crit[[1]][1])
  }
  # Criteria for one component only
  #  if(naxis==2)
  #  {
  #    critRGCCA=c(critRGCCA,fit.rgcca$crit[[1]][1]+fit.rgcca$crit[[2]][1])
  #
  #  }
  # if(naxis==3)
  # {
  #   critRGCCA=c(critRGCCA,fit.rgcca$crit[[1]][1]+fit.rgcca$crit[[2]][1]+fit.rgcca$crit[[3]][1])
  # }
  
  if(superblock)
  {
    # Getting back the global weights
    w = fit.rgcca$astar[[J+1]]
    moy = matrix(attr(scaledConcatenedBlocks, "scaled:center"), nrow = NROW(scaledConcatenedBlocks), ncol=sum(nvar), byrow = TRUE)
    stdev = matrix(attr(scaledConcatenedBlocks, "scaled:scale"),
                   nrow = NROW(scaledConcatenedBlocks),ncol=sum(nvar), byrow = TRUE)
    
    # TO DO for more than two components required
    if(naxis>1)
    {
      gamma=NULL
      sigma=NULL
      for(k in 1:naxis)
      { # regression of the relevant block on the y
        gamma=cbind(gamma,apply(scaledConcatenedBlocks, 1, function(x) lm(x~0+w[,k])$coefficients[1]))
        residuals=apply(scaledConcatenedBlocks, 1, function(x) (lm(x~0+w)$residuals))
        sigma=cbind(sigma,sqrt(sum(residuals^2/(J*nsuj))))
        # if(noise){eps=matrix(rnorm(dim(moy)[1]*dim(moy)[2],m=0,sd=sigma),nrow=dim(moy)[1],ncol=dim(moy)[2])}else{eps=0}
      }
      Xhat =(gamma%*%t(w))*stdev + moy
      concatenedXhat= Xhat
    }
    # if only 1 component is required
    if(naxis==1)
    {
      
      w1=w[,1]
      gamma=apply(scaledConcatenedBlocks, 1, function(x) lm(x~0+w)$coefficients[1])
      #=t(as.matrix(w))%*%t(scaledConcatenedBlocks)/sum(w*w)
      # sigma est la somme des résidus
      residuals=apply(scaledConcatenedBlocks, 1, function(x) (lm(x~0+w)$residuals))
      sigma=sqrt(sum(residuals^2/(J*nsuj)))
      centeredXhat=gamma%*%t(w)
      #print(dim(centeredXhat))
      #if(noise){eps=matrix(rnorm(dim(moy)[1]*dim(moy)[2],m=0,sd=sigma),nrow=dim(moy)[1],ncol=dim(moy)[2])}else{eps=0}
      Xhat =(centeredXhat)*stdev + moy
      print("sd")
      print(stdev[1,1:5])
      print("moy")
      print(moy[1,1:5])
      print(Xhat[1:5,1:5])
      
      concatenedCenteredXhat=centeredXhat
      concatenedXhat= Xhat
    }
    #Russett[1:10,]
  }
  if(!superblock)
  {
    Xhat=list();centeredXhat=list()
    gamma=NULL
    stdev=moy=sigma=list()
    for(j in 1:J)
    {
      w = fit.rgcca$astar[[j]]
      w=w[,1:naxis]
      moy[[j]] = matrix(attr(scaledConcatenedBlocks, "scaled:center")[debutBlock[j]:finBlock[j]],
                        nrow = NROW(scaledConcatenedBlocks), ncol=length(debutBlock[j]:finBlock[j]), byrow = TRUE)
      stdev[[j]]  = matrix(attr(scaledConcatenedBlocks, "scaled:scale")[debutBlock[j]:finBlock[j]],
                           nrow = NROW(scaledConcatenedBlocks),ncol=length(debutBlock[j]:finBlock[j]), byrow = TRUE)
      
      if(naxis==1)
      {
        gamma=apply(scaledConcatenedBlocks[,debutBlock[j]:finBlock[j]], 1, function(x) lm(x~0+w)$coefficients[1])
        residuals=apply(scaledConcatenedBlocks[,debutBlock[j]:finBlock[j]], 1, function(x) (lm(x~0+w)$residuals))
        sigma[[j]]=sqrt(sum(residuals^2/(J*nsuj)))
        
        # mean and standard deviation of each variables are saved
        centeredXhat[[j]]=matrix(gamma,ncol=naxis)%*%matrix(w,nrow=naxis)
        #if(noise){eps=matrix(rnorm(dim(moy)[1]*dim(moy)[2],m=0,sd=sigma),nrow=dim(moy)[1],ncol=dim(moy)[2])}else{eps=0}
        Xhat[[j]] = (centeredXhat[[j]])*stdev[[j]]  + moy[[j]]
      }
      # if(naxis>1)
      # {
      #   gamma=NULL
      #   residuals=NULL
      #   for(k in 1:naxis)
      #   { # regression of the relevant block on the y
      #     gamma=cbind(gamma,apply(scaledConcatenedBlocks, 1, function(x) lm(x~w[,k])$coefficients[2]))
      #     residuals=cbind(residuals,apply(scaledConcatenedBlocks, 1, function(x) lm(x~w[,k])$residuals))
      #   }
      #   Xhat[[j]] =(gamma%*%t(w))*stdev[[j]]  + moy[[j]]
      # }
    }
    concatenedCenteredXhat=Reduce(cbind,centeredXhat)
    concatenedXhat= Reduce(cbind,Xhat)
  }
  
  # concatened blocks is imputed by concatenedXhat  and Alist, concatenedBlocks and  scaledConcatenedBlock are reaffected
  concatenedBlocks[indNA] = concatenedXhat[indNA]
  scaledConcatenedBlocks = scale2(concatenedBlocks, scale=scale,bias = TRUE)
  
  # Criteria of convergence are calculated
  #----------------------------------------
  diff[[i]] <-concatenedBlocks - concatenedXhat
  # diff[[i]]=scaledConcatenedBlocks-concatenedCenteredXhat
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
  
  # When the algorithm stops, some noise can be added to the model
  #----------------------------------------------------------------
  
  i<- i + 1
  print(i)
  # Reaffecting Alist
  #----------------------
  Alist=list()
  for(j in 1:length(nvar))
  {
    if(j==1){sel=1:nvar[1]}else{debut=sum(nvar[1:(j-1)])+1;fin=debut+(nvar[j]-1);sel=debut:fin}
    Alist[[j]]=as.matrix(concatenedBlocks[,sel])
    colnames( Alist[[j]])=colnames(A[[j]])
  }
  names(Alist)=names(A)
}

if(graph)
{
  vec=rep(NA,i-1)
  for(k in 1:(i-1)){vec[k]=criterion[[k]]  }
  # x11()
  png("objective.png")
  plot(vec[-1], xlab = "iteration", ylab = "objective",pch=16)
  dev.off()
  
  png("criterion.png")
  plot(critRGCCA, xlab = "iteration", ylab = "RGCCA criterion",pch=16)
  dev.off()
}
return(list(A=Alist,crit=unlist(critRGCCA),stab=unlist(criterion),obj=unlist(objective),stdev=stdev,sigma=sigma,moy=moy,indNA=indNA2))
}
