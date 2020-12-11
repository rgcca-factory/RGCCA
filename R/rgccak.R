#' The function rgccak() is called by rgccad() and does not have to be used by the user.
#' The function rgccak() computes the RGCCA block components, outer weight vectors, etc., 
#' for each block and each dimension. Depending on the dimensionality of each block \eqn{X_j , j = 1, ..., J}, 
#' the primal (when \eqn{n > p_j}) or the dual (when \eqn{n < p_j}) algorithm is used (see Tenenhaus et al. 2015) 
#' @inheritParams select_analysis
#' @inheritParams rgccaNa
#' @inheritParams rgccad
#' @param A  A list that contains the \eqn{J} blocks of variables from which block components are constructed.
#' It could be eiher the original matrices (\eqn{X_1, X_2, ..., X_J}) or the residual matrices (\eqn{X_{h1}, X_{h2}, ..., X_{hJ}}).
#' @param na.rm If TRUE, NIPALS algorithm taking missing values into account is run. RGCCA is run only on available data.
#' @param estimateNA to choose between "no","first","iterative","superblock","new","lebrusquet") TO BE DEVELOPPED
#' @param initImpute 'rand' or 'mean': initialization for optimization method
#' @return \item{Y}{A \eqn{n * J} matrix of RGCCA outer components}
#' @return \item{Z}{A \eqn{n * J} matrix of RGCCA inner components}
#' @return \item{a}{A list of \eqn{J} elements. Each element of \eqn{a} is a matrix that contains the outer weight vectors for each block.}
# #' @return \item{converg}{Speed of convergence of the algorithm to reach the tolerance.}
#' @return \item{AVE}{A list of numerical values giving the indicators of model quality based on the Average Variance Explained (AVE): AVE(for each block), AVE(outer model), AVE(inner model).}
#' @return \item{call}{Call of the function}
#' @return \item{crit}{A vector of integer that contains for each component the values of the analysis criteria across iterations.}
#' @return \item{tau}{Either a 1*J vector or a \eqn{\mathrm{max}(ncomp) \times J} matrix containing the values
#' of the regularization parameters .Tau varies from 0 (maximizing the correlation) to 1 (maximizing the covariance).
#' If tau = "optimal" the regularization paramaters are estimated for each block and each dimension using the Schafer and Strimmer (2005)
#' analytical formula . If tau is a \eqn{1\times J} vector, tau[j] is identical across the dimensions of block \eqn{\mathbf{X}_j}.
#' If tau is a matrix, tau[k, j] is associated with \eqn{\mathbf{X}_{jk}} (\eqn{k}th residual matrix for block \eqn{j}). It can be estimated by using \link{rgcca_permutation}.}
#' @references Tenenhaus M., Tenenhaus A. and Groenen PJF (2017), Regularized generalized canonical correlation analysis: A framework for sequential multiblock component methods, Psychometrika, in press
#' @references Tenenhaus A., Philippe C., & Frouin V. (2015). Kernel Generalized Canonical Correlation Analysis. Computational Statistics and Data Analysis, 90, 114-131.
#' @references Tenenhaus A. and Tenenhaus M., (2011), Regularized Generalized Canonical Correlation Analysis, Psychometrika, Vol. 76, Nr 2, pp 257-284.
#' @references Schafer J. and Strimmer K., (2005), A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics. Statist. Appl. Genet. Mol. Biol. 4:32.
#' @title Internal function for computing the RGCCA parameters (RGCCA block components, outer weight vectors, etc.).
#' @importFrom MASS ginv
#' @importFrom stats cor rnorm
#' @importFrom graphics plot
#' @importFrom Deriv Deriv
rgccak=function (A, C, tau = "optimal", scheme = "centroid",verbose = FALSE, init = "svd", bias = TRUE, tol = 1e-08,na.rm=TRUE,estimateNA="no",scale=TRUE,scale_block=TRUE,initImpute="rand")
{  
    call=list(A=A, C=C, scheme = scheme,verbose = verbose, init = init, bias = bias, tol = tol,na.rm=na.rm,estimateNA=estimateNA,scale=scale,scale_block=scale_block,initImpute=initImpute)
        
     if(mode(scheme) != "function") 
    {
    if(!scheme %in% c("horst","factorial","centroid")){stop_rgcca("Please choose scheme as 'horst','factorial','centroid' or as a convex function")}
    if(scheme=="horst"){ g <- function(x) x}
    if(scheme=="factorial"){ g <- function(x)  x^2}  
    if(scheme=="centroid"){g <- function(x) abs(x)}
   
     }else {g<-scheme}
    
    
    J <- length(A) # nombre de blocs
    n <- NROW(A[[1]]) # nombre d'individus
    pjs <- sapply(A, NCOL) # nombre de variables par bloc
    Y <- matrix(0, n, J)
    if (!is.numeric(tau)) # cas ou on estime le tau de maniere intelligente (a creuser)
        tau = sapply(A, tau.estimate) # d apres Schafer and Strimmer
    

    A0 <-A
    A <- lapply(A, as.matrix)
    # Initialisation of missing values
     if(estimateNA=="lebrusquet")
     {
       #  A =imputeColmeans(A) 
         
         
         
         A = lapply(A,function(x)
                {
                    res=apply(x,2,function(t)
                            {
                                m=mean(t,na.rm=TRUE)
                                s=sd(t,na.rm=TRUE)
                                t[is.na(t)]=rnorm(sum(is.na(t),mean=m,sd=s))
                                return(t)
                            });
                            return(res)
                 })
       #  A=lapply(A,scale2,bias=TRUE)
       #  if(scale_block){A=lapply(A,function(x){return(x/sqrt(NCOL(x)))})}
        # liste de blocs
     }
    if(estimateNA=="superblock")
    { # unscaled A ! ! ! 
        if(initImpute=="rand")
        {
            A = lapply(A,function(x)
            {
                res=apply(x,2,function(t)
                {
                    m=mean(t,na.rm=TRUE)
                    s=sd(t,na.rm=TRUE)
                    t[is.na(t)]=rnorm(sum(is.na(t)),mean=m,sd=s)
                    return(t)
                });
                return(res)
            })    
        }
        if(initImpute=="colMeans")
        {
            A =imputeColmeans(A) 
        }
        Binit=A
        if (scale == TRUE) 
        {
            A1 = lapply(A, function(x) scale2(x,scale=TRUE, bias = bias)) # le biais indique si on recherche la variance biaisee ou non
            if(scale_block)
            {
                A = lapply(A1, function(x) {y=x/sqrt(NCOL(x));return(y)} )
            }
            else
            {
                A=A1
            }
            # on divise chaque bloc par la racine du nombre de variables pour avoir chaque poids pour le meme bloc
        }
        if (scale == FALSE)
        { 
            
            A1 = lapply(A, function(x) scale2(x, scale=FALSE, bias = bias)) 
            if(scale_block)
            {
                A = lapply(A1, function(x) {covarMat=cov2(x,bias=bias);varianceBloc=sum(diag(covarMat)); return(x/sqrt(varianceBloc))})
            }
            else
            {
                A=A1
            }
            
        }
        means=lapply(A1,function(x){M=matrix(rep(attributes(x)$'scaled:center' ,dim(x)[1]),dim(x)[1],dim(x)[2],byrow=TRUE);return(M)})
        stdev=lapply(A1,function(x){M=matrix(rep(attributes(x)$'scaled:scale' ,dim(x)[1]),dim(x)[1],dim(x)[2],byrow=TRUE);return(M)})
        if(scale_block){stdev=lapply(stdev,function(x){return(x/sqrt(NCOL(x)))})}
        
    }
       a <- alpha <- M <- Minv <- K <- list() # initialisation variables internes
    which.primal <- which((n >= pjs) == 1) # on raisonne differement suivant la taille du bloc
    which.dual <- which((n < pjs) == 1)

    if (init == "svd") { #initialisation intelligente dans les differents cas (a creuser)
        for (j in which.primal) {
            a[[j]] <- initsvd(A[[j]]) # pas la
        }
        for (j in which.dual) { 
            alpha[[j]] <- initsvd(A[[j]])
            K[[j]] <-pm( A[[j]] , t(A[[j]]),na.rm=na.rm) #A*t(A) plutot que t(A)*A
        }
    }
    else if (init == "random") {
        for (j in which.primal) {
            a[[j]] <- rnorm(pjs[j]) # on initialise suivant la loi normale
        }
        for (j in which.dual) {
            alpha[[j]] <- rnorm(n)
            K[[j]] <- pm(A[[j]] , t(A[[j]]),na.rm=na.rm) #A*t(A) plutot que t(A)*A
        }
    }
    else {
        stop_rgcca("init should be either random or by SVD.")
    }
   
    N = ifelse(bias, n, n - 1)
	# premiers reglages avant la boucle : initialisation du premier Y (correspondant a la fonction a maximiser)
    for (j in which.primal) 
    {
        ifelse(tau[j] == 1,
        {
            a[[j]] <- drop(1/sqrt(t(a[[j]]) %*% a[[j]])) * a[[j]] # calcul de la premiere composante (les aj sont les wj) : on les norme dans ce cas : c'eest la condition |w|=1
            if(a[[j]][1]<0){a[[j]]=-a[[j]]}
            Y[, j] <- pm(A[[j]] , a[[j]],na.rm=na.rm) # projection du bloc sur la premiere composante
        }, 
        {
           # M[[j]] <- ginv(tau[j] * diag(pjs[j]) + ((1 - tau[j])/(N)) * (pm(t(A[[j]]) , A[[j]],na.rm=na.rm))) #calcul de la fonction a minimiser ?
            #-taking NA into account in the N
            nmat=matrix(N,pjs[j],pjs[j])
            nacol=unique(which(is.na(A[[j]]),arr.ind=T)[,"col"])
            if(length(nacol)!=0)
            {
                sujToRemove=is.na(t(A[[j]][,nacol]))%*%is.na(A[[j]][,nacol])
                nmat[nacol,nacol]=nmat[nacol,nacol]-sujToRemove
            }
             # 
            # if(bias)
            # {
            #     nmat= t(!is.na(A[[j]]))%*%(!is.na(A[[j]]))
            # }
            # else
            # {
            #     nmat= t(!is.na(A[[j]]))%*%(!is.na(A[[j]]))-1
            # }
            
            nmat[nmat==0]=NA
            M[[j]] <- ginv(tau[j] * diag(pjs[j]) + ((1 - tau[j])) *nmat^(-1)* (pm(t(A[[j]]) , A[[j]],na.rm=na.rm))) #calcul de la fonction a minimiser ?
            a[[j]] <- drop(1/sqrt(t(a[[j]])%*% M[[j]]%*%a[[j]]) )* ( M[[j]] %*% a[[j]]) # calcul premiere composante (a creuser)
            if(a[[j]][1]<0){a[[j]]=-a[[j]]}
            Y[, j] <-pm( A[[j]] ,a[[j]],na.rm=na.rm) # projection du bloc sur la premiere composante
        })
    }
    for (j in which.dual)
    {
        ifelse(tau[j] == 1, {
            alpha[[j]] = drop(1/sqrt(t(alpha[[j]]) %*% K[[j]] %*%  alpha[[j]])) * alpha[[j]]
            a[[j]] =pm( t(A[[j]]), alpha[[j]],na.rm=na.rm)
            if(a[[j]][1]<0){a[[j]]=-a[[j]]}
            Y[, j] =pm( A[[j]], a[[j]],na.rm=na.rm)
        }, {
          
           # M[[j]] = tau[j] * diag(n) + (1 - tau[j])/(N) * K[[j]]  # contraire de la matrice de covariace
           
           #----taking NA into account in the N
             
             # if(bias)
             # {
             #    nmat=matrix()
             #    nmat= t(!is.na(A[[j]]))%*%(!is.na(A[[j]]))
             # }
             # else
             # {
             #    nmat= t(!is.na(A[[j]]))%*%(!is.na(A[[j]]))-1
             # }
             #            
                        
         
             #nmat[nmat==0]=NA
             
             # nmat=matrix(N,n,n)
             # which()
          
            # M[[j]] <- tau[j] * diag(n) + ((1 - tau[j])) *nmat^(-1)* K[[j]] #calcul de la fonction a minimiser ?
            M[[j]] <- tau[j] * diag(n) + ((1 - tau[j])) *1/N* K[[j]] #calcul de la fonction a minimiser ?
            Minv[[j]] = ginv(M[[j]])
            alpha[[j]] = drop(1/sqrt(t(alpha[[j]])%*% M[[j]]%*% K[[j]]%*% alpha[[j]])) * alpha[[j]]
            a[[j]] =pm( t(A[[j]]), alpha[[j]],na.rm=na.rm) 
            if(a[[j]][1]<0){a[[j]]=-a[[j]]}
            Y[, j] = pm(A[[j]] ,a[[j]],na.rm=na.rm) 
        })
    }

			# ajout de na.rm=TRUE
    crit_old <- sum(C * g(cov2(Y, bias = bias)),na.rm=na.rm)# critere d'arret: h(cov(Y))
    iter = 1
    crit = numeric()
    Z = matrix(0, NROW(A[[1]]), J)
    a_old = a
    
    dg = Deriv::Deriv(g, env = parent.frame())# on derive la fonction g
   
    repeat 
    { # on rentre dans la boucle a proprement parler
      Yold <- Y #valeur de f

       for (j in which.primal)
      { # on parcourt les blocs pour estimer wj = a[[j]] : c'est le rouage de la pres

          dgx = dg(cov2(Y[, j], Y, bias = bias))# covariance entre les differents blocs: dgx indique + - 1
          if(tau[j] == 1)
          { # si tau = 1
             Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * matrix(rep(dgx, n), n, J, byrow = TRUE) * Y,na.rm=na.rm)
		     a[[j]] = drop(1/sqrt(pm(pm(t(Z[, j]) ,A[[j]],na.rm=na.rm) ,  pm( t(A[[j]]) ,Z[, j],na.rm=na.rm),na.rm=na.rm))) *pm (t(A[[j]]), Z[,  j],na.rm=na.rm)  
		     if(a[[j]][1]<0){a[[j]]=-a[[j]]}
		     # Y[, j] =pm( A[[j]], a[[j]],na.rm=na.rm) #Nouvelle estimation de j
		
#------------------ si on estime les donnees manquantes dans le cas ou tau=1
			      if(estimateNA %in% c("first","iterative","superblock","new","lebrusquet"))
			      { 
	    	        if(estimateNA=="superblock")
			        {
			             for(k in 1:dim(A[[j]])[2])
		                {
			                 missing=is.na(A0[[j]][,k])
                            if(sum(missing)>0)
                            {
                                S_k=sum(missing)*sum(A[[j]][!missing,k]^2,na.rm=TRUE)/sum(!missing)# Variance part allowed for missing data
                                x_miss=sqrt(S_k)*a[[j]][k]*missing*Z[,j]/sqrt(sum((a[[j]][k]*missing*Z[,j])^2,na.rm=TRUE))
                                Binit[[j]][,k]=(means[[j]][,k]*missing+stdev[[j]][,k]*missing*x_miss)+(rep(1,length(missing))-missing)*Binit[[j]][,k]
                                
                            }
                      }
	    	   
			            if (scale == TRUE) 
			            {
			                A1[[j]]= scale2(Binit[[j]],scale=TRUE, bias = bias)# le biais indique si on recherche la variance biaisee ou non
			                if(scale_block)
			                {
			                    A[[j]]= A1[[j]] /sqrt(NCOL(A[[j]]))
			                }
			                else
			                {
			                    A[[j]]=A1[[j]]
			                }
			                # on divise chaque bloc par la racine du nombre de variables pour avoir chaque poids pour le meme bloc
			            }
			            if (scale == FALSE)
			            { 
			                
			                A1[[j]] = scale2(Binit[[j]], scale=FALSE, bias = bias)
			                if(scale_block)
			                {   covarMat=cov2(A1[[j]],bias=bias);varianceBloc=sum(diag(covarMat))
			                    A[[j]] = A1[[j]]/sqrt(varianceBloc)
			                }
			                else
			                {
			                    A[[j]]=A1[[j]]
			                }
			                
			            }
	    	         
			            means[[j]]=matrix(rep(attributes(A1[[j]])$'scaled:center' ,dim(A1[[j]])[1]),dim(A1[[j]])[1],dim(A1[[j]])[2],byrow=TRUE)
			            stdev[[j]]=matrix(rep(attributes(A1[[j]])$'scaled:scale' ,dim(A1[[j]])[1]),dim(A1[[j]])[1],dim(A1[[j]])[2],byrow=TRUE)
			            if(scale_block){ stdev[[j]]=stdev[[j]]*sqrt(NCOL(A1[[j]]))}
			            
			        }
			       # Ainter=A
    			     if(estimateNA=="lebrusquet")
    			     {
        		        for(k in 1:NCOL(A[[j]]))
        		        {
			           
			                
			                 missing=is.na(A0[[j]][,k])
			                if(sum(missing)!=0)
			                { 
			                    
			                    #if(k==1)
			                    #{
			                        title=paste("bloc",j,",var",k,"iter",iter)
			                        png(filename=paste(title,".png",sep=""))
			                        newx_k=leb(x_k=A[[j]][,k],missing,z=Z[,j],scale_block=TRUE,weight=sqrt(pjs[j]),argmax=ifelse(a[[j]][k]>0,TRUE,FALSE),graph=FALSE,main=title,abscissa=A0[[j]][,k])
			                        dev.off()
			                    #}
    			                 #else
    			                 #{
    			                 #    newx_k=leb(x_k=A[[j]][,k],missing,z=Z[,j],scale_block=TRUE,weight=sqrt(pjs[j]),argmax=ifelse(a[[j]][k]>0,TRUE,FALSE),graph=FALSE,main=title)
    			                     
    			                 #}
			                    # on affecte le nouveau k
			                    A[[j]][,k]=newx_k
			                   # Ainter[[j]][,k]=newx_k
			        
			                }
			                
			            } 
			            
			        # 
			        #  # if(estimateNA=="first"){missing=is.na(A[[j]][,k])} # n'est valable qu'au premier
			        #   if(estimateNA=="first"){missing=is.na(A0[[j]][,k])} # n'est valable qu'au premier
			        #   #if(estimateNA=="iterative"){missing=is.na(A[[j]][,k])} 
			        #   if(estimateNA=="iterative"){  missing=is.na(A0[[j]][,k])}
			        #   if(sum(missing)>0)
			        #   {
			        #      #if(estimateNA=="first")
			        #     #{
			        #     #  A[[j]][missing,k]=a[[j]][k]*Z[missing,j];
			        #     #  A[[j]][,k]=scale(A[[j]][,k])
			        #     #}
			        #     if(estimateNA=="first")
			        #     {
			        #       #A[[j]][missing,k]=a[[j]][k]*Z[missing,j];
			        # 
			        #       Ainter= A0[[j]][,k]
			        #       Ainter[missing]=a[[j]][k]*Z[missing,j]
			        			        #       A[[j]][,k]=scale(Ainter)

			        #     }
			        #     if(estimateNA=="new")
			        #     {
			        #       ones=c(1,length(missing))
			        #       K=diag(!missing)
			        #       resol=solveLagrangien(K,w,x_obs,n)
			        #       alpha=resol[1]
			        #       mu=resol[2]
			        #       lambda=resol[3]
			        #       nu=resol[4]
			        #       x_miss=(a[[j]][k]*Z[missing,j]+nu%*%ones)/lambda
			        #     }
			        #    
			        #    if(estimateNA=="iterative")
			        #    {
			        #      A[[j]][missing,k]=scale(a[[j]][k]*Z[missing,j])
			        #    }
			        #    
			       #   }
			        #  
			        # #A[[j]][,k]=scale(A[[j]][,k])
			        # #	 print( A[[j]][,k])
			        }
			        
			      }
		     Y[, j] =pm( A[[j]], a[[j]],na.rm=na.rm)
			       
           }else
			      { # si tau different de 1
              Z[, j] = rowSums(matrix(rep(C[j, ], n), n,  J, byrow = TRUE) * matrix(rep(dgx, n), n,  J, byrow = TRUE) * Y,na.rm=na.rm)
             a[[j]] = drop(1/sqrt(pm(pm(t(Z[, j]) ,A[[j]],na.rm=na.rm) , pm( pm( M[[j]] , t(A[[j]]),na.rm=na.rm) , Z[, j],na.rm=na.rm),na.rm=na.rm))) * pm(M[[j]],pm( t(A[[j]]) ,Z[, j]))
             if(a[[j]][1]<0){a[[j]]=-a[[j]]}
             Y[, j] = pm(A[[j]] ,a[[j]],na.rm=na.rm)
          }
       }

      for (j in which.dual)
      {
          dgx = dg(cov2(Y[, j], Y, bias = bias))
          ifelse(tau[j] == 1, 
            {
              Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * matrix(rep(dgx, n), n, J, byrow = TRUE) * Y,na.rm=na.rm)
              alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Z[, j])) * Z[, j]
              a[[j]] =pm( t(A[[j]]) , alpha[[j]],na.rm=na.rm)
              if(a[[j]][1]<0){a[[j]]=-a[[j]]}
              Y[, j] =pm( A[[j]], a[[j]],na.rm=na.rm)
           }, 
           {
            Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * matrix(rep(dgx, n), n, J, byrow = TRUE) * Y,na.rm=na.rm)
          #  alpha[[j]] = drop(1/sqrt(pm(pm(pm(t(Z[, j]), K[[j]] ),  Minv[[j]]) , Z[, j]))) * pm(Minv[[j]] , Z[,  j])
			alpha[[j]] = drop(1/sqrt(t(Z[, j])%*% K[[j]] %*% Minv[[j]]%*% Z[, j])) * (Minv[[j]] %*% Z[,  j])
                   
		   a[[j]] =pm( t(A[[j]]) , alpha[[j]],na.rm=na.rm)
		   if(a[[j]][1]<0){a[[j]]=-a[[j]]}
            Y[, j] =pm( A[[j]], a[[j]],na.rm=na.rm)
          })
      }

      crit[iter] <- sum(C * g(cov2(Y, bias = bias)),na.rm=na.rm)
      if (verbose & (iter%%1) == 0) 
      {
          cat(" Iter: ", formatC(iter, width = 3, format = "d"), " Fit:", formatC(crit[iter], digits = 8, width = 10, format = "f"), " Dif: ", formatC(crit[iter] - crit_old, digits = 8, width = 10, format = "f"),   "\n")
      }
       stopping_criteria = c(drop(crossprod(Reduce("c", mapply("-", a, a_old)))), crit[iter] - crit_old)
     
       if (any(stopping_criteria < tol) | (iter > 1000)) # critere d'arret de la boucle
        {break}  
      crit_old = crit[iter]
      a_old <- a
      iter <- iter + 1
    }
    if (iter > 1000) 
        warning("The RGCCA algorithm did not converge after 1000 iterations.")
    if (iter < 1000 & verbose) 
        cat("The RGCCA algorithm converged to a stationary point after",  iter - 1, "iterations \n")
    if (verbose) 
    {
        plot(crit[1:iter], xlab = "iteration", ylab = "criteria")
    }
    AVEinner <- sum(C * cor(Y)^2/2)/(sum(C)/2)
    call$tau=tau
    
    if(estimateNA!="no")
    {
        result <- list(Y = Y, a = a, crit = crit, AVE_inner = AVEinner, A=A,call=call,tau=tau)
    }
    else{result <- list(Y = Y, a = a, crit = crit, AVE_inner = AVEinner,call=call,tau=tau)}
    return(result)
}