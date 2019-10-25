#' The function rgccak() is called by rgcca() and does not have to be used by the user. 
#' The function rgccak() computes the RGCCA block components, outer weight vectors, etc., 
#' for each block and each dimension. Depending on the dimensionality of each block \eqn{X_j , j = 1, ..., J}, 
#' the primal (when \eqn{n > p_j}) or the dual (when \eqn{n < p_j}) algorithm is used (see Tenenhaus et al. 2015) 
#' @param A  A list that contains the \eqn{J} blocks of variables. Either the blocks (\eqn{X_1, X_2, ..., X_J}) or the residual matrices (\eqn{X_{h1}, X_{h2}, ..., X_{hJ}}).
#' @param C  A design matrix that describes the relationships between blocks. (Default: complete design).
#' @param tau A \eqn{1 * J} vector that contains the values of the shrinkage parameters \eqn{\tau_j}, \eqn{ j=1, ..., J}. (Default: \eqn{\tau_j = 1}, \eqn{ j=1, ..., J}).
#' If tau = "optimal" the shrinkage intensity paramaters are estimated using the Schafer and Strimmer (2005) 
#' analytical formula. 
#' @param scheme The value is "horst", "factorial", "centroid" or any diffentiable convex scheme function g designed by the user (default: "centroid").
#' @param verbose  Will report progress while computing if verbose = TRUE (default: TRUE).
#' @param init The mode of initialization to use in the RGCCA algorithm. The alternatives are either by Singular Value Decompostion or random (default : "svd").
#' @param bias A logical value for either a biaised or unbiaised estimator of the var/cov.
#' @param tol Stopping value for convergence.
#' @param na.rm If TRUE, NIPALS algorithm taking missing values into account is run. RGCCA is run only on available data.
#' @return \item{Y}{A \eqn{n * J} matrix of RGCCA outer components}
#' @return \item{Z}{A \eqn{n * J} matrix of RGCCA inner components}
#' @return \item{a}{A list of outer weight vectors}
#' @return \item{crit}{The values of the objective function to be optimized in each iteration of the iterative procedure.}
# #' @return \item{converg}{Speed of convergence of the algorithm to reach the tolerance.}
#' @return \item{AVE}{Indicators of model quality based on the Average Variance Explained (AVE): 
#' AVE(for one block), AVE(outer model), AVE(inner model).}
#' @return \item{C}{A design matrix that describes the relationships between blocks (user specified).}
#' @return \item{tau}{\eqn{1 * J} vector containing the value for the tau penalties applied to each of the \eqn{J} blocks of data (user specified)}
#' @return \item{scheme}{The scheme chosen by the user (user specified).}
#' @references Tenenhaus M., Tenenhaus A. and Groenen PJF (2017), Regularized generalized canonical correlation analysis: A framework for sequential multiblock component methods, Psychometrika, in press
#' @references Tenenhaus A., Philippe C., & Frouin V. (2015). Kernel Generalized Canonical Correlation Analysis. Computational Statistics and Data Analysis, 90, 114-131.
#' @references Tenenhaus A. and Tenenhaus M., (2011), Regularized Generalized Canonical Correlation Analysis, Psychometrika, Vol. 76, Nr 2, pp 257-284.
#' @references Schafer J. and Strimmer K., (2005), A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics. Statist. Appl. Genet. Mol. Biol. 4:32.
#' @title Internal function for computing the RGCCA parameters (RGCCA block components, outer weight vectors, etc.).
#' @export rgccak
#' @importFrom MASS ginv
#' @importFrom stats cor rnorm
#' @importFrom graphics plot
rgccak=function (A, C, tau = "optimal", scheme = "centroid",verbose = FALSE, init = "svd", bias = TRUE, tol = 1e-08,na.rm=TRUE,estimateNA="no",sameBlockWeight=TRUE) 
{
# A liste de matrices "blocs" dans un ordre precis (cf matrice connexion)
# C matrice de connexion, 
# tau = 0 ou 1
# scheme : fonction g
# scale: transformations appliquees aux blocs
# verbose : affichage
# init : initialisation
# biais : covariance estimee avec ou sans biais
# tol: critere d'arret de l'algorithme
     if(mode(scheme) != "function") 
    {
    if(!scheme %in% c("horst","factorial","centroid")){stop("Please choose scheme as 'horst','factorial','centroid' or as a convex function")}
    if(scheme=="horst"){ g <- function(x) x}
    if(scheme=="factorial"){ g <- function(x)  x^2}  
    if(scheme=="centroid"){g <- function(x) abs(x)}
      
    }else {g<-scheme}
    A0 <-A
    A <- lapply(A, as.matrix)
    if(estimateNA=="lebrusquet")
    {
        A =imputeColmeans(A) 
       # liste de blocs
    }
    # initialisation du A à la moyenne

   # A=lapply(A,function(M){M[is.na(M)]<-0;return(M)})
    J <- length(A) # nombre de blocs
    n <- NROW(A[[1]]) # nombre d'individus
    pjs <- sapply(A, NCOL) # nombre de variables par bloc
    Y <- matrix(0, n, J)
    if (!is.numeric(tau)) # cas ou on estime le tau de maniere intelligente (a creuser)
        tau = sapply(A, tau.estimate) # d'apres Schafer and Strimmer
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
        stop("init should be either random or by SVD.")
    }
   
    N = ifelse(bias, n, n - 1)
	# premiers reglages avant la boucle : initialisation du premier Y (correspondant a la fonction a maximiser)
    for (j in which.primal) 
    {
     	 ifelse(tau[j] == 1,
     	  {
            a[[j]] <- drop(1/sqrt(t(a[[j]]) %*% a[[j]])) * a[[j]] # calcul de la premiere composante (les aj sont les wj) : on les norme dans ce cas : c'eest la condition |w|=1
            # Remplir les valeurs de A manquantes avec la strategie
            Y[, j] <- pm(A[[j]] , a[[j]],na.rm=na.rm) # projection du bloc sur la premiere composante
           
           
        }, 
        {
           # M[[j]] <- ginv(tau[j] * diag(pjs[j]) + ((1 - tau[j])/(N)) * (pm(t(A[[j]]) , A[[j]],na.rm=na.rm))) #calcul de la fonction a minimiser ?
         
            
            #-taking NA into account in the N
            nmat=ifelse(bias,t(!is.na(A[[j]]))%*%(!is.na(A[[j]])),t(!is.na(A[[j]]))%*%(!is.na(A[[j]]))-1)
            nmat[nmat==0]=NA
            M[[j]] <- ginv(tau[j] * diag(pjs[j]) + ((1 - tau[j])) *nmat^(-1)* (pm(t(A[[j]]) , A[[j]],na.rm=na.rm))) #calcul de la fonction a minimiser ?
            #-----------------------
           
            a[[j]] <- drop(1/sqrt(t(a[[j]])%*% M[[j]]%*%a[[j]]) )* ( M[[j]] %*% a[[j]]) # calcul premiere composante (a creuser)
            Y[, j] <-pm( A[[j]] ,a[[j]],na.rm=na.rm) # projection du bloc sur la premiere composante
        })
    }
    for (j in which.dual)
    {
        ifelse(tau[j] == 1, {
            alpha[[j]] = drop(1/sqrt(t(alpha[[j]]) %*% K[[j]] %*%  alpha[[j]])) * alpha[[j]]
            a[[j]] =pm( t(A[[j]]), alpha[[j]],na.rm=na.rm)
            Y[, j] =pm( A[[j]], a[[j]],na.rm=na.rm)
        }, {
          
           # M[[j]] = tau[j] * diag(n) + (1 - tau[j])/(N) * K[[j]]  # contraire de la matrice de covariace
           
           #----taking NA into account in the N
            nmat=ifelse(bias,t(!is.na(A[[j]]))%*%(!is.na(A[[j]])),t(!is.na(A[[j]]))%*%(!is.na(A[[j]]))-1)
            nmat[nmat==0]=NA
            M[[j]] <- tau[j] * diag(n) + ((1 - tau[j])) *nmat^(-1)* K[[j]] #calcul de la fonction a minimiser ?
            #-----------------------
            
             Minv[[j]] = ginv(M[[j]])
            alpha[[j]] = drop(1/sqrt(t(alpha[[j]])%*% M[[j]]%*% K[[j]]%*% alpha[[j]])) * alpha[[j]]
            a[[j]] =pm( t(A[[j]]), alpha[[j]],na.rm=na.rm)
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
    #  print(paste("iter",iter))
       for (j in which.primal)
      { # on parcourt les blocs pour estimer wj = a[[j]] : c'est le rouage de la pres
         
          dgx = dg(cov2(Y[, j], Y, bias = bias))# covariance entre les differents blocs: dgx indique + - 1
          if(tau[j] == 1)
          { # si tau = 1
             Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * matrix(rep(dgx, n), n, J, byrow = TRUE) * Y,na.rm=na.rm)
		     a[[j]] = drop(1/sqrt(pm(pm(t(Z[, j]) ,A[[j]],na.rm=na.rm) ,  pm( t(A[[j]]) ,Z[, j],na.rm=na.rm),na.rm=na.rm))) *pm (t(A[[j]]), Z[,  j],na.rm=na.rm)  
		      Y[, j] =pm( A[[j]], a[[j]],na.rm=na.rm) #Nouvelle estimation de j
		
#------------------ si on estime les donnees manquantes dans le cas ou tau=1
			      if(estimateNA %in% c("first","iterative","superblock","new","lebrusquet"))
			      { 
			         

			   #    print(paste("block",j))
			        if(estimateNA=="superblock"){missing=is.na(A0[[j]][,k])}
			        if(estimateNA=="superblock")
			        {
			            if(j==length(A))
			            {
			                for(k in 1:dim(A[[j]])[2])
			                {
			                    A[[j]][missing,k]=scale(a[[j]][k]*Z[missing,j])
			                }
			            }
			        }
			        for(k in 1:NCOL(A[[j]]))
			        {
			            if(estimateNA=="lebrusquet")
			            {
			                missing=is.na(A0[[j]][,k])
			                if(sum(missing)!=0)
			                { 
			                     sum(C * g(cov2(Y, bias = bias)),na.rm=na.rm)
			                  
			                    newx_k=leb(x_k=A[[j]][,k],missing,z=Z[,j],sameBlockWeight=sameBlockWeight,weight=sqrt(pjs[j]),argmax=ifelse(a[[j]][k]>0,TRUE,FALSE))
			                    A[[j]][,k]=newx_k
			                    
			                    # mean(xres)
			                    # t(xres)%*%xres
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
			        #       if(j==1 && k==2) 
			        #       {
			        #      #  print("avant scaling observed")
			        #       # print(A[[j]][!missing,k][1])
			        #        print("avt scaling missing")
			        #        print(Ainter[missing])
			        #       }
			        #       A[[j]][,k]=scale(Ainter)
			        #       if(j==1&& k==2)
			        #       {
			        #        # print("apres scaling observed")
			        #        # print(A[[j]][!missing,k][1])
			        #         print("apres scaling missing")
			        #         print(A[[j]][missing,k][1])
			        #       }
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
			        #    # print("apres boucle")
			        #    # print(A[[j]][!missing,k][1])
			        #   }
			        #  
			        # #A[[j]][,k]=scale(A[[j]][,k])
			        # #	 print( A[[j]][,k])
			        }
			        
			      }
			       
           }else
			      { # si tau different de 1
              Z[, j] = rowSums(matrix(rep(C[j, ], n), n,  J, byrow = TRUE) * matrix(rep(dgx, n), n,  J, byrow = TRUE) * Y,na.rm=na.rm)
             a[[j]] = drop(1/sqrt(pm(pm(t(Z[, j]) ,A[[j]],na.rm=na.rm) , pm( pm( M[[j]] , t(A[[j]]),na.rm=na.rm) , Z[, j],na.rm=na.rm),na.rm=na.rm))) * pm(M[[j]],pm( t(A[[j]]) ,Z[, j]))
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
              Y[, j] =pm( A[[j]], a[[j]],na.rm=na.rm)
           }, 
           {
            Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * matrix(rep(dgx, n), n, J, byrow = TRUE) * Y,na.rm=na.rm)
          #  alpha[[j]] = drop(1/sqrt(pm(pm(pm(t(Z[, j]), K[[j]] ),  Minv[[j]]) , Z[, j]))) * pm(Minv[[j]] , Z[,  j])
			alpha[[j]] = drop(1/sqrt(t(Z[, j])%*% K[[j]] %*% Minv[[j]]%*% Z[, j])) * (Minv[[j]] %*% Z[,  j])
                   
		   a[[j]] =pm( t(A[[j]]) , alpha[[j]],na.rm=na.rm)
            Y[, j] =pm( A[[j]], a[[j]],na.rm=na.rm)
          })
      }
    # print(lapply(A,head))

      crit[iter] <- sum(C * g(cov2(Y, bias = bias)),na.rm=na.rm)
      if (verbose & (iter%%1) == 0) 
      cat(" Iter: ", formatC(iter, width = 3, format = "d"), " Fit:", formatC(crit[iter], digits = 8, width = 10, format = "f"), " Dif: ", formatC(crit[iter] - crit_old, digits = 8, width = 10, format = "f"),   "\n")
      stopping_criteria = c(drop(crossprod(Reduce("c", mapply("-", a, a_old)))), crit[iter] - crit_old)
     
      # if (any(stopping_criteria < tol) | (iter > 1000)) # critere d'arret de la boucle
           if (any(abs(stopping_criteria) < tol) | (iter > 1000)) # critere d'arret de la boucle
               
          break
      crit_old = crit[iter]
      a_old <- a
      iter <- iter + 1
    }
    if (iter > 1000) 
        warning("The RGCCA algorithm did not converge after 1000 iterations.")
    if (iter < 1000 & verbose) 
        cat("The RGCCA algorithm converged to a stationary point after",  iter - 1, "iterations \n")
    if (verbose) 
        plot(crit[1:iter], xlab = "iteration", ylab = "criteria")
    AVEinner <- sum(C * cor(Y)^2/2)/(sum(C)/2)
	#	AVEinner=diag(cov(res$Y[[1]]))/sum(diag(cov(A[[1]] )))
    if(estimateNA!="no")
    {
        result <- list(Y = Y, a = a, crit = crit, AVE_inner = AVEinner,  C = C, tau = tau, scheme = scheme,A=A)
    }
    else{result <- list(Y = Y, a = a, crit = crit, AVE_inner = AVEinner,  C = C, tau = tau, scheme = scheme)}
    
       return(result)
}