#' impute with the rgccac function
#' @param A  A list that contains the \eqn{J} blocks of variables \eqn{X_1, X_2, ..., X_J}.
#' @param scale  If scale = TRUE, each block is standardized to zero means and unit variances and then divided by the square root of its number of variables (default: TRUE).
#' @param verbose  If verbose = TRUE, the progress will be report while computing (default: TRUE).
#' @param init The mode of initialization to use in RGCCA algorithm. The alternatives are either by Singular Value Decompostion ("svd") or random ("random") (Default: "svd").
#' @param bias A logical value for biaised or unbiaised estimator of the var/cov (default: bias = TRUE).
#' @param tol The stopping value for convergence.
#' @param estimationType type of estimation for missing values
#' @param initNA type of initialzation for missing values ("nipals", "mean"or "rand")
#' @param ni number maximal of iterations
#' @export imputec
#' @examples 
#'  data(Russett)
#' X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
#' X_ind = as.matrix(Russett[,c("gnpr","labo")]);
#' X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
#' X_polit[1,1]=NA
#' A = list(X_agric, X_ind, X_polit);
#'  result.rgccac1 = imputec(A, init="svd", scale = FALSE)
#'  

imputec=function(A, scale = TRUE , init="svd", bias = TRUE, tol = 1e-5, verbose=TRUE,estimationType="singleBlock",initNA="mean",ni=10)
{
    # checking parameters
    match.arg(initNA,c("mean","rand","nipals"))
    match.arg(estimationType,c("singleBlock","allInclusive","covMat"))
    
    A0=A
    missing_ind=which(apply(is.na(do.call(cbind,A0)),1,sum)>0) 
    #First imputation
    if(!is.null(missing_ind))
    {
        #TODO adding other initializations
        A=imputeColmeans(A0)
    }else{warnings("No missing data.")}
    
    # Loop up to the convergence of the algorithm
    loop=1
    difference=objective=old=criterion=list()
    old[[1]]=Inf
    continue=TRUE
    rgccaCrit=c()
    while (continue)
    {
        # this step includes a scaling step that should be taken into account
        res.rgccac=rgccac(A=A, scale = scale , init=init, bias = bias, tol = tol, verbose=verbose) 
        rgccaCrit=c(rgccaCrit,res.rgccac$crit[length(res.rgccac$crit)])
        stdevs=res.rgccac$scale_scale
        means=res.rgccac$scale_center
        lambda=res.rgccac$lambda
        J=length(A)
        n=NROW(A[[1]])
        residual_variance=res.rgccac$residual_variance
        # All-inclusive
        if(estimationType=="allInclusive")  # formula 94
        {
            estimation=list()
            nbvarTotal=sum(sapply(A,NCOL))
            nu=matrix(NA,n,J)
           
            Lambda=matrix(0,nbvarTotal,J)
            start_block=c(1,1+cumsum(sapply(A,NCOL)))
            for(j in 1:J)
            {
               Lambda[start_block[j]:(start_block[j+1]-1),j]=res.rgccac$lambda[[j]]
            }
            Phat=res.rgccac$Phat
            Acentered=lapply(A,function(x){return(scale(x,scale=FALSE))})
            Atot=do.call(cbind,Acentered)
            
            Sigma=1/n*t(Atot)%*%Atot
            missing_ind_j=which(apply(is.na(A0[[j]]),1,sum)>0) 
          
            for (j in 1:J)
            {
                stdevj=stdevs[[j]]
                meanj=means[[j]]
                estimation[[j]]=matrix(NA,n,NCOL(A[[j]]))
                for(i in 1:n)
                {
                
                    nu[i,j]= as.vector(Phat%*%t(Lambda)%*%solve(Sigma)%*%matrix(Atot[i,],ncol=1))[j]
                    estimation[[j]][i,]=lambda[[j]]*nu[i,j]
                    missing_ind_j=which(apply(is.na(A0[[j]]),1,sum)>0) 
                    if(i %in% missing_ind_j)
                    {
                        missing_var=which(is.na(A0[[j]][i,]))
                        if(!scale)
                        {
                            A[[j]][i,missing_var]=(meanj+estimation[[j]][i,])[missing_var]
                        }
                        #   if(scale)
                        #   {
                        #       A[[j]][i,missing_var]=(meanj+stdevj*estimation[[j]][i,])[missing_var]
                        #   }
                        
                    }
                }
            }
        }
        # Single-block
        if(estimationType=="singleBlock")  # formula 82
        {
            
            nu=matrix(NA,NROW(A[[1]]),J)
            start_block=c(1,1+cumsum(sapply(A,NCOL)))
            estimation=list()
            for (j in 1:J)
            {
                estimation[[j]]=matrix(NA,NROW(A[[j]]),NCOL(A[[j]]))
                invthetaj= diag(1/residual_variance[start_block[j]:(start_block[j+1]-1)]) # !! Negative residual variance
                # /!\ invthetaj should be defined positive
                if(sum(invthetaj<0)>0)
                {
                    invthetaj[invthetaj<0]=0
                    warning("the matrix was not positive... Negative coefficients were transformed as 0")
                }
                
                missing_ind_j=which(apply(is.na(A0[[j]]),1,sum)>0) 
                stdevj=stdevs[[j]]
                meanj=means[[j]]
                
                for(i in 1:n)
                {
                    nu[i,j]=matrix(lambda[[j]],nrow=1)%*%invthetaj%*%scale(A[[j]],scale=FALSE)[i,]/(1+matrix(lambda[[j]],nrow=1)%*%invthetaj%*%matrix(lambda[[j]],ncol=1))
                    #c\heck: matrix(lambda[[j]],nrow=1)%*%solve(cov(A[[j]]))%*%scale(A[[j]],scale=FALSE)[i,]
                    estimation[[j]][i,]=lambda[[j]]*nu[i,j]
                    
                    if(i %in% missing_ind_j)
                    {
                        missing_var=which(is.na(A0[[j]][i,]))
                        if(!scale)
                        {
                            A[[j]][i,missing_var]=(meanj+estimation[[j]][i,])[missing_var]
                        }
                        #   if(scale)
                        #   {
                        #       A[[j]][i,missing_var]=(meanj+stdevj*estimation[[j]][i,])[missing_var]
                        #   }
                        
                    }
                    
                }
                # head(estimation[[j]])
                # to compare to 
                #if(j==1) print(head(A[[1]]))
            }
        }
        
        concatenedBlocks=do.call(cbind,lapply(A0,function(x) scale(x,scale=FALSE)))
        concatenedXhat=do.call(cbind, estimation)
        indNA=which(is.na(concatenedBlocks),arr.ind=T)
        
        # Criteria of convergence are calculated
        #----------------------------------------
        difference[[loop]] <- concatenedBlocks - concatenedXhat
        difference[[loop]][indNA] <- 0
        objective[[loop]] <-
            sqrt(sum(difference[[loop]] ^ 2) / (
                dim(concatenedBlocks)[1] *
                    dim(concatenedBlocks)[2] - sum(is.na(concatenedBlocks))
            ))
        criterion[[loop]] <- abs(1 - objective[[loop]] / old[[loop]])
        old[[loop + 1]] <- objective[[loop]]
        if (!is.nan(criterion[[loop]])) 
        {
            if (criterion[[loop]] < tol && (loop > 2)) 
            {
                continue <- FALSE
            }
            if (objective[[loop]] < tol)
            {
                continue <- FALSE
                if (verbose) 
                {
                    cat(
                        "The algorithm converged to a stationary point after",
                        loop - 1,
                        "iterations \n"
                    )
                }
            }
        }
        if (loop > ni) 
        {
            continue <- FALSE
            warning(
                paste(
                    "The RGCCA imputation algorithm did not converge after ",
                    ni,
                    " iterations"
                )
            )
        }
        loop <- loop + 1
        print(paste0("iteration in algo:",loop))
    }
    return(list(A=A,obj=unlist(objective),crit=unlist(criterion),rgccaCrit=rgccaCrit))
}