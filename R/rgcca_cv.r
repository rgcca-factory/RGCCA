#' rgcca_cv
#'
#' This function is dedicated to supervised approaches (with response parameter in rgcca not NULL). By default, it is the last block in blocks.
#' It allows to test a grid of parameters (sparsity, ncomp or tau) with cross-validation approaches. 
#' This cross-validation is based on RMSE for quantitative response. 
#' @inheritParams rgcca_crossvalidation
#' @inheritParams rgcca
#' @param parallelization if TRUE parallelization is run, if FALSE, no parallelisation is run. If NULL (default) parallelization is always used except for Windows in case of length(nperm)<10
#' @param par "sparsity", "tau" or "ncomp".
#' @param par_value Grid of values to be tested. #TODO
#' @param type_cv  type of crossvalidation. Default to "regression" #TODO
#' @param n_cv number of cross-validation. Default to 1 (with kfold option and 5)
#' @param one_value_per_cv If TRUE, the k values obtained for each k-fold steps are averaged 
#' @param rgcca_res a result of rgcca (from whom all parameters will be imported)
#' @export
#' @examples
#' data("Russett")
#' blocks <- list(
#'     agriculture = Russett[, seq(3)],
#'     industry = Russett[, 4:5],
#'     politic = Russett[, 6:11])
#'     res=rgcca_cv(blocks,response=3, type="rgcca",par="tau",par_value=c(0,0.2,0.3),n_cv=1,n_cores=1)
#'     
#'@importFrom utils txtProgressBar setTxtProgressBar
rgcca_cv=function( blocks,
          type = "rgcca",
          response=NULL,
          par = "tau",
          par_value = NULL,
          validation = "kfold",
          type_cv = "regression",
          fit = "lm",
          k=5,
          n_cv = 1,
          one_value_per_cv=FALSE,
          n_cores = parallel::detectCores() - 1,
          quiet = TRUE,
          superblock=FALSE,
          scale=TRUE,
          scale_block=TRUE,
          tol=1e-6,
          scheme="factorial",
          method="nipals",
          rgcca_res=NULL,
          parallelization=NULL,
          tau=rep(1,length(blocks)),
          ncomp=rep(1,length(blocks)),
          sparsity=rep(1,length(blocks)),
          init="svd",
          bias=TRUE,
          ...)
{
    if(!missing(blocks)&class(blocks)=="rgcca"){rgcca_res=blocks}
    if(class(rgcca_res)=="rgcca")
    {
        message("All parameters were imported by a rgcca object provided in the blocks/rgcca_res parameter")
        scale_block=rgcca_res$call$scale_block
        scale=rgcca_res$call$scale
        scheme=rgcca_res$call$scheme
        response=rgcca_res$call$response
        tol=rgcca_res$call$tol
        method=rgcca_res$call$method
        bias=rgcca_res$call$bias
        blocks<-rgcca_res$call$raw
        superblock=rgcca_res$call$superblock
        connection=rgcca_res$call$connection
        tau=rgcca_res$call$tau
        ncomp=rgcca_res$call$ncomp
        sparsity=rgcca_res$call$sparsity
    }
    
    if(is.null(response)){ stop("response is required for rgcca_cv (it is an integer comprised between 1 and the number of blocks) ")}
    if(superblock){stop_rgcca("Cross-validation is only possible without superblock")}
    if(validation=="loo")
    {
        k=dim(blocks[[1]])[1]
        if(n_cv!=1)
        {
            cat("n_cv value was replaced by 1 (is not relevant for loo option)")
        };n_cv=1
    }

    check_integer("n_cores", n_cores, 0)
    match.arg(par, c("tau", "sparsity", "ncomp"))
    min_spars <- NULL

    if (length(blocks) < 1)
        stop_rgcca("Crossvalidation required a number of blocks larger than 1.")
    
    ncols <- sapply(blocks, NCOL)
    
    set_spars <- function(max = 1) {
        if (length(max) == 1)
            f <- quote(max)
        else
            f <- quote(max[x])
        sapply(seq(min_spars), function(x) seq(eval(f), min_spars[x], len = 10))
    }
    set_penalty <- function () {
        if(par == "sparsity"){
            if(type!="sgcca"){cat("As par=='sparsity', the type parameter was replaced by 'sgcca'")}
            type <- "sgcca"
            min_spars <<- sapply(ncols, function(x) 1 / sqrt(x))
        }else{
            if(type=="sgcca"){cat("As par!='sparsity', the type parameter was replaced by 'rgcca'")}
            type <- "rgcca"
            min_spars <<- sapply(ncols, function(x) 0)
        }
        
        if (is.null(par_value))
            par_value <- set_spars()
        else if (class(par_value) %in% c("data.frame", "matrix"))
            par_value <- t(sapply(seq(NROW(par_value)), function(x) check_tau(par_value[x, ], blocks, type = type)))
        else{
            if (any(par_value < min_spars))
                stop_rgcca(paste0("par_value should be upper than : ", paste0(round(min_spars, 2), collapse = ",")))
            par_value <- check_tau(par_value, blocks, type = type)
            par_value <- set_spars(max = par_value)
        }
  
        colnames(par_value) <- names(blocks)
        return(list(par, par_value))
    }
    
    switch(
        par,
        "ncomp" = {
            if (!class(par_value) %in% c("data.frame", "matrix")) {
                if (is.null(par_value) || any(par_value > ncols)) {
                    ncols[ncols > 5] <- 5
                    par_value <- ncols
                }else
                    par_value <- check_ncomp(par_value, blocks)
                par_value <- lapply(par_value, function(x) seq(x))
                par_value <- expand.grid(par_value)
            }else
                par_value <- t(sapply(seq(NROW(par_value)), function(x) check_ncomp(par_value[x, ], blocks, 1)))
            par <- list(par, par_value)
        },
        "sparsity" = par <- set_penalty(),
        "tau" = par <- set_penalty()
    )
  
    print(paste("Cross-validation for", par[[1]], "is in progression..."))
    pb <- txtProgressBar(max=dim(par[[2]])[1])
    n_rep=ifelse(one_value_per_cv,n_cv,n_cv*k)
    res=matrix(NA,dim(par[[2]])[1],n_cv*k);rownames(res)=apply(round(par[[2]],digits=2),1,paste,collapse="-");
    for(i in 1:dim(par[[2]])[1])
        {
            if(par[[1]]=="ncomp")
            {
                rgcca_res=rgcca(blocks=blocks, type=type,response=response,ncomp=par[[2]][i,],superblock=superblock,scale=scale,scale_block=scale_block,scheme=scheme,tol=tol,method=method,tau=tau, sparsity=sparsity,bias=bias,init=init)
            }
            if(par[[1]]=="sparsity")
            {
                rgcca_res=rgcca(blocks=blocks, type="sgcca",response=response,sparsity=par[[2]][i,],superblock=superblock,scale=scale,scale_block=scale_block,scheme=scheme,tol=tol,method=method, ncomp=ncomp,bias=bias,init=init)
            }
            if(par[[1]]=="tau")
            {
                rgcca_res=rgcca(blocks=blocks, type=type,response=response,tau=par[[2]][i,],superblock=superblock,scale=scale,scale_block=scale_block,scheme=scheme,tol=tol,method=method, ncomp=ncomp,bias=bias,init=init)
            }
 
            res_i=c()
            for(n in 1:n_cv)
            {
                if(one_value_per_cv)
                {
                    res_cv_n=rgcca_crossvalidation(
                        rgcca_res,
                        validation = validation,
                        model = type_cv,
                        fit = fit,
                        # new_scaled = TRUE,
                        k = k,
                        n_cores =n_cores,
                        parallelization=parallelization
                    )
                    if(!is.null(res_cv_n))
                    {
                        res_i=c(res_i,res_cv_n$scores)
                    }
                    else
                    {
                        print("One cross-validation was stopped (at least one variable had null standard deviations)")
                    }
               
                }
                else
                {
                    res_cv_n=rgcca_crossvalidation(
                        rgcca_res,
                        validation = validation,
                        model = type_cv,
                        fit = fit,
                        #new_scaled = TRUE,
                        k = k,
                        n_cores =n_cores,
                        parallelization=parallelization)
                    
                    res_i= c(res_i,res_cv_n$list_scores)
                }
              
            }
  
            res[i,]=as.numeric(res_i)
        
            #Sys.sleep(0.5); 
            setTxtProgressBar(pb, i)
            mat_cval=res
            rownames(mat_cval)=1:NROW(mat_cval)
        }
    
    Sys.sleep(1)
    call=list(n_cv=n_cv,
              response=response,
              par = par,
              par_value = par_value,
              validation = validation,
              type_cv = type_cv,
              fit = fit,
              k=k,
              one_value_per_cv=one_value_per_cv,
              superblock=FALSE,
              scale=scale,
              scale_block=scale_block,
              tol=tol,
              scheme=scheme,
              method=method,
              blocks=blocks
    )
    par2=par[[2]]
    rownames(par2) = 1:NROW(par2)
    colnames(par2)=names(blocks)
    
    res2=list(cv=mat_cval,call=call,bestpenalties=par[[2]][which.min(apply(mat_cval,1,mean)),],penalties=par2)
    class(res2)="cval"
    return(res2)
    
}
