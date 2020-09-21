#' Tuning RGCCA parameters
#' 
#' Run through a set of parameters (sparsity or number of selected components) with permutation to select the one maximizing RGCCA criterion 
#' The sparsity parameter is tuned with only one component per block.
#' @inheritParams set_connection
#' @inheritParams bootstrap
#' @inheritParams rgcca
#' @inheritParams plot2D
#' @param par_type A character giving the parameter to tune among "sparsity" or "tau".
#' @param par_length If par_value = NULL, an integer indicating the number of sets of parameters to be tested. The parameters are uniformly distributed.
#' @param par_value  If par_type="sparsity", a matrix, a vector or an integer containing sets of constraint 
#' variables to be tested, one row by combination. By default, sgcca.permute takes 10 sets between 
#' min values ($1/sqrt(ncol)$) and 1. If par_type="ncomp", a matrix, a vector or an integer containing sets of number of 
#' components, one row by set. By default, sgcca.permute takes as many 
#' combinations as the maximum number of columns in each block. If par_type="tau",... #TODO
#' @param n_run An integer giving the number of permutation tested for each set of constraint
#' @return A permutation object, which is a list containing :
#' @return \item{pval}{The p-value}
#' @return \item{zstat}{The Z statistic}
#' @return \item{bestpenalties}{Penalties corresponding to the best Z-statistic}
#' @return \item{permcrit}{RGCCA criteria obtained with the permutation set}
#' @return \item{crit}{A vector of integer giving the values of the RGCCA criteria across iterations.}
#' @examples
#' data("Russett")
#' A = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' res = rgcca_permutation(A, n_run = 5, n_cores = 1)
#' rgcca_permutation(A, par_type = "sparsity", par_value = 0.8, n_run = 2,
#'  n_cores = 1)
#' rgcca_permutation(A, par_type = "sparsity", par_value = c(0.6, 0.75, 0.5), 
#' n_run = 2, n_cores = 1)
#' rgcca_permutation(A, par_type = "sparsity", 
#' par_value = matrix(c(0.6, 0.75, 0.5), 3, 3, byrow = TRUE),
#'  n_run = 2, n_cores = 1)
#' rgcca_permutation(A, par_type = "tau", par_value = 0.8, n_run = 2, 
#' n_cores = 1)
#' rgcca_permutation(A, par_type = "tau", par_value = c(0.6, 0.75, 0.5),
#'  n_run = 2, n_cores = 1)
#' rgcca_permutation(A, par_type = "tau", par_value = 
#' matrix(c(0.6, 0.75, 0.5), 3, 3, byrow = TRUE),  n_run = 2, n_cores = 1)
#' print(res)
#' @export
rgcca_permutation <- function(
    blocks,
    par_type = "tau",
    par_value = NULL,
    par_length=10,
    n_run = 20,
    n_cores = parallel::detectCores() - 1,
    quiet = TRUE,
    type = "rgcca",
    scale = TRUE,
    scale_block = TRUE,
    connection = matrix(1,length(blocks),length(blocks)) - diag(length(blocks)),
    scheme = "factorial",
    ncomp = rep(1, length(blocks)),
    tau = rep(1, length(blocks)),
    sparsity = rep(1, length(blocks)),
    init = "svd",
    bias = TRUE,
    tol = 1e-08,
    response = NULL,
    superblock = FALSE,
    method = "nipals",
    rgcca_res=NULL,
    parallelization=NULL
    ) 
    {
    if(!missing(blocks)&class(blocks)=="rgcca"){rgcca_res=blocks}
    if(class(rgcca_res)=="rgcca")
    {
        message("All parameters were imported from a rgcca object.")
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
   
    # call <- as.list(formals(rgcca_permutation))
    call=list(type=type, par_type = par_type, par_value = par_value, n_run=n_run, quiet=quiet,connection=connection,method=method,tol=tol,scheme=scheme,scale=scale,scale_block=scale_block,blocks=blocks,superblock=superblock)
    check_integer("n_run", n_run)
    check_integer("n_cores", n_cores, 0)
    match.arg(par_type, c("tau", "sparsity"))
    min_spars <- NULL

    if (length(blocks) < 1)
        stop_rgcca("Permutation required a number of blocks larger than 1.\n")

    if(!superblock)
    {
        ncols <- sapply(blocks, NCOL)
    }
    else
    {
        ncol_block=sapply(blocks, NCOL)
        ncols <- c(ncol_block,sum(ncol_block))
        names(ncols)=c(names(ncol_block),"superblock")
    }
 
    set_spars <- function(max = 1) {
        if (length(max) == 1)
            f <- quote(max)
        else
            f <- quote(max[x])
        sapply(seq(min_spars), function(x) seq(eval(f), min_spars[x], len = par_length))
    }
    set_penalty <- function () {
        if (par_type == "sparsity") {
            if (!tolower(type) %in% c("spls", "spca", "sgcca"))
                warning("The sparsity is chosen but the analyse was not sparse. By default, a SGCCA will be performed.")
            type <<- "sgcca"
            min_spars <<- sapply(ncols, function(x) 1 / sqrt(x))
        }else{
            if (tolower(type) %in% c("spls", "spca", "sgcca"))
                warning("The tau is chosen but the analyse is sparse. By default, a RGCCA will be performed.")
            type <<- "rgcca"
            min_spars <<- sapply(ncols, function(x) 0)
        }

        if (is.null(par_value))
            par_value <- set_spars()
        else if (class(par_value) %in% c("data.frame", "matrix"))
        {
            if(par_type=="tau")
            {
                par_value <- t(sapply(seq(NROW(par_value)), function(x) check_tau(par_value[x, ], blocks, type = type,superblock=superblock)))
            }
                     
        }
                else{
            if (any(par_value < min_spars))
                stop_rgcca(paste0("par_value should be upper than : ", paste0(round(min_spars, 2), collapse = ",")))
             if(par_type=="tau"){par_value <- check_tau(par_value, blocks, type = type,superblock=superblock)}
           if(par_type=="sparsity"){par_value <- set_spars(max = par_value)}
        }

        if(superblock){coln=c(names(blocks),"superblock")}
        else{coln=names(blocks)}
        if(is.null(dim(par_value))){par_value=matrix(par_value,nrow=1)}
        colnames(par_value) <- coln
        return(list(par_type, par_value))
    }

    switch(
        par_type,
    #     "ncomp" = {
    #     if (!class(par_value) %in% c("data.frame", "matrix")) {
    #         if (is.null(par_value) || any(par_value > ncols)) {
    #             ncols[ncols > 5] <- 5
    #             par_value <- ncols
    #         }else
    #             par_value <- check_ncomp(par_value, blocks)
    #         par_value <- lapply(par_value, function(x) seq(x))
    #         par_value <- expand.grid(par_value)
    #     }else
    #         par_value <- t(sapply(seq(NROW(par_value)), function(x) check_ncomp(par_value[x, ], blocks, 1)))
    #     par <- list(par_type, par_value)
    # },
    "sparsity" = par <- set_penalty(),
    "tau" = par <- set_penalty()
    )

    message("Permutation in progress...\n", appendLF = FALSE)

    varlist <- c(ls(getNamespace("RGCCA")))
    # get the parameter dot-dot-dot
    # args_values <- list(...)
    # args_names <- names(args_values)
    # n <- args_values
    # if (!is.null(n))
    #     n <- seq(length(args_values))
    # for (i in n) {
    #     if (!is.null(args_names[i])) {
    #         # dynamically asssign these values
    #         assign(args_names[i], args_values[[i]])
    #         # send them to the clusters to parallelize
    #         varlist <- c(varlist, args_names[i])
    #         # without this procedure rgcca_crossvalidation(rgcca_res, blocks = blocks2)
    #         # or rgcca_crossvalidation(rgcca_res, blocks = lapply(blocks, scale)
    #         # does not work.
    #     }
    # }

    pb <- txtProgressBar(max=dim(par[[2]])[1])
    crits=means=sds=rep(NA,nrow(par[[2]]))
    permcrit=matrix(NA,nrow(par[[2]]),n_run)
    for(i in 1:nrow(par[[2]]))
    {
      
        crits[i] <- rgcca_permutation_k(
            blocks,
            connection=connection,
            par = par[[1]],
            par_value=par[[2]][i,],
            perm = FALSE,
            type = type,
            n_cores = 1,
            quiet=quiet,
            superblock=superblock,
            scale=scale,
            scale_block=scale_block,
            scheme=scheme,
            tol=tol,
            ncomp=ncomp,
            sparsity=sparsity
        )
          
      
            res<- parallelize(
                varlist,
                seq(n_run), 
                function(x)
                    rgcca_permutation_k(
                        blocks = blocks,
                        par = par[[1]],
                        par_value=par[[2]][i,],
                        type = type,
                        n_cores = n_cores,
                        quiet = quiet,
                        superblock=superblock,
                        scheme=scheme,
                        tol=tol,
                        scale=scale,
                        scale_block=scale_block,
                        connection=connection,
                        method=method
                    ),
                n_cores = n_cores,
                envir = environment(),
                applyFunc = "parSapply",
                parallelization=parallelization
            )   
       
        
        permcrit[i,] =res
        means[i]=mean(permcrit[i,],na.rm=T)
        sds[i]=sd(permcrit[i,],na.rm=T)
        setTxtProgressBar(pb, i)
    }
    cat("\n")
 

    par <- par[[2]]
    pvals <- sapply(
        seq(NROW(par)),
        function(k)
            mean(permcrit[k, ] >= crits[k]))
    zs <- sapply(
        seq(NROW(par)),
        function(k){
            z <- (crits[k] - mean(permcrit[k, ])) / (sd(permcrit[k, ]))
            if (is.na(z) || z == "Inf")
                z <- 0
            return(z)
        })
    rownames(par) = 1:NROW(par)
    if(superblock){coln=c(names(blocks),"superblock")}
    else{coln=names(blocks)}
    colnames(par)=coln
    structure(
        list(
            call=call,
            zstat = zs,
            bestpenalties = par[which.max(zs), ],
            permcrit = permcrit,
            means=means,
            sds=sds,
            crit = crits,
            pvals = pvals,
            penalties = par
        ),
        class = "permutation"
    )
}
