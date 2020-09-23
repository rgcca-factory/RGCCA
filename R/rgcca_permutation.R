#' Tuning RGCCA parameters with permutation
#' 
#' To tune the sparsity coefficient (if the model is sparse) or tau 
#' (otherwise), we observe the deviation between the model and a set of models 
#'  where the lines of each block are permuted. The model with the best 
#' combination of parameters is the one with the highest deviation with the 
#' RGCCA criteria.
#' @inheritParams set_connection
#' @inheritParams bootstrap
#' @inheritParams rgcca
#' @inheritParams plot2D
#' @param par_type A character giving the parameter to tune among "sparsity" or "tau".
#' @param par_length An integer indicating the number of sets of parameters to be tested (if perm.value = NULL). The parameters are uniformly distributed.
#' @param par_value A matrix (n*p, with p the number of blocks and n the number 
#' of combinations to be tested), a vector (of p length) or a numeric value 
#' giving sets of penalties (tau for RGCCA, sparsity for SGCCA) to be tested, 
#' one row by combination. By default, it takes 10 sets between min values (0
#'  for RGCCA and $1/sqrt(ncol)$ for SGCCA) and 1.
#' @param n_run An integer giving the number of permutation tested for each set of constraint
#' @return \item{zstat}{The Z statistic is the difference, for each combination, between the non-permuted R/SGCCA criterion and the mean of the permuted R/SGCCA criteria divided by the standard deviation of the permuted R/SGCCA criteria.}
#' @return \item{bestpenalties}{Penalties giving the best Z-statistic for each blocks}
#' @return \item{permcrit}{A matrix of R/SGCCA criteria for each combination and each permutation}
#' @return \item{means}{Mean of the permutated R/SGCCA criteria for each combination}
#' @return \item{sds}{Standard deviation of the permutated R/SGCCA criteria for each combination}
#' @return \item{crit}{R/SGCCA criterion for each combination}
#' @return \item{pval}{The p-value is the fraction of the permuted R/SGCCA criteria, for each combination, that is greater than or equal to the non-permuted R/SGCCA criterion}
#' @return \item{penalties}{A matrix giving, for each blocks, the penalty combinations (tau or sparsity)}
#' @examples
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' res = rgcca_permutation(blocks, n_run = 5, n_cores = 1)
#' rgcca_permutation(blocks, par_type = "sparsity", par_value = 0.8, n_run = 2,
#'  n_cores = 1)
#' rgcca_permutation(blocks, par_type = "sparsity", par_value = c(0.6, 0.75, 0.5), 
#' n_run = 2, n_cores = 1)
#' rgcca_permutation(blocks, par_type = "sparsity", 
#' par_value = matrix(c(0.6, 0.75, 0.5), 3, 3, byrow = TRUE),
#'  n_run = 2, n_cores = 1)
#' rgcca_permutation(blocks, par_type = "tau", par_value = 0.8, n_run= 2, 
#' n_cores = 1)
#' rgcca_permutation(blocks, par_type = "tau", par_value = c(0.6, 0.75, 0.5),
#'  n_run = 2, n_cores = 1)
#' rgcca_permutation(blocks, par_type = "tau", par_value = 
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
    tol = 1e-8,
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
    check_integer("par_length", n_run)
    check_integer("par_value", n_run, min = 0)
    check_integer("n_cores", n_cores, min = 0)
    match.arg(par_type, c("tau", "sparsity"))
    if (!is.null(parallelization))
        check_boolean("parallelization", parallelization)
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
    set_penalty <- function () 
    {
        # Selecting the minimal value
        if (par_type == "sparsity") 
        {
            if (!tolower(type) %in% c("spls", "spca", "sgcca"))
                warning("The sparsity is chosen but the analyse was not sparse. By default, a SGCCA will be performed.")
            type <<- "sgcca"
            min_spars <<- sapply(ncols, function(x) 1 / sqrt(x))
        }else
        {
            if (tolower(type) %in% c("spls", "spca", "sgcca"))
                warning("The tau is chosen but the analyse is sparse. By default, a RGCCA will be performed.")
            type <<- "rgcca"
            min_spars <<- sapply(ncols, function(x) 0)
        }
        
        
        if (is.null(par_value))
            par_value <- set_spars()
        else if (class(par_value) %in% c("data.frame", "matrix")) #when a matrix is entered
        {
            if(par_type=="tau")
            {
                par_value <- t(sapply(seq(NROW(par_value)), function(x) check_tau(par_value[x, ], blocks, type = type,superblock=superblock)))
            }
                     
        }
        else #when a vector is entered
        {
            if (any(par_value < min_spars))
                stop_rgcca(paste0("par_value should be upper than : ", paste0(round(min_spars, 2), collapse = ",")))
             if(par_type=="tau")
             {
              par_value <- check_tau(par_value, blocks, type = type,superblock=superblock)
              par_value <- set_spars(max = par_value)
             }
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
