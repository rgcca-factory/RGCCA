#' Tuning RGCCA parameters
#' 
#' Run through a set of parameters (sparsity or number of selected components) with permutation to select the one maximizing RGCCA criterion 
#' The sparsity parameter is tuned with only one component per block.
#' @inheritParams set_connection
#' @inheritParams bootstrap
#' @inheritParams rgcca
#' @param perm.par "sparsity","tau" or "ncomp".
#' @param perm.value  If perm.par="sparsity", a matrix, a vector or an integer containing sets of constraint 
#' variables to be tested, one row by combination. By default, sgcca.permute takes 10 sets between 
#' min values ($1/sqrt(ncol)$) and 1. If perm.par="ncomp", a matrix, a vector or an integer containing sets of number of 
#' components, one row by set. By default, sgcca.permute takes as many 
#' combinations as the maximum number of columns in each block. If perm.par="tau",... #TODO
#' @param nperm Number of permutation tested for each set of constraint
#' @return A object permutation, which is a list containing :
#' @return \item{pval}{Pvalue}
#' @return \item{zstat}{Statistic Z}
#' @return \item{bestpenalties}{Penalties corresponding to the best Z-statistic}
#' @return \item{permcrit}{RGCCA criteria obtained with permutation set}
#' @return \item{crit}{ RGCCA criterion for the original dataset}
#' @examples
#' data("Russett")
#' A = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' res = rgcca_permutation(A, nperm = 5, n_cores = 1)
#' rgcca_permutation(A, perm.par = "sparsity", perm.value = 0.8, nperm = 2,
#'  n_cores = 1)
#' rgcca_permutation(A, perm.par = "sparsity", perm.value = c(0.6, 0.75, 0.5), 
#' nperm = 2, n_cores = 1)
#' rgcca_permutation(A, perm.par = "sparsity", 
#' perm.value = matrix(c(0.6, 0.75, 0.5), 3, 3, byrow = TRUE),
#'  nperm = 2, n_cores = 1)
#' rgcca_permutation(A, perm.par = "tau", perm.value = 0.8, nperm = 2, 
#' n_cores = 1)
#' rgcca_permutation(A, perm.par = "tau", perm.value = c(0.6, 0.75, 0.5),
#'  nperm = 2, n_cores = 1)
#' rgcca_permutation(A, perm.par = "tau", perm.value = 
#' matrix(c(0.6, 0.75, 0.5), 3, 3, byrow = TRUE),  nperm = 2, n_cores = 1)

#' print(res)
#' @export
rgcca_permutation <- function(
    blocks,
    perm.par = "tau",
    perm.value = NULL,
    nperm = 20,
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
    ...) 
    {
    if(class(blocks)=="permutation")
    {
        blocks<-blocks$call$blocks
    }
    # call <- as.list(formals(rgcca_permutation))
    call=list(type=type, perm.par = perm.par, perm.value = perm.value, nperm=nperm, quiet=quiet,connection=connection,method=method,tol=tol,scheme=scheme,scale=scale,scale_block=scale_block,blocks=blocks,superblock=superblock)
    check_integer("nperm", nperm)
    check_integer("n_cores", n_cores, 0)
    match.arg(perm.par, c("tau", "sparsity"))
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
        sapply(seq(min_spars), function(x) seq(eval(f), min_spars[x], len = 10))
    }
    set_penalty <- function () {

        if(perm.par == "sparsity"){
            if(type!="sgcca"){cat("As par=='sparsity', the type parameter was replaced by 'sgcca'\n")}
            type <<- "sgcca"
            min_spars <<- sapply(ncols, function(x) 1 / sqrt(x))
        }else{
            if(type!="rgcca"){cat("As par!='sparsity', the type parameter was replaced by 'rgcca'\n")}
            type <<- "rgcca"
            min_spars <<- sapply(ncols, function(x) 0)
        }

        if (is.null(perm.value))
            perm.value <- set_spars()
        else if (class(perm.value) %in% c("data.frame", "matrix"))
        {
            if(perm.par=="tau")
            {
                perm.value <- t(sapply(seq(NROW(perm.value)), function(x) check_tau(perm.value[x, ], blocks, type = type,superblock=superblock)))
            }
                     
        }
                else{
            if (any(perm.value < min_spars))
                stop_rgcca(paste0("perm.value should be upper than : ", paste0(round(min_spars, 2), collapse = ",")))
             if(perm.par=="tau"){perm.value <- check_tau(perm.value, blocks, type = type,superblock=superblock)}
           if(perm.par=="sparsity"){perm.value <- set_spars(max = perm.value)}
        }

        if(superblock){coln=c(names(blocks),"superblock")}
        else{coln=names(blocks)}
        if(is.null(dim(perm.value))){perm.value=matrix(perm.value,nrow=1)}
        colnames(perm.value) <- coln
        return(list(perm.par, perm.value))
    }

    switch(
        perm.par,
    #     "ncomp" = {
    #     if (!class(perm.value) %in% c("data.frame", "matrix")) {
    #         if (is.null(perm.value) || any(perm.value > ncols)) {
    #             ncols[ncols > 5] <- 5
    #             perm.value <- ncols
    #         }else
    #             perm.value <- check_ncomp(perm.value, blocks)
    #         perm.value <- lapply(perm.value, function(x) seq(x))
    #         perm.value <- expand.grid(perm.value)
    #     }else
    #         perm.value <- t(sapply(seq(NROW(perm.value)), function(x) check_ncomp(perm.value[x, ], blocks, 1)))
    #     par <- list(perm.par, perm.value)
    # },
    "sparsity" = par <- set_penalty(),
    "tau" = par <- set_penalty()
    )

    message("Permutation in progress...\n", appendLF = FALSE)

    varlist <- c(ls(getNamespace("RGCCA")))
    # get the parameter dot-dot-dot
    args_values <- list(...)
    args_names <- names(args_values)
    n <- args_values
    if (!is.null(n))
        n <- seq(length(args_values))
    for (i in n) {
        if (!is.null(args_names[i])) {
            # dynamically asssign these values
            assign(args_names[i], args_values[[i]])
            # send them to the clusters to parallelize
            varlist <- c(varlist, args_names[i])
            # without this procedure rgcca_crossvalidation(rgcca_res, blocks = blocks2)
            # or rgcca_crossvalidation(rgcca_res, blocks = lapply(blocks, scale)
            # does not work.
        }
    }

    pb <- txtProgressBar(max=dim(par[[2]])[1])
    crits=means=sds=rep(NA,nrow(par[[2]]))
    permcrit=matrix(NA,nrow(par[[2]]),nperm)
    for(i in 1:nrow(par[[2]]))
    {
      
        crits[i] <- rgcca_permutation_k(
            blocks,
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
            ...
        )
          
        if(nperm>20)
        {
            res<- parallelize(
                varlist,
                seq(nperm), 
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
                        ...
                    ),
                n_cores = n_cores,
                envir = environment(),
                applyFunc = "parSapply"
            )   
        }
        if(nperm<=20)
        {
            res=sapply(seq(nperm),function(k)  {rgcca_permutation_k(
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
                ...
            )})
        }
         
        
        permcrit[i,] =res
        means[i]=mean(permcrit[i,],na.rm=T)
        sds[i]=sd(permcrit[i,],na.rm=T)
        setTxtProgressBar(pb, i)
    }
 

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
