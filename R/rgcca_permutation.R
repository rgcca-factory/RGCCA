#' Run through a set of constraint parameters sparsitys to select the best with permutation
#' Only one component per block for the time being
#' 
#' @inheritParams set_connection
#' @inheritParams bootstrap
#' @inheritParams rgcca
#' @param p_spars A matrix, a vector or an integer containing sets of constraint 
#' variables, one row by combination. By default, sgcca.permute takes 10 sets between 
#' min values ($1/sqrt(ncol)$) and 1
#' @param p_ncomp A matrix, a vector or an integer containing sets of number of 
#' components, one row by set. By default, sgcca.permute takes as many 
#' combinations as the maximum number of columns in each block
#' @param nperm Number of permutation tested for each set of constraint
#' @return A list containing :
#' @return \item{pval}{Pvalue}
#' @return \item{zstat}{Statistic Z}
#' @return \item{bestpenalties}{Penalties corresponding to the best Z-statistic}
#' @return \item{permcrit}{RGCCA criteria obtained with permutation set}
#' @return \item{crit}{ RGCCA criterion for the original dataset}
#' @examples
#' data("Russett")
#' A = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_permutation(A, nperm = 2, n_cores = 1)
#' rgcca_permutation(A, p_ncomp = TRUE, nperm = 2, n_cores = 1)
#' rgcca_permutation(A, p_spars = 0.8, nperm = 2, n_cores = 1)
#' rgcca_permutation(A, p_spars = c(0.6, 0.75, 0.5), nperm = 2, n_cores = 1)
#' rgcca_permutation(A, p_spars = matrix(c(0.6, 0.75, 0.5), 3, 3, byrow = TRUE), nperm = 2, n_cores = 1)
#' rgcca_permutation(A, p_ncomp = 2, nperm = 2, n_cores = 1)
#' rgcca_permutation(A, p_ncomp = c(2,2,3), nperm = 2, n_cores = 1)
#' rgcca_permutation(A, p_ncomp = matrix(c(2,2,3), 3, 3, byrow = TRUE), nperm = 2, n_cores = 1)
#' @export
rgcca_permutation <- function(
    blocks,
    type = "rgcca",
    p_spars = TRUE,
    p_ncomp = FALSE,
    nperm = 20,
    n_cores = parallel::detectCores() - 1,
    quiet=TRUE,
    ...) {
    call=list(type=type,p_spars=p_spars,p_ncomp=p_ncomp,nperm=nperm,quiet=quiet)
    check_integer("nperm", nperm)
    check_integer("n_cores", n_cores, 0)

    if (any(p_ncomp == FALSE) && any(p_spars == FALSE))
        stop("Select one parameter among 'p_spars' or 'p_ncomp' to optimize. By default, p_spars is selected.")

    if (length(blocks) < 1)
        stop("Permutation required a number of blocks larger than 1.")

    ncols <- sapply(blocks, NCOL)
    min_spars <- sapply(ncols, function(x) 1 / sqrt(x))

    set_spars <- function(max = 1) {
        if (length(max) == 1)
            f <- quote(max)
        else
            f <- quote(max[x])
        sapply(seq(min_spars), function(x) seq(eval(f), min_spars[x], len = 10))
    }

    if (!any(p_ncomp == FALSE)) {
        if (!class(p_ncomp) %in% c("data.frame", "matrix")) {
            if (isTRUE(p_ncomp) || any(p_ncomp > ncols)) {
                ncols[ncols > 5] <- 5
                p_ncomp <- ncols
            }else
                p_ncomp <- check_ncomp(p_ncomp, blocks)
            p_ncomp <- lapply(p_ncomp, function(x) seq(x))
            p_ncomp <- expand.grid(p_ncomp)
        }else
            p_ncomp <- t(sapply(seq(NROW(p_ncomp)), function(x) check_ncomp(p_ncomp[x, ], blocks, 1)))
        par <- list("ncomp", p_ncomp)
    }

    if (!any(p_spars == FALSE)) {
        if (identical(p_spars, TRUE))
            p_spars <- set_spars()
        else if (class(p_spars) %in% c("data.frame", "matrix"))
            p_spars <- t(sapply(seq(NROW(p_spars)), function(x) check_tau(p_spars[x, ], blocks, type = "sgcca")))
        else{
            if (any(p_spars < min_spars))
                stop(paste0("p_spars should be upper than 1 / sqrt(NCOL(blocks)) : ", paste0(round(min_spars, 2), collapse = ",")))
            p_spars <- check_tau(p_spars, blocks, type = "sgcca")
            p_spars <- set_spars(max = p_spars)
        }

        colnames(p_spars) <- names(blocks)
        par <- list("sparsity", p_spars)
        type <- "sgcca"
    }

    crits <- rgcca_permutation_k(
        blocks,
        par = par,
        perm = FALSE,
        type = type,
        n_cores = 1,
        quiet=quiet,
        ...
    )

    cat("Permutation in progress...")
    # To uncomment when it is tested
    # if (Sys.info()["sysname"] == "Windows") {
    # 
    #     e <- environment()
    #     cl <- parallel::makeCluster(n_cores)
    # 
    #     parallel::clusterExport(
    #         cl,
    #         c(
    #             "blocks",
    #             "p_spars",
    #             "nperm",
    #             "C",
    #             "ncomp",
    #             "scheme",
    #             "out",
    #             "crit",
    #             "crits",
    #             "tol"
    #         ),
    #         envir = e
    #     )

        # /!\ To be uncomment (packaging)
        # parallel::clusterEvalQ(cl, library(devtools))

        # tryCatch({
        #     parallel::clusterEvalQ(cl, load_all("RGCCA/R/."))
        # }, error = function(e) {
        #     warning("error : probably an issue with the localisation of RGCCA functions")
        # })
# /!\ End to be uncomment (packaging)
        # Close cluster even if there is an error or a warning with rgcca_permutation_k
    #     permcrit <- tryCatch({
    #         parallel::parSapply(cl, seq(nperm), function(x)
    #             rgcca_permutation_k(
    #                 blocks = blocks,
    #                 par = par,
    #                 type = type,
    #                 ...
    #             ))
    #     }, error = function(e) {
    #         warning("an error occured with rgcca_permutation_k")
    #         return(NULL)
    #     })
    # 
    #     parallel::stopCluster(cl)
    # 
    #     if (is.null(permcrit))
    #         return(NULL)
    # 
    # } else {
        permcrit <- as.matrix(simplify2array(parallel::mclapply(
            seq(nperm),
            function(x){
                res <- rgcca_permutation_k(
                    blocks = blocks,
                    par = par,
                    type = type,
                    n_cores = 1,
                    quiet=quiet,
                    ...
                )
                return(res)
                },
            mc.cores = n_cores)))
    # }

    cat("OK.\n", append = TRUE)

    par <- par[[2]]

    pvals <- sapply(
        seq(NROW(par)),
        function(i)
            mean(permcrit[i, ] >= crits[i]))
    zs <- sapply(
        seq(NROW(par)),
        function(i){
            z <- (crits[i] - mean(permcrit[i, ])) / (sd(permcrit[i, ]))
            if (is.na(z) || z == "Inf")
                z <- 0
            return(z)
        })

    structure(
        list(
            call=call,
            pvals = pvals,
            zstat = zs,
            bestpenalties = par[which.max(zs), ],
            permcrit = permcrit,
            crit = crits,
            penalties = par
        ),
        class = "permutation"
    )
}
