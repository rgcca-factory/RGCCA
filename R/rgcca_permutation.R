#' Run through a set of constraint parameters c1s to select the best with permutation
#' Only one component per block for the time being
#' 
#' @inheritParams set_connection
#' @inheritParams bootstrap
#' @inheritParams rgcca.analyze
#' @param p_c1 A matrix, a vector or an integer containing sets of constraint 
#' variables, one row by set. By default, sgcca.permute takes 10 sets between 
#' min values ($1/sqrt(ncol)$) and 1
#' @param p_ncomp A matrix, a vector or an integer containing sets of number of 
#' components, one row by set. By default, sgcca.permute takes as many 
#' combinations as the maximum number of columns in each block
#' @param nperm Number of permutation tested for each set of constraint
#' @param ... Others RGCCA parameters #TODO
#' @return A list containing :
#' @return \item{pval}{Pvalue}
#' @return \item{zstat}{Statistic Z}
#' @return \item{bestpenalties}{Penalties corresponding to the best Z-statistic}
#' @return \item{permcrit}{RGCCA criteria obtained with permutation set}
#' @return \item{crit}{ RGCCA criterion for the original dataset}
#' data("Russett")
#' A = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_permutation(A, nperm = 2, n_cores = 1)
#' rgcca_permutation(A, p_c1 = TRUE, nperm = 2, n_cores = 1)
#' rgcca_permutation(A, p_c1 = 0.8, nperm = 2, n_cores = 1)
#' rgcca_permutation(A, p_c1 = c(0.6, 0.75, 0.5), nperm = 2, n_cores = 1)
#' rgcca_permutation(A, p_c1 = matrix(c(0.6, 0.75, 0.5), 3, 3, byrow = T), nperm = 2, n_cores = 1)
#' rgcca_permutation(A, p_ncomp = 2, nperm = 2, n_cores = 1)
#' rgcca_permutation(A, p_ncomp = c(2,2,3), nperm = 2, n_cores = 1)
#' rgcca_permutation(A, p_ncomp = matrix(c(2,2,3), 3, 3, byrow = T), nperm = 2, n_cores = 1)
#' @export
rgcca_permutation <- function(
    blocks,
    type = "rgcca",
    p_c1 = FALSE,
    p_ncomp = TRUE,
    nperm = 20,
    n_cores = parallel::detectCores() - 1,
    ...) {

    if (any(p_ncomp == FALSE) && any(p_c1 == FALSE))
        stop("Select one parameter among 'p_c1' or 'p_ncomp' to optimize. By default, p_ncomp is selected.")

    ncols <- sapply(blocks, NCOL)
    min_c1s <- sapply(ncols, function(x) 1 / sqrt(x))

    set_c1s <- function(max = 1) {
        if (length(max) == 1)
            f <- quote(max)
        else
            f <- quote(max[x])
        sapply(seq(ncols), function(x) seq(eval(f), min_c1s[x], len = 10))
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

    if (!any(p_c1 == FALSE)) {
        if (identical(p_c1, TRUE))
            p_c1 <- set_c1s()
        else if (class(p_c1) %in% c("data.frame", "matrix"))
            p_c1 <- t(sapply(seq(NROW(p_c1)), function(x) check_tau(p_c1[x, ], blocks, type = "sgcca")))
        else{
            if (any(p_c1 < min_c1s))
                stop(paste0("p_c1 should be upper than 1 / sqrt(NCOL(blocks)) : ", paste0(round(min_c1s, 2), collapse = ",")))
            p_c1 <- check_tau(p_c1, blocks, type = "sgcca")
            p_c1 <- set_c1s(max = p_c1)
        }

        colnames(p_c1) <- names(blocks)
        par <- list("c1", p_c1)
        type <- "sgcca"
    }

    crits <- rgcca_permutation_k(
        blocks,
        par = par,
        perm = FALSE,
        type = type,
        n_cores = 1,
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
    #             "p_c1",
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
        permcrit <- simplify2array(parallel::mclapply(
            seq(nperm),
            function(x){
                res <- rgcca_permutation_k(
                    blocks = blocks,
                    par = par,
                    type = type,
                    n_cores = 1,
                    ...
                )
                return(res)
                },
            mc.cores = n_cores))
    # }

    cat("OK.\n", append = TRUE)

    par <- par[[2]]

    pvals <- sapply(
        seq(NROW(par)), 
        function(i)
            mean(permcrit[i, ] >= crits[i]))
    zs <- sapply(
        seq(NROW(par)), 
        function(i)
            (crits[i] - mean(permcrit[i, ])) / (sd(permcrit[i, ])))

    list(
        pvals = pvals,
        zstat = zs,
        bestpenalties = par[which.max(zs), ],
        permcrit = permcrit,
        crit = crits,
        penalties = par
    )
}
