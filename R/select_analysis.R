#' Define the parameters associated with each multi-block component
#' method of the literature.
#'
#' @param method A character string indicating the multi-block component
#' method to consider: rgcca, sgcca, pca, spca, pls, spls, cca,
#' ifa, ra, gcca, maxvar, maxvar-b, maxvar-a, mcoa,cpca-1, cpca-2,
#' cpca-4, hpca, maxbet-b, maxbet, maxdiff-b, maxdiff, maxvar-a,
#' sabscor, ssqcor, ssqcor, ssqcov-1, ssqcov-2, ssqcov, sumcor,
#' sumcov-1, sumcov-2, sumcov, sabscov, sabscov-1, sabscov-2.
#' @inheritParams plot_var_2D
#' @inheritParams set_connection
#' @param blocks List of blocks.
#' @param response Numerical value giving the position of the response block. When
#' the response argument is filled the supervised mode is automatically
#' activated.
#' @param connection Symmetric matrix (J*J) that describes the relationships
#' between blocks. Elements of the connection matrix must be positive ; but
#' usually equal to 1 if block \eqn{j} and block \eqn{k} are connected, and 0
#' otherwise.
#' @param penalty Vector of length J (or character string for 'optimal'
#' setting) indicating the values of the tuning parameters.
#' @param ncomp Vector of length J indicating the number of block components
#' for each block.
#' @param scheme Character string or a function giving the scheme function for
#' covariance maximization among "horst" (the identity function), "factorial"
#'  (the squared values), "centroid" (the absolute values). The scheme function
#'  can be any continously differentiable convex function and it is possible to
#'  design explicitely the sheme function (e.g. function(x) x^4) as argument of
#'  rgcca function.  See (Tenenhaus et al, 2017) for details.
#' @param verbose Logical value indicating whether the warnings are displayed.
#' @param quiet Logical value indicating if warning messages are reported.
#' @return \item{blocks}{List of blocks.}
#' @return \item{scheme}{Character string or a function giving the scheme
#' function used for covariance maximization.}
#' @return \item{penalty}{Vector of length J (or character string for
#' 'optimal' setting) indicating the values of the tuning parameters.}
#' @return \item{ncomp}{Vector of length J indicating the number of block
#' components or each block.}
#' @return \item{connection}{Symmetric matrix (J*J) that describes the
#' relationships between blocks.}
#' @return \item{superblock}{Logical value indicating if superblock is
#' included in the analysis.}

select_analysis <- function(
    blocks,
    connection = 1 - diag(length(blocks)),
    penalty = rep(1, length(blocks)),
    group_sparsity = NULL,
    ncomp = rep(1, length(blocks)),
    scheme = "centroid",
    superblock = TRUE,
    method  = "rgcca",
    verbose = TRUE,
    quiet = FALSE,
    response = NULL){

    J <- length(blocks)
    msg_superblock <- "A superblock is considered."
    msg_type <- paste0("By using a ", toupper(method), ", ")
    warn.method.value <- warn.method.par <- warn.msg.super <- character(0)

    if (quiet)
        verbose <- FALSE

    ### SETTINGS ###

    warnParam <- function(param, x) {
        warn.method.par <<- c(warn.method.par, paste(deparse(substitute(param))))
        warn.method.value <<- c(warn.method.value, toString(x))
    }

    setPenalty <- function(x) {
        warnParam(penalty, x)
        return(x)
    }

    setScheme <- function(x) {
        warnParam(scheme, x)
        return(x)
    }

    set_connection <- function(x) {
        warnParam(connection, paste(deparse(substitute(x))))
        return(x)
    }

    warnSuper <- function(x) {
        if (class(x) %in% c("matrix", "data.frame") &&
            NCOL(x) < (length(blocks)) &&
            is.null(response)){
            warn.msg.super <<- c(warn.msg.super, deparse(substitute(x)))
            return(cbind(x, 1))
        }else if (length(x) < (length(blocks)) && is.null(response)) {
            warn.msg.super <<- c(warn.msg.super, deparse(substitute(x)))
            if(deparse(substitute(x)) == "ncomp")
                return(c(x, max(x)))
            else
                return(c(x, 1))
        } else{
            return(x)
        }
    }

    setSuperblock <- function(verbose = TRUE) {
        blocks <<- c(blocks, superblock = list(Reduce(cbind, blocks)))
        superblock <<- TRUE
        connection <<- NULL
        ncomp <<- warnSuper(ncomp)
    }

    set2Block <- function(method) {
        check_nblocks(blocks, method)

        scheme <<- setScheme("horst")
        connection <<- set_connection(1 - diag(2))
    }

    ### CHECK TYPES ###

    if (length(grep("[sr]gcca", tolower(method))) == 1) {
        if (superblock) {
            setSuperblock(FALSE)
            penalty <- warnSuper(penalty)
        } else
            superblock <- FALSE
    } else
        superblock <- FALSE

    if (length(grep("^s?pca$", tolower(method))) == 1) {
        check_nblocks(blocks, method)

        scheme <- setScheme("horst")
        setSuperblock()
        if (tolower(method) == "pca")
            penalty <- setPenalty(c(1, 1))
    }

    # 2 Blocks cases
    else if (tolower(method) %in% c("cca", "ra", "ifa", "pls", "spls")) {
        set2Block(method)

        if (tolower(method) == "cca")
            penalty <- setPenalty(c(0, 0))

        else if (tolower(method) %in% c("ifa", "pls"))
            penalty <- setPenalty(c(1, 1))

        else if (tolower(method) == "ra")
            penalty <- setPenalty(c(1, 0))

    }

    # Design matrix of 1 values everywhere
    else if (tolower(method) %in% c("sumcor",
                                "ssqcor",
                                "sabscor",
                                "sumcov-1",
                                "maxbet",
                                "ssqcov-1",
                                "maxbet-b",
                                "sabscov-1")) {
        connection <- set_connection(matrix(1, J, J))

        # COR models
        if (tolower(method) %in% c("sumcor", "ssqcor", "sabscor")) {
            penalty <- setPenalty(rep(0, J))

            switch(
                tolower(method),
                "sumcor" = {
                    scheme <- setScheme("horst")
                },
                "ssqcor" = {
                    scheme <- setScheme("factorial")
                },
                "sabscor" = {
                    scheme <- setScheme("centroid")
                }
            )
        }

        # COV models
        else if (tolower(method) %in% c(
            "sumcov-1",
            "maxbet",
            "ssqcov-1",
            "maxbet-b",
            "sabscov-1"
        )) {
            penalty <- setPenalty(rep(1, J))

            if (tolower(method) %in% c("sumcov-1", "maxbet"))
                scheme <- setScheme("horst")

            else if (tolower(method) %in% c("ssqcov-1", "maxbet-b"))
                scheme <- setScheme("factorial")

            else if (tolower(method) %in% c("sabscov-1"))
                scheme <- setScheme("centroid")

        }

        # Design matrix with 1 values everywhere except
        # on the diagonal equals to 0
    }

    else if (tolower(method) %in% c("sumcov",
                                  "sumcov-2",
                                  "maxdiff",
                                  "ssqcov",
                                  "ssqcov-2",
                                  "maxdiff-b",
                                  "sabscov-2")) {
        connection <- set_connection(1 - diag(J))

        if (tolower(method) %in% c("sumcov", "sumcov-2", "maxdiff")) {
            scheme <- setScheme("horst")
            penalty <- setPenalty(rep(1, J))
        }

        else if (tolower(method) %in% c("ssqcov", "ssqcov-2", "maxdiff-b")) {
            scheme <- setScheme("factorial")
            penalty <- setPenalty(rep(1, J))
        }

    }

    # Models with a superblock
    else if (tolower(method) %in% c("gcca", "maxvar", "maxvar-b",
                                  "cpca-1", "cpca-2", "maxvar-a", "mcoa",
                                  "cpca-4", "hpca")){
        setSuperblock()

        if (tolower(method) %in% c("gcca", "maxvar", "maxvar-b")) {
            scheme <- setScheme("factorial")
            penalty <- setPenalty(rep(0, J + 1))
        }

        else if (tolower(method) == "cpca-1") {
            scheme <- function(x) x
            penalty <- setPenalty(c(rep(1, J), 0))
        }

        else if (tolower(method) %in% c("maxvar-a", "cpca-2", "mcoa")){
            scheme <- setScheme("factorial")
            penalty <- setPenalty(c(rep(1, J), 0))
        }

        else if (tolower(method) %in% c("hpca", "cpca-4")) {
            scheme <- function(x) x^4
            penalty <- setPenalty(c(rep(1, J), 0))
        }
        }

    ### WARNINGS ###
    n = length(warn.method.par)
    if (verbose & n > 0) {
        setPlural = function(x = warn.method.par,
                             y = warn.method.value,
                             sep = " and "){
            warn.method.par <<- paste0(x, collapse = sep)
            warn.method.value <<- paste0(y, collapse = sep)
        }

        if (n > 1) {
            grammar = "s were respectively"
            if (n == 2)
                setPlural()
            else{
                warn.method = c(warn.method.par[n], warn.method.value[n])
                setPlural(warn.method.par[-n], warn.method.value[-n], ", ")
                setPlural(c(warn.method.par, warn.method[1]),
                        c(warn.method.value, warn.method[2]))
            }
        } else
            grammar <- " was"

        msg <- paste0(warn.method.par, " parameter",
                      grammar, " set to ", warn.method.value)

        if (superblock & tolower(method) != "pca")
            msg <- paste0(msg, " and ", msg_superblock)

        warning(paste0(msg_type, msg , "."))
    }

    if (verbose & superblock) {
        if (n < 0) paste0(msg_superblock, msg_superblock)
    }

    if (!quiet & length(warn.msg.super) > 0) {
        if (length(warn.msg.super) > 1) {
            warn.msg.super <- paste(warn.msg.super, collapse = " and ")
            grammar <- "were those"
        } else
            grammar <- "was the one"
    }
    
    ##TODO : add some checks for group_sparsity 

    return(list(scheme = scheme,
                penalty = penalty,
                group_sparsity = group_sparsity,
                ncomp = ncomp,
                connection = connection,
                superblock = superblock)
           )
}
