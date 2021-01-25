#' Define the parameters associated with each multi-block component
#' method of the literature.
#'
#' @param type A character string indicating the multi-block component
#' method to consider: rgcca, sgcca, pca, spca, pls, spls, cca, 
#' ifa, ra, gcca, maxvar, maxvar-b, maxvar-a, mcoa,cpca-1, cpca-2, 
#' cpca-4, hpca, maxbet-b, maxbet, maxdiff-b, maxdiff, maxvar-a,  
#' sabscor, ssqcor, ssqcor, ssqcov-1, ssqcov-2, ssqcov, sumcor, 
#' sumcov-1, sumcov-2, sumcov, sabscov, sabscov-1, sabscov-2.
#' @inheritParams plot_var_2D
#' @inheritParams set_connection
#' @param blocks A list of blocks
#' @param response An integer giving the position of the response block. When 
#' the response argument is filled the supervised mode is automatically 
#' activated.  
#' @param connection A symmetric matrix (J*J) that describes the relationships 
#' between blocks. Elements of the connection matrix must be positive ; but 
#' usually equal to 1 if block \eqn{j} and block \eqn{k} are connected, and 0 
#' otherwise.
#' @param penalty A vector of length J (or character string for 'optimal' 
#' setting) indicating the values of the tuning parameters.
#' @param ncomp A vector of length J indicating the number of block components 
#' for each block.
#' @param scheme A character string or a function giving the scheme function for 
#' covariance maximization among "horst" (the identity function), "factorial"
#'  (the squared values), "centroid" (the absolute values). The scheme function 
#'  can be any continously differentiable convex functin and it is possible to 
#'  design explicitely the sheme function (e.g. function(x) x^4) as argument of 
#'  rgcca function.  See (Tenenhaus et al, 2017) for details.
#' @param verbose A logical value indicating whether the warnings are displayed
#' @param quiet A boolean hidding the warnings
#' @return \item{blocks}{A list of blocks.}
#' @return \item{scheme}{A character string or a function giving the scheme 
#' function for covariance maximization.}
#' @return \item{penalty}{A vector of length J (or character string for 
#' 'optimal' setting) indicating the values of the tuning parameters.}
#' @return \item{ncomp}{A vector of length J indicating the number of block 
#' components or each block.}
#' @return \item{connection}{A symmetric matrix (J*J) that describes the 
#' relationships between blocks}
#' @return \item{superblock}{A logical value indicating if superblock is 
#' included in the analysis}

select_analysis <- function(
    blocks,
    connection = 1 - diag(length(blocks)),
    penalty = rep(1, length(blocks)),
    ncomp = rep(1, length(blocks)),
    scheme = "centroid",
    superblock = TRUE,
    type  = "rgcca",
    verbose = TRUE,
    quiet = FALSE,
    response = NULL){

    J <- length(blocks)
    msg_superblock <- "a superblock is used"
    msg_type <- paste0("By using a ", toupper(type), ", ")
    warn.type.value <- warn.type.par <- warn.msg.super <- character(0)

    if (quiet)
        verbose <- FALSE

    ### SETTINGS ###

    warnParam <- function(param, x) {
        warn.type.par <<- c(warn.type.par, paste(deparse(substitute(param))))
        warn.type.value <<- c(warn.type.value, toString(x))
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

    set2Block <- function(type) {
        if (length(blocks) != 2)
            check_nblocks(blocks, type)

        scheme <<- setScheme("horst")
        connection <<- set_connection(1 - diag(2))
    }

    ### CHECK TYPES ###

    if (length(grep("[sr]gcca", tolower(type))) == 1) {
        if (superblock) {
            setSuperblock(FALSE)
            penalty <- warnSuper(penalty)
        } else
            superblock <- FALSE
    } else
        superblock <- FALSE

    if (length(grep("^s?pca$", tolower(type))) == 1) {
        if (length(blocks) != 1)
            check_nblocks(blocks, type)

        scheme <- setScheme("horst")
        setSuperblock()
        if (tolower(type) == "pca")
            penalty <- setPenalty(c(1, 1))
    }

    # 2 Blocks cases
    else if (tolower(type) %in% c("cca", "ra", "ifa", "pls", "spls")) {
        set2Block(type)

        if (tolower(type) == "cca")
            penalty <- setPenalty(c(0, 0))

        else if (tolower(type) %in% c("ifa", "pls"))
            penalty <- setPenalty(c(1, 1))

        else if (tolower(type) == "ra")
            penalty <- setPenalty(c(1, 0))

    }

    # Design matrix of 1 values everywhere
    else if (tolower(type) %in% c("sumcor",
                                "ssqcor",
                                "sabscor",
                                "sumcov-1",
                                "maxbet",
                                "ssqcov-1",
                                "maxbet-b",
                                "sabscov-1")) {
        connection <- set_connection(matrix(1, J, J))

        # COR models
        if (tolower(type) %in% c("sumcor", "ssqcor", "sabscor")) {
            penalty <- setPenalty(rep(0, J))

            switch(
                tolower(type),
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
        else if (tolower(type) %in% c(
            "sumcov-1",
            "maxbet",
            "ssqcov-1",
            "maxbet-b",
            "sabscov-1"
        )) {
            penalty <- setPenalty(rep(1, J))

            if (tolower(type) %in% c("sumcov-1", "maxbet"))
                scheme <- setScheme("horst")

            else if (tolower(type) %in% c("ssqcov-1", "maxbet-b"))
                scheme <- setScheme("factorial")

            else if (tolower(type) %in% c("sabscov-1"))
                scheme <- setScheme("centroid")

        }

        # Design matrix with 1 values everywhere except 
        # on the diagonal equals to 0
    }

    else if (tolower(type) %in% c("sumcov",
                                  "sumcov-2",
                                  "maxdiff", 
                                  "ssqcov", 
                                  "ssqcov-2", 
                                  "maxdiff-b",
                                  "sabscov-2")) {
        connection <- set_connection(1 - diag(J))

        if (tolower(type) %in% c("sumcov", "sumcov-2", "maxdiff")) {
            scheme <- setScheme("horst")
            penalty <- setPenalty(rep(1, J))
        }

        else if (tolower(type) %in% c("ssqcov", "ssqcov-2", "maxdiff-b")) {
            scheme <- setScheme("factorial")
            penalty <- setPenalty(rep(1, J))
        }

    }

    # Models with a superblock
    else if (tolower(type) %in% c("gcca", "maxvar", "maxvar-b", 
                                  "cpca-1", "cpca-2", "maxvar-a", "mcoa", 
                                  "cpca-4", "hpca")){
        setSuperblock()

        if (tolower(type) %in% c("gcca", "maxvar", "maxvar-b")) {
            scheme <- setScheme("factorial")
            penalty <- setPenalty(rep(0, J + 1))
        }
        
        else if (tolower(type) == "cpca-1") {
            scheme <- function(x) x
            penalty <- setPenalty(c(rep(1, J), 0))
        }

        else if (tolower(type) %in% c("maxvar-a", "cpca-2", "mcoa")){
            scheme <- setScheme("factorial")
            penalty <- setPenalty(c(rep(1, J), 0))
        }
        
        else if (tolower(type) %in% c("hpca", "cpca-4")) {
            scheme <- function(x) x^4
            penalty <- setPenalty(c(rep(1, J), 0))
        }
        }

    ### WARNINGS ###
    n = length(warn.type.par)
    if (verbose & n > 0) {
        setPlural = function(x = warn.type.par,
                             y = warn.type.value,
                             sep = " and "){
            warn.type.par <<- paste0(x, collapse = sep)
            warn.type.value <<- paste0(y, collapse = sep)
        }

        if (n > 1) {
            grammar = "s were respectively"
            if (n == 2)
                setPlural()
            else{
                warn.type = c(warn.type.par[n], warn.type.value[n])
                setPlural(warn.type.par[-n], warn.type.value[-n], ", ")
                setPlural(c(warn.type.par, warn.type[1]),
                        c(warn.type.value, warn.type[2]))
            }
        } else
            grammar <- " was"

        msg <- paste0(warn.type.par, " parameter",
                      grammar, " set to ", warn.type.value)

        if (superblock & tolower(type) != "pca")
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

    return(list(scheme = scheme,
                penalty = penalty,
                ncomp = ncomp,
                connection = connection,
                superblock = superblock)
           )
}
