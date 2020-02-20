#' Performs a r/sgcca
#'
#' Performs a r/sgcca with predefined parameters
#' @inheritParams select_analysis
#' @inheritParams rgccaNa
#' @inheritParams sgccaNa
#' @param scale A boolean scaling the blocks
#' @param init A character among "svd" (Singular Value Decompostion) or "random"
#' for algorithm initialization
#' @param bias A boolean for a biased variance estimator
#' @param verbose A boolean to display the progress of the analysis
#' @param response An integer giving the index of a block considered as a 
#' response among a list of blocks
#' @param tol An integer for the stopping value for convergence
#' @param sameBlockWeight If TRUE, all blocks are weighted by their own variance: all the blocks have the same weight
#' @return A RGCCA object
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca(blocks)
#' @export
#' @import ggplot2
#' @importFrom grDevices dev.off rgb colorRamp pdf colorRampPalette
#' @importFrom graphics plot
#' @importFrom stats cor quantile runif sd na.omit p.adjust pnorm qnorm weights
#' @importFrom utils read.table write.table packageVersion installed.packages head
#' @importFrom scales hue_pal
# @importFrom optparse OptionParser make_option parse_args
# @importFrom plotly layout ggplotly style plotly_build %>% plot_ly add_trace
# @importFrom visNetwork visNetwork visNodes visEdges
# @importFrom igraph graph_from_data_frame V<- E<-
#' @importFrom methods is
rgcca <- function(
    blocks,
    connection = 1 - diag(length(blocks)),
    response = NULL,
    superblock = TRUE,
    tau = rep(1, length(blocks)),
    sparsity = rep(1, length(blocks)),
    ncomp = rep(2, length(blocks)),
    type = "rgcca",
    verbose = FALSE,
    scheme = "factorial",
    scale = TRUE,
    init = "svd",
    bias = TRUE,
    tol = 1e-08,
    quiet = TRUE,
    sameBlockWeight = TRUE,
    method = "complete",
    knn.k = "all",
    knn.output = "weightedMean",
    knn.klim = NULL,
    knn.sameBlockWeight = TRUE,
    pca.ncp = 1) {

    if (tolower(type) %in% c("sgcca", "spca", "spls")) {
        if (!missing(tau) && missing(sparsity))
           stop(paste0("sparsity parameter required for ", tolower(type), "(instead of tau)."))
        gcca <- sgccaNa
        par <- "sparsity"
        penalty <- sparsity
    }else{
        if (!missing(sparsity) & missing(tau))
           stop(paste0("tau parameter required for ", tolower(type), "(instead of sparsity)."))
        gcca <- rgccaNa
        par <- "tau"
        penalty <- tau
    }
    
    match.arg(init, c("svd", "random"))
    match.arg(knn.output, c("mean", "random", "weightedMean" ))
    check_method(type)
    if (!is.null(response))
        check_blockx("response", response, blocks)
    check_integer("tol", tol, float = TRUE, min = 0)

    for (i in c("pca.ncp", "knn.klim")) {
        if (!(i == "knn.klim" && is.null(get(i))))
            check_integer(i, get(i))
    }

    if (!knn.k %in% c("all", "auto"))
        check_integer("knn.k", knn.k)

    for (i in c("superblock","verbose", "scale", "bias", "quiet", "knn.sameBlockWeight"))
        check_boolean(i, get(i))

    choices <- list(c("horst", "factorial", "centroid"))
    choice <- c(scheme)
    for (i in length(choices)) {
        if (!choice[i] %in% (choices[[i]]) && !is.function(choice[i]))
            stop(paste0(choice[i], " must be one of '", paste(choices[[i]], collapse = ", "), "' or a function."))
    }

    penalty <- elongate_arg(penalty, blocks)
    ncomp <- elongate_arg(ncomp, blocks)

    opt <- select_analysis(
        blocks = blocks,
        connection = connection,
        penalty = penalty,
        ncomp = ncomp,
        scheme = scheme,
        superblock = superblock,
        type  = type,
        quiet = quiet,
        response = response
    )

    opt$blocks <- scaling(blocks, scale,sameBlockWeight = sameBlockWeight)
    opt$superblock <- check_superblock(response, opt$superblock, !quiet)
    opt$blocks <- set_superblock(opt$blocks, opt$superblock, type, !quiet)

    if(opt$superblock && any(opt$tau)=="optimal")
    {
        stop("Optimal tau is not available with superblock option.")
    }

    if (!is.null(response)) {
        # || tolower(type) == "ra"
        response <- check_blockx("response", response, opt$blocks)
        pars <- c("blocks", "ncomp", "penalty")
        for (i in seq(length(pars)))
            opt[[pars[i]]] <- c(opt[[pars[i]]][-response], opt[[pars[i]]][response])
    }


    if (!is.matrix(opt$connection) || !is.null(response))
        opt$connection <- set_connection(
            opt$blocks,
            (opt$superblock | !is.null(response))
        )

    check_connection(opt$connection, opt$blocks)
    opt$penalty <- check_tau(opt$penalty, opt$blocks, type)
    opt$ncomp <- check_ncomp(opt$ncomp, opt$blocks)

    warn_on <- FALSE

    if (any(sapply(opt$blocks, NCOL) > 1000)) {
            # if( (type <-<- "sgcca" && tau > 0.3) || type !<- "sgcca" )
            warn_on <- TRUE
    }

    if (warn_on || !quiet)
        message("RGCCA in progress ...")

    func <- quote(
        gcca(
            blocks = opt$blocks,
            connection = opt$connection,
            ncomp = opt$ncomp,
            verbose = verbose,
            scheme = opt$scheme,
            scale = scale,
            init = init,
            bias = bias,
            tol = tol,
            sameBlockWeight = sameBlockWeight,
            method = method,
            knn.k = knn.k,
            knn.output = knn.output,
            knn.klim = knn.klim,
            knn.sameBlockWeight = knn.sameBlockWeight,
            pca.ncp = pca.ncp,
            prescaling = FALSE
        )
    )

    func[[par]] <- opt$penalty
    func_out <- eval(as.call(func))$rgcca

    for (i in c("a", "astar", "Y")) {
        names(func_out[[i]]) <- names(opt$blocks)
        for (j in seq(length(opt$blocks))) {
            if (i %in%  c("a", "astar") && NCOL(opt$blocks[[j]]) == 1)
                row.names(func_out[[i]][[j]]) <- colnames(opt$blocks[[j]])
        }
    }
    names(func_out$AVE$AVE_X) <- names(opt$blocks)

    class(func_out) <- tolower(type)
    func_out$call <- list(
        blocks = opt$blocks,
        connection = opt$connection,
        superblock = opt$superblock,
        ncomp = opt$ncomp,
        scheme = opt$scheme
    )

    func_out$call[[par]] <- opt$penalty

    for (i in c(
        "scale",
        "init",
        "bias",
        "tol",
        "verbose",
        "response",
        "sameBlockWeight",
        "method",
        "knn.k",
        "knn.output",
        "knn.klim",
        "pca.ncp",
        "type"
    ))
        func_out$call[[i]] <- as.list(environment())[[i]]

    # adding potential modified A to the list of outputs 
    # (if imputed or restricted -only complete)
    if (method != "nipals")
        func_out$usedBlocks <- func_out$A

    class(func_out) <- "rgcca"
    invisible(func_out)
}
