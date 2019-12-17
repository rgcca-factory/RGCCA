#' Performs a r/sgcca
#'
#' Performs a r/sgcca with predefined parameters
#' @inheritParams select_analysis
#' @param scale A boolean scaling the blocks
#' @param init A character among "svd" (Singular Value Decompostion) or "random"
#' for alorithm initialization
#' @param bias A boolean for a biased variance estimator
#' @param verbose A boolean to display the progress of the analysis
#' @param response An integer giving the index of a block considered as a 
#' response among a list of blocks
#' @param tol An integer for the stopping value for convergence
#' @return A RGCCA object
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca.analyze(blocks)
#' @export
#' @import RGCCA
#' @import ggplot2
#' @importFrom grDevices dev.off rgb colorRamp pdf colorRampPalette
#' @importFrom graphics plot
#' @importFrom stats cor quantile runif sd na.omit p.adjust pnorm qnorm weights
#' @importFrom utils read.table write.table packageVersion installed.packages head
#' @importFrom scales hue_pal
#' @importFrom optparse OptionParser make_option parse_args
#' @importFrom plotly layout ggplotly style plotly_build %>% plot_ly add_trace
#' @importFrom visNetwork visNetwork visNodes visEdges
#' @importFrom igraph graph_from_data_frame V<- E<-
#' @importFrom methods is
rgcca.analyze <- function(
    blocks,
    connection = 1 - diag(length(blocks)),
    response = NULL,
    superblock = TRUE,
    tau = rep(1, length(blocks)),
    ncomp = rep(2, length(blocks)),
    type = "rgcca",
    verbose = TRUE,
    scheme = "factorial",
    scale = TRUE,
    init = "svd",
    bias = TRUE, 
    tol = 1e-08) {

    tau <- elongate_arg(tau, blocks)
    ncomp <- elongate_arg(ncomp, blocks)

    opt <- select_analysis(
        blocks = blocks,
        connection = connection,
        tau = tau,
        ncomp = ncomp,
        scheme = scheme,
        superblock = superblock,
        type  = type
    )

    opt$blocks <- scaling(blocks, scale)
    superblock <- check_superblock(response, opt$superblock)
    opt$blocks <- set_superblock(opt$blocks, opt$superblock, type)

    if (!is.null(response)) {
        # || tolower(type) == "ra"
        response <- check_blockx("response", response, opt$blocks)
        par <- c("blocks", "ncomp", "tau")
        for (i in seq(length(par)))
            opt[[par[i]]] <- c(opt[[par[i]]][-response], opt[[par[i]]][response])
    }

    if (!is.matrix(opt$connection) || !is.null(response))
        opt$connection <- set_connection(
            opt$blocks,
            (opt$superblock | !is.null(response))
        )

    check_connection(opt$connection, opt$blocks)
    opt$tau <- check_tau(opt$tau, opt$blocks, type)
    opt$ncomp <- check_ncomp(opt$ncomp, opt$blocks)

    warn_on <- FALSE

    if (any(sapply(opt$blocks, NCOL) > 1000)) {
            # if( (type <-<- "sgcca" && tau > 0.3) || type !<- "sgcca" )
            warn_on <- TRUE
    }

    if (warn_on & verbose)
        message("RGCCA in progress ...")

    if (tolower(type) == "sgcca") {
        gcca <- RGCCA::sgcca
        par <- "c1"
    } else{
        gcca <- RGCCA::rgcca
        par <- "tau"
    }

    func <- quote(
        gcca(
            A = opt$blocks,
            C = opt$connection,
            ncomp = opt$ncomp,
            verbose = FALSE,
            scheme = opt$scheme,
            scale = FALSE,
            init = init,
            bias = bias,
            tol = tol
        )
    )
    func[[par]] <- opt$tau

    func_out <- eval(as.call(func))
    for (i in c("a", "astar", "Y"))
        names(func_out[[i]]) <- names(opt$blocks)
    names(func_out$AVE$AVE_X) <- names(opt$blocks)
    func_out$blocks <- opt$blocks
    func_out$superblock <- opt$superblock
    for (i in c("scale", "init", "bias", "tol", "verbose"))
        func_out[[i]] <- as.list(environment())[[i]]

    class(func_out) <- tolower(type)
    invisible(func_out)
}
