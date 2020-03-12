#' Regularized (or Sparse) Generalized Canonical Correlation Analysis (R/SGCCA) 
#' 
#' RGCCA is a generalization
#' of regularized canonical correlation analysis to three or more sets of variables. SGCCA extends RGCCA to address the issue of variable selection
#' @details
#' Given \eqn{J} matrices \eqn{\mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}} that represent 
#' \eqn{J} sets of variables observed on the same set of \eqn{n} individuals. The matrices 
#' \eqn{\mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}} must have the same number of rows, 
#' but may (and usually will) have different numbers of columns. The aim of RGCCA is to study 
#' the relationships between these \eqn{J} blocks of variables. It constitutes a general 
#' framework for many multi-block data analysis methods. It combines the power of 
#' multi-block data analysis methods (maximization of well identified criteria) 
#' and the flexibility of PLS path modeling (the researcher decides which blocks 
#' are connected and which are not). Hence, the use of RGCCA requires the construction 
#' (user specified) of a design matrix, (\eqn{\mathbf{C}}), that characterize 
#' the connections between blocks. Elements of the (symmetric) design matrix \eqn{\mathbf{C} = (c_{jk})} 
#' is equal to 1 if block \eqn{j} and block \eqn{k} are connected, and 0 otherwise.
#' The objective is to find a fixed point of the stationary equations related to the RGCCA optimization 
#' problem. The function rgcca() implements a monotonically convergent algorithm (i.e. the bounded
#' criteria to be maximized increases at each step of the iterative procedure) that is very 
#' similar to the PLS algorithm proposed by Herman Wold. Moreover, depending on the 
#' dimensionality of each block \eqn{\mathbf{X}_j}, \eqn{j = 1, \ldots, J}, the primal (when \eqn{n > p_j}) algorithm or 
#' the dual (when \eqn{n < p_j}) algorithm is used (see Tenenhaus et al. 2013). 
#' Moreover, by deflation strategy, rgcca() allow to compute several RGCCA block
#' components (specified by ncomp) for each block. Block components of each block are guaranteed to 
#' be orthogonal with the use of the deflation. The so-called symmetric deflation is considered in
#' this implementation, i.e. each block is deflated with respect to its own component.
#' It should be noted that the numbers of components per block can differ from one block to another. 
#' SGCCA extends RGCCA to address the issue of variable selection. Specifically, 
#' RGCCA is combined with an L1-penalty that gives rise to Sparse GCCA (SGCCA) Blocks are not necessarily fully connected
#' within the SGCCA framework.
#' The SGCCA algorithm is very similar to the RGCCA algorithm and keeps the same monotone 
#' convergence properties (i.e. the bounded criteria to be maximized increases 
#' at each step of the iterative procedure and hits at convergence a stationary point).
#' Moreover, using a deflation strategy, sgcca() enables the computation of several SGCCA block 
#' components (specified by ncomp) for each block. Block components for each block are guaranteed to be orthogonal 
#' when using this deflation strategy. The so-called symmetric deflation is considered in this implementation,
#' i.e. each block is deflated with respect to its own component. 
#' Moreover, we stress that the numbers of components per block could differ from one block to another. 
#' @inheritParams rgccaNa
#' @inheritParams sgccaNa
#' @inheritParams select_analysis

#' @return A RGCCA object
#' @return \item{Y}{A list of \eqn{J} elements. Each element of \eqn{Y} is a matrix that contains the RGCCA components for the corresponding block.}
#' @return \item{a}{A list of \eqn{J} elements. Each element of \eqn{a} is a matrix that contains the outer weight vectors for each block.}
#' @return \item{astar}{A list of \eqn{J} elements. Each element of astar is a matrix defined as Y[[j]][, h] = A[[j]]\%*\%astar[[j]][, h].}
#' @return \item{tau}{A vector or matrix that contains the values of the shrinkage parameters applied to each block and each dimension (user specified).}
#' @return \item{crit}{A vector that contains the values of the criteria across iterations.}
#' @return \item{mode}{A \eqn{1 \times J} vector that contains the formulation ("primal" or "dual") applied to each of the \eqn{J} blocks within the RGCCA alogrithm} 
#' @return \item{AVE}{indicators of model quality based on the Average Variance Explained (AVE): AVE(for one block), AVE(outer model), AVE(inner model).}
#' @return \item{A}{ blocks used in the calculations. Imputed block if imputation method were chosen}
#' @return \item{call}{Call of the function}
#' @references Tenenhaus A. and Tenenhaus M., (2011), Regularized Generalized Canonical Correlation Analysis, Psychometrika, Vol. 76, Nr 2, pp 257-284.
#' @references Tenenhaus A. et al., (2013), Kernel Generalized Canonical Correlation Analysis, submitted.
#' @references Schafer J. and Strimmer K., (2005), A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics. Statist. Appl. Genet. Mol. Biol. 4:32.
#' @examples
#' #############
#' # Example 1 #
#' #############
#' data(Russett)
#' X_agric =as.matrix(Russett[,c("gini","farm","rent")])
#' X_ind = as.matrix(Russett[,c("gnpr","labo")])
#' X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
#' A = list(X_agric, X_ind, X_polit);names(A)=c("Agri","Indus","Polit")
#' #Define the design matrix (output = C) 
#' C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
#' result.rgcca = rgcca(A,type="rgcca", C, tau = c(1, 1, 1),superblock=FALSE,
#'  scheme = "factorial", scale = TRUE)
#' lab = as.vector(apply(Russett[, 9:11], 1, which.max))
#' plot(result.rgcca,type="ind",i_block=1,i_block_y=2,resp=lab)
#' ############################################
#' # Example 2: RGCCA and multiple components #
#' ############################################
#' result.rgcca = rgcca(A,type="rgcca",connection= C, superblock=FALSE,
#' tau = rep(1, 3), ncomp = c(2, 2, 2),
#'                      scheme = "factorial", verbose = TRUE)
#' plot(result.rgcca,resp=lab)
#' plot(result.rgcca,type="ave")
#' plot(result.rgcca,type="network")
#' plot(result.rgcca,type="weight")
#' ############################################
#' # Example : SGCCA #
#' ############################################
#' result.sgcca = rgcca(A,type="sgcca",connection= C, superblock=FALSE,
#' sparsity = rep(0.8, 3), ncomp = c(2, 2, 2),
#'                      scheme = "factorial", verbose = TRUE)
#' plot(result.sgcca,resp=lab)
#' plot(result.sgcca,type="ave")
#' plot(result.sgcca,type="network")
#' plot(result.sgcca,type="weight")

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
#' @seealso \code{\link[RGCCA]{plot.rgcca}}, \code{\link[RGCCA]{print.rgcca}},
#' \code{\link[RGCCA]{rgcca_crossvalidation}},
#' \code{\link[RGCCA]{rgcca_permutation}}
#' \code{\link[RGCCA]{rgcca_predict}} 
rgcca <- function(
    blocks,
    type = "rgcca",
    scale = TRUE,
    sameBlockWeight = TRUE,
    connection = matrix(1,length(blocks),length(blocks)) - diag(length(blocks)),
    scheme = "factorial",
    ncomp = rep(2, length(blocks)),
    tau = rep(1, length(blocks)),
    sparsity = rep(1, length(blocks)),
    init = "svd",
    bias = TRUE,
    tol = 1e-08,
    response = NULL,
    superblock = TRUE,
    method = "complete",
    verbose = FALSE,
    quiet = TRUE,
    knn.k = "all",
    knn.output = "weightedMean",
    knn.klim = NULL,
    knn.sameBlockWeight = TRUE) {

    if (!missing(sparsity) && missing(type))
        type <- "sgcca"

    if (!missing(connection) && missing(superblock))
        superblock <- FALSE

    if (!missing(response) && missing(superblock))
        superblock <- FALSE

    # if (!missing(superblock) && !(missing(response) || missing(connection)))
        

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

    for (i in c("knn.klim")) {
        if (!(i == "knn.klim" && is.null(get(i))))
            check_integer(i, get(i))
    }

    if (!knn.k %in% c("all", "auto"))
        check_integer("knn.k", knn.k)

    for (i in c("superblock","verbose", "scale", "bias", "quiet", "knn.sameBlockWeight"))
        check_boolean(i, get(i))

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

    if (opt$superblock && any(opt$tau) == "optimal")
        stop("Optimal tau is not available with superblock option.")

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

    if (warn_on && !quiet)
        message("Analysis in progress ...")
    
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
            pca.ncp =1,
            prescaling = FALSE,
            quiet=quiet
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
        "type"
    ))
        func_out$call[[i]] <- as.list(environment())[[i]]
   

    class(func_out) <- "rgcca"
    invisible(func_out)
}
