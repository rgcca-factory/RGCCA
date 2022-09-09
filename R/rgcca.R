#' Regularized (or Sparse, or Multiway) Generalized Canonical Correlation
#' Analysis (S/M/RGCCA)
#'
#' RGCCA is a generalization of regularized canonical correlation analysis to
#' three or more sets of variables. SGCCA extends RGCCA to address the issue of
#' variable selection. MGCCA extends RGCCA to address the issue of tensor
#' structured data.
#'
#' @details
#' Given J matrices \eqn{\mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}} that
#' represent \eqn{J} sets of variables observed on the same set of \eqn{n}
#' individuals. The matrices \eqn{\mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}}
#' must have the same number of rows, but may (and usually will) have different
#' numbers of columns. The aim of RGCCA is to study the relationships between
#' these \eqn{J} blocks of variables. It constitutes a general framework for
#' many multi-block data analysis methods (see Tenenhaus and Tenenhaus, 2011 ;
#' Tenenhaus et al. 2017). It combines the power of multi-block data analysis
#' methods (maximization of well identified criteria) and the flexibility of
#' PLS path modeling (the researcher decides which blocks are connected and
#' which are not). Hence, the use of RGCCA requires the construction (user
#' specified) of a design matrix, (\eqn{\mathbf{C}}), that characterizes the
#' connections between blocks. Elements of the (symmetric) design matrix
#' \eqn{\mathbf{C} = (c_{jk})} are positive (and usually equal to 1 if block
#' \eqn{j} and block \eqn{k} are connected, and 0 otherwise). The function
#' rgcca() implements a monotone global convergent algorithm - i.e. the
#' bounded criteria to be maximized increases at each step of the iterative
#' procedure and hits, at convergence a stationary point of the RGCCA
#' optimization problem. Moreover, depending on the dimensionality of each
#' block \eqn{\mathbf{X}_j}, \eqn{j = 1, \ldots, J}, the primal (when
#' \eqn{n > p_j}) algorithm or the dual (when \eqn{n < p_j}) algorithm is used
#' (see Tenenhaus et al. 2015). At last, a deflation strategy is used to compute
#' several RGCCA block components (specified by ncomp) for each block. Block
#' components of each block are guaranteed to be orthogonal. The so-called
#' symmetric deflation is implemented (i.e. each block is deflated with respect
#' to its own component). It should be noted that the numbers of components
#' per block can differ from one block to another.
#' SGCCA extends RGCCA to address the issue of variable selection
#' (Tenenhaus et al, 2014).
#' Specifically, RGCCA is combined with an L1-penalty that gives rise to Sparse
#' GCCA (SGCCA). The SGCCA algorithm is very similar to the RGCCA algorithm and
#' keeps the same convergence properties (i.e. the bounded criteria to be
#' maximized increases at each step of the iterative procedure and hits at
#' convergence a stationary point). Moreover, using a deflation strategy,
#' sgcca() enables the computation of several SGCCA orthogonal block components
#' (specified by ncomp) for each block.
#' MGCCA extends RGCCA to address the issue of tensor structured data.
#' Specifically, RGCCA is combined with a Kronecker constraint that gives rise
#' to Multiway GCCA (MGCCA). The MGCCA algorithm is very similar to the RGCCA
#' algorithm and keeps the same convergence properties (i.e. the bounded
#' criteria to be maximized increases at each step of the iterative procedure
#' and hits at convergence a stationary point). Moreover, using a deflation
#' strategy, mgcca() enables the computation of several MGCCA orthogonal
#' block components (specified by ncomp) for each block. MGCCA can handle
#' blocks of variables structured as higher order arrays.
#' The rgcca() function can handle missing values using a NIPALS type algorithm
#' (non-linear iterative partial least squares algorithm) described in
#' (Tenenhaus et al, 2005). Guidelines describing how to use RGCCA in practice
#' are provided in (Garali et al., 2017).
#' @inheritParams rgccaNa
#' @inheritParams sgccaNa
#' @inheritParams mgccaNa
#' @inheritParams select_analysis
#' @return A RGCCA object
#' @return \item{Y}{A list of \eqn{J} elements. Each element of the list \eqn{Y}
#' is a matrix that contains the RGCCA block components for the corresponding
#' block.}
#' @return \item{a}{A list of \eqn{J} elements. Each element of the list \eqn{a}
#' is a matrix of block weight vectors for the corresponding block.}
#' @return \item{astar}{A list of \eqn{J} elements. Each element of astar is a
#' matrix defined as Y[[j]][, h] = A[[j]]\%*\%astar[[j]][, h].}
#' @return \item{tau}{Either a vector of length J or a matrix of dimension
#' \eqn{\mathrm{max}(ncomp) \times J} containing the values of the shrinkage
#' parameters. tau varies from 0 (maximizing the correlation) to 1 (maximizing
#' the covariance). If tau = "optimal" the shrinkage paramaters are estimated
#' for each block and each dimension using the Schafer and Strimmer (2005)
#' analytical formula. If tau is a vector of length J, tau[j] is identical
#' across the dimensions of block \eqn{\mathbf{X}_j}. If tau is a matrix,
#' tau[k, j] is associated with \eqn{\mathbf{X}_{jk}} (\eqn{k}th residual
#' matrix of block \eqn{j}). tau can be also estimated using
#' \link{rgcca_permutation}.}
#' @return \item{crit}{A list of vector of length max(ncomp). Each vector of
#' the list is related to one specific deflation stage and reports the values
#' of the criterion for this stage across iterations.}
#' @return \item{primal_dual}{A \eqn{1 \times J} vector that contains the
#' formulation ("primal" or "dual") applied to each of the \eqn{J} blocks
#' within the RGCCA alogrithm.}
#' @return \item{AVE}{A list of numerical values giving the indicators of model
#' quality based on the Average Variance Explained (AVE): AVE(for each block),
#' AVE(outer model), AVE(inner model).}
#' @return \item{A}{A list that contains the J blocks of variables X1, X2, ...,
#' XJ. Block Xj is a matrix of dimension n x p_j where p_j is the number of
#' variables in X_j. These blocks are imputed when an imputation strategy is
#' selected.}
#' @return \item{call}{Call of the function}
#' @references Garali I, Adanyeguh IM, Ichou F, Perlbarg V, Seyer A, Colsch B,
#' Moszer I, Guillemot V, Durr A, Mochel F, Tenenhaus A. A strategy for
#' multimodal data integration: application to biomarkers identification
#' in spinocerebellar ataxia. Briefings in Bioinformatics. 2018 Nov 27;19(6):1356-1369.
#' @references Tenenhaus M., Tenenhaus A. and Groenen P. J. (2017). Regularized
#' generalized canonical correlation analysis: a framework for sequential
#' multiblock component methods. Psychometrika, 82(3), 737-777.
#' @references Tenenhaus A., Philippe C. and Frouin, V. (2015). Kernel
#' generalized canonical correlation analysis. Computational Statistics and
#' Data Analysis, 90, 114-131.
#' @references Tenenhaus A., Philippe C., Guillemot V., Le Cao K. A., Grill J.
#' and Frouin, V., Variable selection for generalized canonical correlation
#' analysis, Biostatistics, vol. 15, no. 3, pp. 569-583, 2014.
#' @references Tenenhaus A. and Tenenhaus M., (2011). Regularized Generalized
#' Canonical Correlation Analysis, Psychometrika, Vol. 76, Nr 2, pp 257-284.
#' @references Schafer J. and Strimmer K. (2005). A shrinkage approach to
#' large-scale covariance matrix estimation and implications for functional
#' genomics. Statistical Applications in Genetics and Molecular Biology 4:32.
#' @examples
#' ####################
#' # Example 1: RGCCA #
#' ####################
#' # Create the dataset
#' data(Russett)
#' blocks = list(agriculture = Russett[, seq(3)],
#'     industry = Russett[, 4:5],
#'     politic = Russett[, 6:11])
#'
#' # Blocks are fully connected, factorial scheme and tau =1 for all blocks is
#' # used by default
#' fit.rgcca = rgcca(blocks=blocks, type = "rgcca", connection = 1-diag(3),
#'                   scheme = "factorial", tau = rep(1, 3))
#' print(fit.rgcca)
#' plot(fit.rgcca, type = "weight", block = 3)
#' politic = as.vector(apply(Russett[, 9:11], 1, which.max))
#' plot(fit.rgcca, type = "ind", block = 1:2, comp = rep(1, 2), resp = politic)
#'
#' ############################################
#' # Example 2: RGCCA and multiple components #
#' ############################################
#' fit.rgcca = rgcca(blocks, type = "rgcca",
#'                   connection = C, superblock = FALSE,
#'                   tau = rep(1, 3), ncomp = c(2, 2, 2),
#'                   scheme = "factorial", verbose = TRUE)
#'
#' politic = as.vector(apply(Russett[, 9:11], 1, which.max))
#' plot(fit.rgcca, type = "ind", block = 1:2,
#'      comp = rep(1, 2), resp = politic)
#'
#' plot(fit.rgcca, type = "ave")
#' plot(fit.rgcca, type = "network")
#' plot(fit.rgcca, type = "weight", block = 1)
#' plot(fit.rgcca, type = "cor")
#'
#' ##################################
#' # Example 3: Sparse GCCA (SGCCA) #
#' ##################################
#'
#' # Tune the model to find the best sparsity coefficients (all the blocks are
#' # connected together)
#' perm.out = rgcca_permutation(blocks, n_cores = 1,
#'                              par_type = "sparsity", n_perms = 10)
#' print(perm.out)
#' plot(perm.out)
#'
#' fit.sgcca = rgcca(blocks, sparsity = perm.out$bestpenalties)
#' plot(fit.sgcca, type = "network")
#' plot(fit.sgcca, type = "ave")
#'
#' # Select the most significant variables
#' b = bootstrap(fit.sgcca, n_cores = 1, n_boot = 100)
#' plot(b, n_cores = 1)
#'
#' ##############################
#' # Example 3: Supervised mode #
#' ##############################
#' # Tune the model for explaining the politic block
#' # (politic connected to the two other blocks)
#' cv.out = rgcca_cv(blocks, response = 3, ncomp = 2, n_cores = 1)
#' print(cv.out)
#' plot(cv.out)
#'
#' fit.rgcca = rgcca(blocks, response = 3, ncomp = 2,
#'                   tau = cv.out$bestpenalties)
#' plot(fit.rgcca, type = "both")
#'
#' b = bootstrap(fit.rgcca, n_cores = 1, n_boot = 10)
#' plot(b, n_cores = 1)
#'
#' ##########################
#' # Example 4: Sparse GCCA #
#' ##########################
#'
#'
#' @export
#' @import ggplot2
#' @importFrom grDevices dev.off rgb colorRamp pdf colorRampPalette
#' @importFrom graphics plot
#' @importFrom stats cor quantile runif sd na.omit p.adjust pnorm qnorm weights
#' @importFrom utils read.table write.table
#' @importFrom scales hue_pal
#' @importFrom stats model.matrix
#' @importFrom methods is
#' @seealso \code{\link[RGCCA]{plot.rgcca}}, \code{\link[RGCCA]{print.rgcca}},
#' \code{\link[RGCCA]{rgcca_cv_k}},
#' \code{\link[RGCCA]{rgcca_permutation}}
#' \code{\link[RGCCA]{rgcca_predict}}

#TODO: Use out$blcoks instead of out$calls$blocks
rgcca <- function(
    blocks,
    type = "rgcca",
    scale = TRUE,
    scale_block = TRUE,
    prescaling = FALSE,
    connection = 1 - diag(length(blocks)),
    scheme = "factorial",
    ncomp = rep(1, length(blocks)),
    tau = rep(1, length(blocks)),
    sparsity = rep(1, length(blocks)),
    ranks = rep(1, length(blocks)),
    regularisation_matrices = NULL,
    kronecker_covariance = F,
    init = "svd",
    bias = TRUE,
    tol = 1e-08,
    response = NULL,
    superblock = FALSE,
    method = "nipals",
    verbose = FALSE,
    quiet = TRUE,
    penalty_coef = 0,
    n_run = 1,
    n_cores = 1)
{

    if(class(blocks)=="permutation")
    {
        message("All the parameters were imported from the fitted rgcca_permutation object")
        scale_block = blocks$call$scale_block
        scale = blocks$call$scale
        scheme = blocks$call$scheme
        connection = blocks$call$connection
        tol = blocks$call$tol
        method = blocks$call$method
        superblock = blocks$call$superblock
        if(blocks$call$par_type == "tau") tau = blocks$bestpenalties
        if(blocks$call$par_type == "ncomp") ncomp = blocks$bestpenalties
        if(blocks$call$par_type == "sparsity") sparsity = blocks$bestpenalties
        superblock <- blocks$call$superblock
        blocks <- blocks$blocks
    }
    if(class(blocks)=="cval")
    {
        message("All the parameters were imported from the fitted cval")
        scale_block = blocks$call$scale_block
        scale = blocks$call$scale
        scheme = blocks$call$scheme
        response = blocks$call$response
        tol = blocks$call$tol
        method = blocks$call$method
        if(blocks$call$par_type[[1]] == "tau") tau=blocks$bestpenalties
        if(blocks$call$par_type[[1]] == "ncomp") ncomp=blocks$bestpenalties
        if(blocks$call$par_type[[1]] == "sparsity")
          sparsity=blocks$bestpenalties
        blocks<-blocks$blocks
    }

    if(length(blocks) == 1){
        if(type != "pca")
        {
            type = "pca"
            message("type='rgcca' is not available for one block only and
                    type was converted to 'pca'.")
        }
    }

    if (any(sapply(blocks, function(x) length(dim(x))) > 2)) {
        if(!type %in% c("mgcca", "ns_mgcca", "gmgcca", "ns_mgcca_penalized", "gmgcca_penalized"))
        {
            message(paste0("type='", type, "' is not available for tensor blocks
                           so type was converted to 'mgcca'."))
            type = "mgcca"
        }
    }

    if (!missing(sparsity) && missing(type))
        type <- "sgcca"

    if (!missing(ranks) && missing(type))
        type <- "mgcca"

    if (!missing(regularisation_matrices) && missing(type))
        type <- "mgcca"

    if (!missing(connection) && missing(superblock))
        superblock <- FALSE

    if (!missing(response) && missing(superblock))
        superblock <- FALSE

    # if (!missing(superblock) && !(missing(response) || missing(connection)))


    if (tolower(type) %in% c("sgcca", "spca", "spls")) {
        if (!missing(tau) && missing(sparsity))
           stop_rgcca(paste0("sparsity parameters required for ",
                             tolower(type), " (instead of tau)."))
        gcca <- sgccaNa
        par <- "sparsity"
        penalty <- sparsity

    }else if (tolower(type) %in% c("mgcca")) {
        gcca <- mgccaNa
        par <- "tau"
        penalty <- tau

    }else if (tolower(type) %in% c("ns_mgcca")) {
      gcca <- ns_mgccaNa
      par <- "tau"
      penalty <- tau

    }else if (tolower(type) %in% c("ns_mgcca_penalized")) {
      gcca <- ns_mgcca_penalizedNa
      par <- "tau"
      penalty <- tau

    }else if (tolower(type) %in% c("gmgcca")) {
      gcca <- gmgccaNa
      par <- "tau"
      penalty <- tau

    }else if (tolower(type) %in% c("gmgcca_penalized")) {
      gcca <- gmgcca_penalizedNa
      par <- "tau"
      penalty <- tau

    } else {
        if (!missing(sparsity) & missing(tau))
           stop_rgcca(paste0("tau parameters required for ",
                             tolower(type), " (instead of sparsity)."))
        gcca <- rgccaNa
        par <- "tau"
        penalty <- tau
    }
    #if (superblock && any(penalty == "optimal"))
    #    stop_rgcca("Optimal tau is not available with superblock option.")

    if (type == "mgcca") {
      if (missing(method)) method <- "complete"
      ranks  <- check_ranks(ranks, blocks)
      regularisation_matrices <- check_reg_matrices(
        regularisation_matrices, blocks)
    }

    if (type %in% c("ns_mgcca", "ns_mgcca_penalized")) {
      if (missing(method)) method  <- "complete"
      ranks                <- check_ranks(ranks, blocks)
      kronecker_covariance <- check_boolean("kronecker_covariance", kronecker_covariance)
    }

    if (type == "gmgcca") {
      if (missing(method)) method  <- "complete"
      ranks                <- check_ranks(ranks, blocks)
    }

    if (type == "gmgcca_penalized") {
      if (missing(method)) method  <- "complete"
    }


    match.arg(init, c("svd", "random"))
    check_method(type)

  # Check blocks size, add NA for missing subjects
    blocks = check_blocks(blocks, add_NAlines=TRUE, n=1,
                          init=TRUE, quiet=quiet)
    if (!is.null(response))
        check_blockx("response", response, blocks)
    check_integer("tol", tol, float = TRUE, min = 0)


    for (i in c("superblock", "verbose", "scale", "bias", "quiet"))
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
    # TODO: should we keep raw? It seems to be used by cv, permutation and predict
    raw = blocks
   if(!is.null(response))
   {
       if(mode(blocks[[response]]) == "character")
       {
              if(length(unique(blocks[[response]])) == 1){
                stop("Only one level in the variable to predict")}
              blocks[[response]] = asDisjonctive(blocks[[response]])
       }

   }

    opt$blocks <- blocks
    opt$superblock <- check_superblock(response, opt$superblock, !quiet)
    opt$blocks <- set_superblock(opt$blocks, opt$superblock, type, !quiet)

    if (!is.null(response)) {
        response <- check_blockx("response", response, opt$blocks)
        }


    if (!is.matrix(opt$connection) || !is.null(response)) {
        opt$connection <- set_connection(
            opt$blocks,
            superblock=opt$superblock,response=response
        )
         opt$connection <- check_connection(opt$connection, opt$blocks)
    } else if (is.matrix(opt$connection)) {
        opt$connection <- check_connection(opt$connection, opt$blocks)
        opt$connection <- opt$connection[names(blocks), names(blocks)]
    }


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
            scale_block = scale_block,
            method = method,
            prescaling = prescaling,
            quiet=quiet,
            n_run = n_run,
            n_cores = n_cores
        )
    )

    if (type == "mgcca") {
      func$regularisation_matrices <- regularisation_matrices
      func$ranks                   <- ranks
    }

    if (type == "ns_mgcca") {
      func$ranks                   <- ranks
      func$kronecker_covariance    <- kronecker_covariance
    }

    if (type == "ns_mgcca_penalized") {
      func$ranks                   <- ranks
      func$kronecker_covariance    <- kronecker_covariance
      func$penalty_coef            <- penalty_coef
    }

    if (type == "gmgcca") {
      func$regularisation_matrices <- regularisation_matrices
      func$ranks                   <- ranks
    }

    if (type == "gmgcca_penalized") {
      func$penalty_coef <- penalty_coef
    }

    func[[par]] <- opt$penalty
    func_out <- eval(as.call(func))$rgcca

    # TODO: check if this is necessary
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
        connection = opt$connection,
        superblock = opt$superblock,
        ncomp = opt$ncomp,
        scheme = opt$scheme,
        raw=raw
    )

    is_optimal <- any(opt$penalty == "optimal")
    func_out$call[["optimal"]] <- is_optimal

    if(is_optimal){
        func_out$call[[par]] <- func_out$tau
    }else
        func_out$call[[par]] <- opt$penalty

    if(!is.null(func_out$tau))
        func_out$tau <- NULL

    if(type == "mgcca") func_out$call$ranks = ranks
    if(type == "ns_mgcca") func_out$call$ranks = ranks
    if(type == "ns_mgcca_penalized") func_out$call$ranks = ranks
    if(type == "gmgcca") func_out$call$ranks = ranks

    for (i in c(
        "scale",
        "init",
        "bias",
        "tol",
        "verbose",
        "response",
        "scale_block",
        "method",
        "type"
    ))
        func_out$call[[i]] <- as.list(environment())[[i]]

    class(func_out) <- "rgcca"
    invisible(func_out)
}
