#' Regularized (or Sparse) Generalized Canonical Correlation Analysis (S/RGCCA)
#'
#' RGCCA is a generalization of regularized canonical correlation analysis to
#' three or more sets of variables. SGCCA extends RGCCA to address the issue of
#' variable selection
#'
#' @details
#' Given J matrices X1, X2, ..., XJ that represent J sets of variables observed
#' on the same set of n individuals. The matrices X1, X2, ..., XJ must have the
#' same number of rows, but may (and usually will) have different
#' numbers of columns. The aim of RGCCA is to study the relationships between
#' these J blocks of variables. It constitutes a general framework for
#' many multi-block data analysis methods (see Tenenhaus and Tenenhaus, 2011 ;
#' Tenenhaus et al. 2017). It combines the power of multi-block data analysis
#' methods (maximization of well identified criteria) and the flexibility of
#' PLS path modeling (the researcher decides which blocks are connected and
#' which are not). Hence, the use of RGCCA requires the construction (user
#' specified) of a design matrix C that characterizes the connections between
#' blocks. Elements of the (symmetric) design matrix C = (c_{jk}) are positive
#' (and usually equal to 1 if block j and block k are connected, and 0
#' otherwise). The rgcca() function implements a monotone global convergent
#' algorithm - i.e. the bounded criteria to be maximized increases at each
#' step of the iterative procedure and hits, at convergence a stationary point
#' of the RGCCA optimization problem. Moreover, depending on the dimensionality
#' of each block Xj, j = 1, \ldots, J, the primal (when n > p_j) algorithm or
#' the dual (when n < p_j) algorithm is used (see Tenenhaus et al. 2015). At
#' last, a deflation strategy is used to compute several RGCCA block components
#' (specified by ncomp) for each block. Block components of each block are
#' guaranteed to be orthogonal. The so-called symmetric deflation is implemented
#' (i.e. each block is deflated with respect to its own component). It should be
#' noted that the numbers of components per block can differ from one block to
#' another. SGCCA extends RGCCA to address the issue of variable selection
#' (Tenenhaus et al, 2014). Specifically, RGCCA is combined with an L1-penalty
#' that gives rise to Sparse GCCA (SGCCA). The SGCCA algorithm is very similar
#' to the RGCCA algorithm and keeps the same convergence properties (i.e. the
#' bounded criteria to be maximized increases at each step of the iterative
#' procedure and hits at convergence a stationary point). Moreover, using a
#' deflation strategy,  sgcca() enables the computation of several SGCCA
#' orthogonal block components (specified by ncomp) for each block. The rgcca()
#' function can handle missing values using a NIPALS type algorithm (non-linear
#' iterative partial least squares algorithm) described in (Tenenhaus et al,
#' 2005). Guidelines describing how to use RGCCA in practice are provided in
#' (Garali et al., 2018).
#' @inheritParams rgccad
#' @inheritParams sgcca
#' @inheritParams select_analysis
#' @param scale Logical value indicating if blocks are standardized.
#' @param scale_block Logical value indicating if each block is divided by
#' the square root of its number of variables.
#' @param NA_method  Character string corresponding to the method used for
#' handling missing values ("nipals", "complete"). (default: "nipals").
#' \itemize{
#' \item{\code{"complete"}}{corresponds to perform RGCCA on the fully observed
#' observations (observations with missing values are removed)}
#' \item{\code{"nipals"}}{corresponds to perform RGCCA algorithm on available
#' data (NIPALS-type algorithm)}}
#' @return A rgcca fitted object
#' @return \item{Y}{List of \eqn{J} elements. Each element of the list \eqn{Y}
#' is a matrix that contains the RGCCA block components for the corresponding
#' block.}
#' @return \item{a}{List of \eqn{J} elements. Each element of the list \eqn{a}
#' is a matrix of block weight vectors for the corresponding block.}
#' @return \item{astar}{List of \eqn{J} elements. Each element of astar is a
#' matrix defined as Y[[j]][, h] = A[[j]]\%*\%astar[[j]][, h].}
#' @return \item{tau}{Regularization parameters used during the analysis.}
#' @return \item{crit}{List of vector of length max(ncomp). Each vector of
#' the list is related to one specific deflation stage and reports the values
#' of the criterion for this stage across iterations.}
#' @return \item{primal_dual}{A \eqn{1 \times J} vector that contains the
#' formulation ("primal" or "dual") applied to each of the \eqn{J} blocks
#' within the RGCCA alogrithm.}
#' @return \item{AVE}{List of numerical values giving the indicators of model
#' quality based on the Average Variance Explained (AVE): AVE(for each block),
#' AVE(outer model), AVE(inner model).}
#' @return \item{A}{List that contains the J blocks of variables X1, X2, ...,
#' XJ. Block Xj is a matrix of dimension n x p_j where p_j is the number of
#' variables in X_j. These blocks are imputed when an imputation strategy is
#' selected.}
#' @return \item{call}{Call of the function.}
#' @references Garali I, Adanyeguh IM, Ichou F, Perlbarg V, Seyer A, Colsch B,
#' Moszer I, Guillemot V, Durr A, Mochel F, Tenenhaus A. (2018) A strategy for
#' multimodal data integration: application to biomarkers identification
#' in spinocerebellar ataxia. Briefings in Bioinformatics. 19(6):1356-1369.
#' @references Tenenhaus M., Tenenhaus A. and Groenen P. J. (2017). Regularized
#' generalized canonical correlation analysis: a framework for sequential
#' multiblock component methods. Psychometrika, 82(3), 737-777.
#' @references Tenenhaus A., Philippe C. and Frouin, V. (2015). Kernel
#' generalized canonical correlation analysis. Computational Statistics and
#' Data Analysis, 90, 114-131.
#' @references Tenenhaus A., Philippe C., Guillemot V., Le Cao K. A., Grill J.
#' and Frouin, V. (2014), Variable selection for generalized canonical
#' correlation analysis, Biostatistics, 15(3), pp. 569-583.
#' @references Tenenhaus A. and Tenenhaus M., (2011). Regularized Generalized
#' Canonical Correlation Analysis, Psychometrika, 76(2), pp 257-284.
#' @references Schafer J. and Strimmer K. (2005). A shrinkage approach to
#' large-scale covariance matrix estimation and implications for functional
#' genomics. Statistical Applications in Genetics and Molecular Biology 4:32.
#' @references Arnaud Gloaguen, Vincent Guillemot, Arthur Tenenhaus.
#' An efficient algorithm to satisfy l1 and l2 constraints.
#' 49emes Journees de Statistique, May 2017, Avignon, France. (hal-01630744)
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
#' fit.rgcca = rgcca(blocks=blocks, method = "rgcca", connection = 1-diag(3),
#'                   scheme = "factorial", tau = rep(1, 3))
#' print(fit.rgcca)
#' plot(fit.rgcca, type = "weight", block = 3)
#' politic = as.vector(apply(Russett[, 9:11], 1, which.max))
#' plot(fit.rgcca, type = "sample", block = 1:2,
#'      comp = rep(1, 2), resp = politic)
#'
#' ############################################
#' # Example 2: RGCCA and multiple components #
#' ############################################
#' fit.rgcca = rgcca(blocks, method = "rgcca",
#'                   connection = C, superblock = FALSE,
#'                   tau = rep(1, 3), ncomp = c(2, 2, 2),
#'                   scheme = "factorial", verbose = TRUE)
#'
#' politic = as.vector(apply(Russett[, 9:11], 1, which.max))
#' plot(fit.rgcca, type = "sample", block = 1:2,
#'      comp = rep(1, 2), resp = politic)
#'
#' plot(fit.rgcca, type = "ave")
#' plot(fit.rgcca, type = "weight", block = 1)
#' plot(fit.rgcca, type = "loadings")
#' \dontrun{
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
#' }
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
rgcca <- function(blocks, method = "rgcca",
                  scale = TRUE, scale_block = TRUE,
                  connection = 1 - diag(length(blocks)),
                  scheme = "factorial",
                  ncomp = rep(1, length(blocks)),
                  tau = rep(1, length(blocks)),
                  sparsity = rep(1, length(blocks)),
                  init = "svd", bias = TRUE, tol = 1e-08,
                  response = NULL, superblock = FALSE,
                  NA_method = "nipals", verbose = FALSE, quiet = TRUE){

    if(class(blocks)=="permutation")
    {
        message("All the parameters were imported from the fitted rgcca_permutation object.")
        scale_block = blocks$call$scale_block
        scale = blocks$call$scale
        scheme = blocks$call$scheme
        connection = blocks$call$connection
        tol = blocks$call$tol
        NA_method = blocks$call$NA_method
        superblock = blocks$call$superblock
        if(blocks$call$par_type == "tau") tau = blocks$bestpenalties
        if(blocks$call$par_type == "ncomp") ncomp = blocks$bestpenalties
        if(blocks$call$par_type == "sparsity") sparsity = blocks$bestpenalties
        superblock <- blocks$call$superblock
        blocks <- blocks$call$blocks
    }
    if(class(blocks)=="cval")
    {
        message("All the parameters were imported from the fitted cval object.")
        scale_block = blocks$call$scale_block
        scale = blocks$call$scale
        scheme = blocks$call$scheme
        response = blocks$call$response
        tol = blocks$call$tol
        NA_method = blocks$call$NA_method
        if(blocks$call$par_type[[1]] == "tau") tau=blocks$bestpenalties
        if(blocks$call$par_type[[1]] == "ncomp") ncomp=blocks$bestpenalties
        if(blocks$call$par_type[[1]] == "sparsity")
          sparsity=blocks$bestpenalties
        blocks<-blocks$call$blocks
    }

    if(length(blocks) == 1){
        if(method != "pca")
        {
            method = "pca"
            message("method='rgcca' is not available for one block only and
                    method was converted to 'pca'.")
        }

    }

    if (!missing(sparsity) && missing(method))
        method <- "sgcca"

    if (!missing(connection) && missing(superblock))
        superblock <- FALSE

    if (!missing(response) && missing(superblock))
        superblock <- FALSE

    if (tolower(method) %in% c("sgcca", "spca", "spls")) {
      if (!missing(tau) && missing(sparsity))
           stop_rgcca(paste0("sparsity parameters required for ",
                             tolower(method), " (instead of tau)."))
        gcca <- sgcca
        par <- "sparsity"
        penalty <- sparsity

    }else{
        if (!missing(sparsity) & missing(tau))
           stop_rgcca(paste0("tau parameters required for ",
                             tolower(method), " (instead of sparsity)."))
        gcca <- rgccad
        par <- "tau"
        penalty <- tau
    }

    match.arg(init, c("svd", "random"))
    check_method(method)
    check_scheme(scheme)

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
        method  = method,
        quiet = quiet,
        response = response
    )
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

    opt$blocks     <- handle_NA(blocks, NA_method = NA_method)
    opt$blocks     <- scaling(opt$blocks,
                              scale = scale,
                              bias = bias,
                              scale_block = scale_block)
    opt$superblock <- check_superblock(response, opt$superblock, !quiet)
    opt$blocks     <- set_superblock(opt$blocks, opt$superblock,
                                     method, !quiet)

    if (NA_method == "nipals" &&
        Reduce("||", lapply(opt$blocks, function(x) any(is.na(x))))) {
      na.rm = TRUE
    } else {
      na.rm = FALSE
    }

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


    opt$penalty <- check_penalty(opt$penalty, opt$blocks, method)
    opt$ncomp <- check_ncomp(opt$ncomp, opt$blocks)

    warn_on <- FALSE

    if (any(sapply(opt$blocks, NCOL) > 1000)) {
            # if( (method <-<- "sgcca" && tau > 0.3) || method !<- "sgcca" )
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
            init = init,
            bias = bias,
            tol = tol,
            quiet = quiet,
            na.rm = na.rm,
            superblock=superblock
        )
    )

    func[[par]] <- opt$penalty
    func_out <- eval(as.call(func))

    for (i in c("a", "astar", "Y")) {
        names(func_out[[i]]) <- names(opt$blocks)
        for (j in seq(length(opt$blocks))) {
            if (i %in%  c("a", "astar") && NCOL(opt$blocks[[j]]) == 1)
                row.names(func_out[[i]][[j]]) <- colnames(opt$blocks[[j]])
        }
    }
    names(func_out$AVE$AVE_X) <- names(opt$blocks)

    class(func_out) <- tolower(method)
    func_out$call <- list(
        blocks = opt$blocks,
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

    for (i in c(
        "scale",
        "init",
        "bias",
        "tol",
        "verbose",
        "response",
        "scale_block",
        "NA_method",
        "method"
    ))
        func_out$call[[i]] <- as.list(environment())[[i]]

    class(func_out) <- "rgcca"
    invisible(func_out)
}
