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
#' @param blocks A list that contains the J blocks of variables
#' \eqn{\mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}}{X1, X2, ..., XJ}.
#' Block \eqn{\mathbf{X}_j}{Xj} is a matrix of dimension
#' \eqn{n \times p_j}{n x p_j} where n is the number of
#' observations and \eqn{p_j} the number of variables.
#' @param method A character string indicating the multi-block component
#' method to consider: rgcca, sgcca, pca, spca, pls, spls, cca,
#' ifa, ra, gcca, maxvar, maxvar-b, maxvar-a, mcoa,cpca-1, cpca-2,
#' cpca-4, hpca, maxbet-b, maxbet, maxdiff-b, maxdiff,
#' sabscor, ssqcor, ssqcov-1, ssqcov-2, ssqcov, sumcor,
#' sumcov-1, sumcov-2, sumcov, sabscov-1, sabscov-2.
#' @param scale Logical value indicating if blocks are standardized.
#' @param scale_block Value indicating if each block is divided by
#' a constant value. If TRUE or "inertia", each block is divided by the
#' sum of eigenvalues of its empirical covariance matrix.
#' If "lambda1", each block is divided by the square root of the highest
#' eigenvalue of its empirical covariance matrix.
#' Otherwise the blocks are not scaled. If standardization is
#' applied (scale = TRUE), the block scaling is applied on the result of the
#' standardization.
#' @param connection  A symmetric matrix (\eqn{J \times J}{J x J}) that
#' describes the relationships between blocks.
#' @param scheme Character string or a function giving the scheme function for
#' covariance maximization among "horst" (the identity function), "factorial"
#'  (the squared values), "centroid" (the absolute values). The scheme function
#'  can be any continously differentiable convex function and it is possible to
#'  design explicitely the scheme function (e.g. function(x) x^4) as argument of
#'  rgcca function.  See (Tenenhaus et al, 2017) for details.
#' @param ncomp Vector of length J indicating the number of block components
#' for each block.
#' @param tau Either a \eqn{1 \times J}{1 x J} vector or a
#' \eqn{\mathrm{max}(ncomp) \times J}{max(ncomp) x J} matrix containing
#' the values of the regularization parameters (default: tau = 1, for each
#' block and each dimension). The regularization parameters varies from 0
#' (maximizing the correlation) to 1 (maximizing the covariance). If
#' tau = "optimal" the regularization parameters are estimated for each block
#' and each dimension using the Schafer and Strimmer (2005) analytical formula.
#' If tau is a \eqn{1 \times J}{1 x J} vector, tau[j] is identical across the
#' dimensions of block \eqn{\mathbf{X}_j}{Xj}. If tau is a matrix, tau[k, j]
#' is associated with \eqn{\mathbf{X}_{jk}}{Xjk} (kth residual matrix for
#' block j). The regularization parameters can also be estimated using
#' \link{rgcca_permutation} or \link{rgcca_cv}.
#' @param sparsity Either a \eqn{1*J} vector or a \eqn{max(ncomp) * J} matrix
#' encoding the L1 constraints applied to the outer weight vectors. The amount
#' of sparsity varies between \eqn{1/sqrt(p_j)} and 1 (larger values of sparsity
#' correspond to less penalization). If sparsity is a vector, L1-penalties are
#' the same for all the weights corresponding to the same block but different
#' components:
#' \deqn{for all h, |a_{j,h}|_{L_1} \le c_1[j] \sqrt{p_j},}
#' with \eqn{p_j} the number of variables of \eqn{X_j}.
#' If sparsity is a matrix, each row \eqn{h} defines the constraints applied to
#' the weights corresponding to components \eqn{h}:
#' \deqn{for all h, |a_{j,h}|_{L_1} \le c_1[h,j] \sqrt{p_j}.} It can be
#' estimated by using \link{rgcca_permutation}.
#' @param init Character string giving the type of initialization to use in
#' the  algorithm. It could be either by Singular Value Decompostion ("svd")
#' or by random initialisation ("random") (default: "svd").
#' @param bias A logical value for biaised (\eqn{1/n}) or unbiaised
#' (\eqn{1/(n-1)}) estimator of the var/cov (default: bias = TRUE).
#' @param tol The stopping value for the convergence of the algorithm.
#' @param response Numerical value giving the position of the response block.
#' When the response argument is filled the supervised mode is automatically
#' activated.
#' @param superblock Boolean indicating the presence of a superblock
#' (deflation strategy must be adapted when a superblock is used).
#' @param NA_method  Character string corresponding to the method used for
#' handling missing values ("nipals", "complete"). (default: "nipals").
#' \itemize{
#' \item{\code{"complete"}}{corresponds to perform RGCCA on the fully observed
#' observations (observations with missing values are removed)}
#' \item{\code{"nipals"}}{corresponds to perform RGCCA algorithm on available
#' data (NIPALS-type algorithm)}}
#' @param verbose Logical value indicating if the progress of the
#' algorithm is reported while computing.
#' @param quiet Logical value indicating if warning messages are reported.
#' @param n_iter_max Integer giving the algorithm's maximum number of
#' iterations.
#' @return A rgcca fitted object
#' @return \item{Y}{List of \eqn{J} elements. Each element of the list \eqn{Y}
#' is a matrix that contains the RGCCA block components for the corresponding
#' block.}
#' @return \item{a}{List of \eqn{J} elements. Each element of the list \eqn{a}
#' is a matrix of block weight vectors for the corresponding block.}
#' @return \item{astar}{List of \eqn{J} elements. Each column of astar[[j]] is a
#' vector such that Y[[j]][, h] = blocks[[j]] \%*\% astar[[j]][, h].}
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
#' blocks <- list(
#'   agriculture = Russett[, seq(3)],
#'   industry = Russett[, 4:5],
#'   politic = Russett[, 6:11]
#' )
#'
#' # Blocks are fully connected, factorial scheme and tau =1 for all blocks is
#' # used by default
#' fit.rgcca <- rgcca(
#'   blocks = blocks, method = "rgcca", connection = 1 - diag(3),
#'   scheme = "factorial", tau = rep(1, 3)
#' )
#' print(fit.rgcca)
#' plot(fit.rgcca, type = "weight", block = 3)
#' politic <- as.vector(apply(Russett[, 9:11], 1, which.max))
#' plot(fit.rgcca,
#'   type = "sample", block = 1:2,
#'   comp = rep(1, 2), resp = politic
#' )
#'
#' ############################################
#' # Example 2: RGCCA and multiple components #
#' ############################################
#' fit.rgcca <- rgcca(blocks,
#'   method = "rgcca",
#'   connection = 1 - diag(3), superblock = FALSE,
#'   tau = rep(1, 3), ncomp = c(2, 2, 2),
#'   scheme = "factorial", verbose = TRUE
#' )
#'
#' politic <- as.vector(apply(Russett[, 9:11], 1, which.max))
#' plot(fit.rgcca,
#'   type = "sample", block = 1:2,
#'   comp = rep(1, 2), resp = politic
#' )
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
#' perm.out <- rgcca_permutation(blocks,
#'   n_cores = 1,
#'   par_type = "sparsity", n_perms = 10
#' )
#' print(perm.out)
#' plot(perm.out)
#'
#' fit.sgcca <- rgcca(blocks, sparsity = perm.out$bestpenalties)
#' plot(fit.sgcca, type = "ave")
#'
#' # Select the most significant variables
#' b <- bootstrap(fit.sgcca, n_cores = 1, n_boot = 100)
#' plot(b, n_cores = 1)
#'
#' ##############################
#' # Example 3: Supervised mode #
#' ##############################
#' # Tune the model for explaining the politic block
#' # (politic connected to the two other blocks)
#' cv.out <- rgcca_cv(blocks, response = 3, ncomp = 2, n_cores = 1)
#' print(cv.out)
#' plot(cv.out)
#'
#' fit.rgcca <- rgcca(blocks,
#'   response = 3, ncomp = 2,
#'   tau = cv.out$bestpenalties
#' )
#' plot(fit.rgcca, type = "both")
#'
#' b <- bootstrap(fit.rgcca, n_cores = 1, n_boot = 10)
#' plot(b, n_cores = 1)
#' }
#'
#' @export
#' @importFrom graphics plot
#' @importFrom stats cor quantile sd p.adjust rnorm pnorm qnorm
#' @importFrom stats model.matrix
#' @importFrom methods is
#' @seealso \code{\link[RGCCA]{plot.rgcca}}, \code{\link[RGCCA]{print.rgcca}},
#' \code{\link[RGCCA]{rgcca_cv}},
#' \code{\link[RGCCA]{rgcca_permutation}}
#' \code{\link[RGCCA]{rgcca_predict}}
rgcca <- function(blocks, method = "rgcca",
                  scale = TRUE, scale_block = "inertia",
                  connection = 1 - diag(length(blocks)),
                  scheme = "factorial",
                  ncomp = rep(1, length(blocks)),
                  tau = rep(1, length(blocks)),
                  sparsity = rep(1, length(blocks)),
                  init = "svd", bias = TRUE, tol = 1e-08,
                  response = NULL,
                  superblock = FALSE,
                  NA_method = "nipals", verbose = FALSE, quiet = TRUE,
                  n_iter_max = 1000) {
  rgcca_args <- as.list(environment())
  ### If specific objects are given for blocks, parameters are imported from
  #   these objects.
  rgcca_args <- get_rgcca_args(blocks, rgcca_args)
  rgcca_args$quiet <- quiet
  rgcca_args$verbose <- verbose

  ### Check parameters
  match.arg(rgcca_args$init, c("svd", "random"))
  rgcca_args$blocks <- check_blocks(rgcca_args$blocks,
    add_NAlines = TRUE, quiet = rgcca_args$quiet,
    response = rgcca_args$response
  )

  blocks <- remove_null_sd(rgcca_args$blocks)$list_m

  check_integer("tol", rgcca_args$tol, float = TRUE, min = 0)
  check_integer("n_iter_max", rgcca_args$n_iter_max, min = 1)
  for (i in c("superblock", "verbose", "scale", "bias", "quiet")) {
    check_boolean(rgcca_args[[i]], get(i))
  }

  rgcca_args$tau <- elongate_arg(rgcca_args$tau, blocks)
  rgcca_args$ncomp <- elongate_arg(rgcca_args$ncomp, blocks)
  rgcca_args$sparsity <- elongate_arg(rgcca_args$sparsity, blocks)

  ### Get last parameters based on the method
  tmp <- select_analysis(rgcca_args, blocks)
  opt <- tmp$opt
  rgcca_args <- tmp$rgcca_args

  ### One hot encode the response block if needed
  disjunction <- NULL
  if (!is.null(rgcca_args$response)) {
    blocks[[rgcca_args$response]] <- as_disjunctive(
      blocks[[rgcca_args$response]]
    )
    disjunction <- attributes(blocks[[rgcca_args$response]])$disjunction
  }
  # Change penalty to 0 if there is a univariate disjunctive block response
  if (isTRUE(disjunction)) {
    if (is.matrix(rgcca_args[[opt$param]])) {
      rgcca_args[[opt$param]][, rgcca_args$response] <- 0
    } else {
      rgcca_args[[opt$param]][rgcca_args$response] <- 0
    }
  }

  ### Apply strategy to deal with NA, scale and prepare superblock
  tmp <- handle_NA(blocks, NA_method = NA_method)
  na.rm <- tmp$na.rm
  blocks <- scaling(tmp$blocks,
    scale = rgcca_args$scale,
    bias = rgcca_args$bias,
    scale_block = rgcca_args$scale_block
  )
  if (rgcca_args$superblock) {
    blocks[["superblock"]] <- Reduce(cbind, blocks)
    colnames(blocks[["superblock"]]) <- paste0(
      "s-", colnames(blocks[["superblock"]])
    )
  }

  ### Call the gcca function
  gcca_args <- rgcca_args[c(
    "connection", "ncomp", "scheme", "init", "bias", "tol",
    "verbose", "quiet", "superblock", "response", "n_iter_max"
  )]
  gcca_args[["na.rm"]] <- na.rm
  gcca_args[["blocks"]] <- blocks
  gcca_args[["disjunction"]] <- disjunction
  gcca_args[[opt$param]] <- rgcca_args[[opt$param]]
  func_out <- do.call(opt$gcca, gcca_args)

  ### Format the output
  func_out <- format_output(func_out, rgcca_args, opt, blocks, disjunction)

  class(func_out) <- "rgcca"
  invisible(func_out)
}
