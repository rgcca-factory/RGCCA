#' SGCCA extends RGCCA to address the issue of variable selection. Specifically,
#' RGCCA is combined with an L1-penalty that gives rise to Sparse GCCA (SGCCA)
#' which is implemented in the function sgcca().
#' Given \eqn{J} matrices \eqn{X_1, X_2, ..., X_J}, that represent
#' \eqn{J} sets of variables observed on the same set of \eqn{n} individuals.
#' The matrices \eqn{X_1, X_2, ..., X_J} must have the same number of rows, but
#' may (and usually will) have different numbers of columns. Blocks are not
#' necessarily fully connected within the SGCCA framework. Hence the use of
#' SGCCA requires the construction (user specified) of a design matrix
#' (\eqn{connection}) that characterizes the connections between blocks.
#' Elements of the symmetric design matrix \eqn{connection = (c_{jk})} are
#' equal to 1 if block \eqn{j} and block \eqn{k} are connected, and 0 otherwise.
#' The SGCCA algorithm is very similar to the RGCCA algorithm and keeps the same
#' monotone convergence properties (i.e. the bounded criteria to be maximized
#' increases at each step of the iterative procedure and hits at convergence
#' a stationary point).
#' Moreover, using a deflation strategy, sgcca() enables the computation of
#' several SGCCA block components (specified by ncomp) for each block. Block
#' components for each block are guaranteed to be orthogonal when using this
#' deflation strategy. The so-called symmetric deflation is considered in this
#' implementation, i.e. each block is deflated with respect to its own
#' component. Moreover, we stress that the numbers of components per block
#' could differ from one block to another.
#' @inheritParams rgcca
#' @inheritParams rgccad
#' @param na.rm If TRUE, runs sgcca only on available data.
#' @return \item{Y}{A list of \eqn{J} elements. Each element of \eqn{Y} is a
#' matrix that contains the analysis components for the corresponding block.}
#' @return \item{a}{A list of \eqn{J} elements. Each element of \eqn{a} is a
#' matrix that contains the outer weight vectors for each block.}
#' @return \item{astar}{A list of \eqn{J} elements. Each column of astar[[j]] is
#' a vector such that Y[[j]][, h] = blocks[[j]] \%*\% astar[[j]][, h].}
#' @return \item{crit}{A vector of integer that contains for each component the
#' values of the analysis criteria across iterations.}
#' @references Tenenhaus, A., Philippe, C., Guillemot, V., Le Cao, K. A.,
#' Grill, J., and Frouin, V. , "Variable selection for generalized canonical
#' correlation analysis.," Biostatistics, vol. 15, no. 3, pp. 569-583, 2014.
#' @title Variable Selection For Generalized Canonical Correlation Analysis
#' (SGCCA)
#' @examples
#' #############
#' # Example 1 #
#' #############
#' \dontrun{
#' # Download the dataset's package at http://biodev.cea.fr/sgcca/.
#' # --> gliomaData_0.4.tar.gz
#'
#' data("ge_cgh_locIGR", package = "gliomaData")
#'
#' blocks <- ge_cgh_locIGR$multiblocks
#' Loc <- factor(ge_cgh_locIGR$y)
#' levels(Loc) <- colnames(ge_cgh_locIGR$multiblocks$y)
#' connection <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
#' tau <- c(1, 1, 0)
#'
#' # rgcca algorithm using the dual formulation for X1 and X2
#' # and the dual formulation for X3
#' blocks[[3]] <- blocks[[3]][, -3]
#' result.rgcca <- rgcca(
#'   blocks = blocks, connection = connection, tau = tau,
#'   ncomp = c(2, 2, 1), scheme = "factorial",
#'   verbose = TRUE, method = "rgcca"
#' )
#' # sgcca algorithm
#' result.sgcca <- rgcca(
#'   blocks = blocks, connection = connection,
#'   sparsity = c(.071, .2, 1), ncomp = c(2, 2, 1),
#'   scheme = "centroid", verbose = TRUE, method = "sgcca"
#' )
#'
#' ############################
#' # plot(y1, y2) for (RGCCA) #
#' ############################
#' layout(t(1:2))
#' plot(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[2]][, 1],
#'   col = "white", xlab = "Y1 (GE)", ylab = "Y2 (CGH)",
#'   main = "Factorial plan of RGCCA"
#' )
#' text(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[2]][, 1],
#'   Loc,
#'   col = as.numeric(Loc), cex = .6
#' )
#' plot(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[1]][, 2],
#'   col = "white", xlab = "Y1 (GE)",
#'   ylab = "Y2 (GE)", main = "Factorial plan of RGCCA"
#' )
#' text(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[1]][, 2],
#'   Loc,
#'   col = as.numeric(Loc), cex = .6
#' )
#'
#' ############################
#' # plot(y1, y2) for (SGCCA) #
#' ############################
#' layout(t(1:2))
#' plot(result.sgcca$Y[[1]][, 1], result.sgcca$Y[[2]][, 1],
#'   col = "white", xlab = "Y1 (GE)",
#'   ylab = "Y2 (CGH)", main = "Factorial plan of SGCCA"
#' )
#' text(result.sgcca$Y[[1]][, 1], result.sgcca$Y[[2]][, 1],
#'   Loc,
#'   col = as.numeric(Loc), cex = .6
#' )
#'
#' plot(result.sgcca$Y[[1]][, 1], result.sgcca$Y[[1]][, 2],
#'   col = "white", xlab = "Y1 (GE)",
#'   ylab = "Y2 (GE)", main = "Factorial plan of SGCCA"
#' )
#' text(result.sgcca$Y[[1]][, 1], result.sgcca$Y[[1]][, 2],
#'   Loc,
#'   col = as.numeric(Loc), cex = .6
#' )
#'
#' # sgcca algorithm with multiple components and different
#' # L1 penalties for each components
#' # (-> sparsity is a matrix)
#' init <- "random"
#' result.sgcca <- sgcca(blocks, connection,
#'   sparsity = matrix(c(.071, .2, 1, 0.06, 0.15, 1),
#'     nrow = 2, byrow = TRUE
#'   ),
#'   ncomp = c(2, 2, 1), scheme = "factorial", bias = TRUE,
#'   init = init, verbose = TRUE
#' )
#' # number of non zero elements per dimension
#' apply(result.sgcca$a[[1]], 2, function(x) sum(x != 0))
#' # (-> 145 non zero elements for a11 and 107 non zero elements for a12)
#' apply(result.sgcca$a[[2]], 2, function(x) sum(x != 0))
#' # (-> 85 non zero elements for a21 and 52 non zero elements for a22)
#' init <- "svd"
#' result.sgcca <- sgcca(blocks, connection,
#'   sparsity = matrix(c(.071, .2, 1, 0.06, 0.15, 1),
#'     nrow = 2, byrow = TRUE
#'   ),
#'   ncomp = c(2, 2, 1), scheme = "factorial",
#'   bias = TRUE,
#'   init = init, verbose = TRUE
#' )
#' }
#' @noRd

sgcca <- function(blocks, connection = 1 - diag(length(blocks)),
                  sparsity = rep(1, length(blocks)),
                  ncomp = rep(1, length(blocks)), scheme = "centroid",
                  init = "svd", bias = TRUE, tol = .Machine$double.eps,
                  verbose = FALSE, na.rm = TRUE,
                  superblock = FALSE, response = NULL,
                  disjunction = NULL, n_iter_max = 1000,
                  comp_orth = TRUE) {
  if (verbose) {
    scheme_str <- ifelse(is(scheme, "function"), "user-defined", scheme)
    cat(
      "Computation of the SGCCA block components based on the",
      scheme_str, "scheme \n"
    )
  }

  ##### Initialization #####
  # ndefl number of deflation per block
  ndefl <- ncomp - 1
  N <- max(ndefl)
  J <- length(blocks)
  pjs <- vapply(blocks, NCOL, numeric(1L))
  nb_ind <- NROW(blocks[[1]])

  crit <- list()
  R <- blocks

  a <- lapply(seq(J), function(b) c())
  Y <- lapply(seq(J), function(b) c())

  if (superblock && comp_orth) {
    P <- c()
  } else {
    P <- lapply(seq(J), function(b) c())
  }

  if (is.vector(sparsity)) {
    sparsity <- matrix(
      rep(sparsity, N + 1),
      nrow = N + 1, J, byrow = TRUE
    )
  }

  ##### Computation of SGCCA components #####
  for (n in seq(N + 1)) {
    if (verbose) {
      cat(paste0(
        "Computation of the SGCCA block components #", n,
        " is under progress...\n"
      ))
    }
    gcca_result <- sgccak(R, connection,
      sparsity = sparsity[n, ], scheme = scheme,
      init = init, bias = bias, tol = tol,
      verbose = verbose, na.rm = na.rm,
      response = response, disjunction = disjunction,
      n_iter_max = n_iter_max
    )

    # Store crit
    crit[[n]] <- gcca_result$crit

    # Store Y, a, factors and weights
    a <- lapply(seq(J), function(b) cbind(a[[b]], gcca_result$a[[b]]))
    Y <- lapply(seq(J), function(b) cbind(Y[[b]], gcca_result$Y[, b]))

    # Deflation procedure
    if (n == N + 1) break
    defl_result <- deflate(gcca_result$a, gcca_result$Y, R, P, ndefl, n,
                           superblock, comp_orth, response, na.rm)
    R <- defl_result$R
    P <- defl_result$P
  }

  ##### Generation of the output #####
  if (N == 0) crit <- unlist(crit)

  astar <- compute_astar(a, P, superblock, comp_orth, N)

  out <- list(Y = Y, a = a, astar = astar, crit = crit)

  class(out) <- "sgcca"
  return(out)
}
