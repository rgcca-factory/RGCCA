#' netSGCCA extends SGCCA to include a Graph Penalty. Specifically,
#' netSGCCA is equivalent to RGCCA combined with an L1-constraint
#' and a Graph Penalty. It is implemented in the function netsgcca().
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
#' Moreover, using a deflation strategy, netsgcca() enables the computation of
#' several netSGCCA block components (specified by ncomp) for each block. Block
#' components for each block are guaranteed to be orthogonal when using this
#' deflation strategy. The so-called symmetric deflation is considered in this
#' implementation, i.e. each block is deflated with respect to its own
#' component. Moreover, we stress that the numbers of components per block
#' could differ from one block to another.
#' @inheritParams rgcca
#' @inheritParams rgccad
#' @param lambda vector of graph penalties
#' @param graph_laplacians list of graph laplacians
#' @param na.rm If TRUE, runs netsgcca only on available data.
#' @return \item{Y}{A list of \eqn{J} elements. Each element of \eqn{Y} is a
#' matrix that contains the analysis components for the corresponding block.}
#' @return \item{a}{A list of \eqn{J} elements. Each element of \eqn{a} is a
#' matrix that contains the outer weight vectors for each block.}
#' @return \item{astar}{A list of \eqn{J} elements. Each column of astar[[j]] is
#' a vector such that Y[[j]][, h] = blocks[[j]] \%*\% astar[[j]][, h].}
#' @return \item{crit}{A vector of integer that contains for each component the
#' values of the analysis criteria across iterations.}
#' @noRd

netsgcca <- function(blocks, connection = 1 - diag(length(blocks)),
                     sparsity = rep(1, length(blocks)),
                     lambda = rep(1, length(blocks)),
                     graph_laplacians,
                     ncomp = rep(1, length(blocks)), scheme = "centroid",
                     init = "svd", bias = TRUE, tol = .Machine$double.eps,
                     verbose = FALSE, na.rm = TRUE,
                     superblock = FALSE, response = NULL,
                     disjunction = NULL, n_iter_max = 1000,
                     comp_orth = TRUE) {
  if (verbose) {
    scheme_str <- ifelse(is(scheme, "function"), "user-defined", scheme)
    cat(
      "Computation of the netSGCCA block components based on the",
      scheme_str, "scheme \n"
    )
  }

  ##### Initialization #####
  # ndefl number of deflation per block
  ndefl <- ncomp - 1
  N <- max(ndefl)
  J <- length(blocks)

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

  ##### Computation of netSGCCA components #####
  for (n in seq(N + 1)) {
    if (verbose) {
      cat(paste0(
        "Computation of the netSGCCA block components #", n,
        " is under progress...\n"
      ))
    }
    gcca_result <- netsgccak(R, connection,
      sparsity = sparsity[n, ],
      lambda = lambda, graph_laplacians = graph_laplacians,
      scheme = scheme,
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

  class(out) <- "netsgcca"
  return(out)
}
