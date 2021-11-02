#' NGCCA extends RGCCA with nonnegative constraints.
#' Given \eqn{J} matrices \eqn{X_1, X_2, ..., X_J}, that represent
#' \eqn{J} sets of variables observed on the same set of \eqn{n} individuals.
#' The matrices \eqn{X_1, X_2, ..., X_J} must have the same number of rows, but
#' may (and usually will) have different numbers of columns. Blocks are not
#' necessarily fully connected within the NGCCA framework. Hence the use of
#' NGCCA requires the construction (user specified) of a design matrix (\eqn{connection})
#' that characterizes the connections between blocks. Elements of the symmetric
#' design matrix \eqn{connection = (c_{jk})} are equal to 1 if block \eqn{j} and block
#' \eqn{k} are connected, and 0 otherwise. The NGCCA algorithm is very similar
#' to the RGCCA algorithm and keeps the same monotone convergence properties
#' (i.e. the bounded criteria to be maximized increases at each step of the
#' iterative procedure and hits at convergence a stationary point).
#' Moreover, using a deflation strategy, ngcca() enables the computation of
#' several NGCCA block components (specified by ncomp) for each block. Block
#' components for each block are guaranteed to be orthogonal when using this
#' deflation strategy. The so-called symmetric deflation is considered in this
#' implementation, i.e. each block is deflated with respect to its own
#' component. Moreover, we stress that the numbers of components per block
#' could differ from one block to another.
#' @inheritParams select_analysis
#' @inheritParams rgccad
#' @return \item{Y}{A list of \eqn{J} elements. Each element of \eqn{Y} is a
#' matrix that contains the analysis components for the corresponding block.}
#' @return \item{a}{A list of \eqn{J} elements. Each element of \eqn{a} is a
#' matrix that contains the outer weight vectors for each block.}
#' @return \item{astar}{A list of \eqn{J} elements. Each element of astar is a
#' matrix defined as Y[[j]][, h] = blocks[[j]]\%*\%astar[[j]][, h]}
#' @return \item{crit}{A vector of integer that contains for each component the
#' values of the analysis criteria across iterations.}
#' @return \item{AVE}{A list of numerical values giving the indicators of model
#' quality based on the Average Variance Explained (AVE): AVE(for each block),
#' AVE(outer model), AVE(inner model).}
#' @references Soon
#' @title Nonnegative Generalized Canonical Correlation Analysis (NGCCA)
#' @examples
#'
#'@export ngcca

ngcca <- function(blocks, connection = 1 - diag(length(blocks)), tau = 1,
                  ncomp = rep(1, length(blocks)), scheme = "centroid",
                  init = "svd", bias = TRUE, tol = .Machine$double.eps,
                  verbose = FALSE,   quiet = FALSE, na.rm = TRUE){

  # ndefl number of deflation per block
  ndefl <- ncomp - 1
  N <- max(ndefl)
  J <- length(blocks)
  pjs <- vapply(blocks, NCOL, numeric(1L))
  nb_ind <- NROW(blocks[[1]])
  AVE_X = list()
  AVE_outer <- rep(NA,max(ncomp))

  Y <- vector(mode = "list", length = J)
  a <- astar <- P <- vector(mode = "list", length = J)
  crit <- list()
  AVE_inner <- rep(NA,max(ncomp))

  for (b in seq_len(J))  {
    a[[b]] <- astar[[b]] <- matrix(NA, pjs[[b]], N + 1)
    Y[[b]] <- matrix(NA, nb_ind, N + 1)
  }

  ###################################################

  if (mode(scheme) != "function") {
    if (verbose) cat("Computation of the NGCCA block components based on the",
                     scheme, "scheme \n")
  }
  if (mode(scheme) == "function" & verbose) {
    cat("Computation of the NGCCA block components based on the g scheme \n")
  }

  ####################################
  # ngcca with 1 component per block #
  ####################################
  ngcca.result <- ngccak(blocks, connection,
                         scheme = scheme, init = init, bias = bias,
                         tol = tol, verbose = verbose, quiet = quiet,
                         na.rm = na.rm)

  for (b in seq_len(J)) {
    Y[[b]][, 1] <- ngcca.result$Y[, b, drop = FALSE]
    a[[b]][, 1] <- ngcca.result$a[[b]]
  }
  astar                      <- a
  AVE_inner[1]               <- ngcca.result$AVE_inner
  crit[[1]]                  <- ngcca.result$crit

  ##############################################
  #               If any ncomp > 1             #
  #      Determination of SGCCA components     #
  ##############################################
  if (N > 0) {
    R <- blocks
    for (b in seq_len(J)) {
      P[[b]] <- matrix(NA, pjs[[b]], N)
    }

    for (n in 2:(N + 1)) {
      if (verbose) message("Computation of the NGCCA block components #", n,
                           " is under progress... \n")

      # Apply deflation
      defla.result <- defl.select(ngcca.result$Y, R, ndefl, n - 1, J, na.rm = na.rm)
      R <- defla.result$resdefl

      ngcca.result <- ngccak(R, connection,
                             scheme = scheme, init = init, bias = bias,
                             tol = tol, verbose = verbose, quiet = quiet,
                             na.rm = na.rm)

      AVE_inner[n] <- ngcca.result$AVE_inner
      crit[[n]] <- ngcca.result$crit


      for (b in seq_len(J))  {
        Y[[b]][, n] <- ngcca.result$Y[, b]
        a[[b]][, n] <- ngcca.result$a[[b]]
        P[[b]][, n - 1] <- defla.result$pdefl[[b]]
        astar[[b]][, n] <- ngcca.result$a[[b]] - astar[[b]][, (1:(n - 1)), drop = F] %*% drop(crossprod(a[[b]][, n], P[[b]][, 1:(n - 1), drop = FALSE]))
      }
    }
  }

  for (b in seq_len(J)) {
    rownames(a[[b]]) = rownames(astar[[b]]) = colnames(blocks[[b]])
    rownames(Y[[b]]) = rownames(blocks[[b]])
    colnames(Y[[b]]) = paste0("comp", 1:max(ncomp))

    #Average Variance Explained (AVE) per block
    AVE_X[[b]] =  apply(cor(blocks[[b]], Y[[b]], use = "pairwise.complete.obs")^2, 2,
                        mean, na.rm = TRUE)
  }

  #AVE outer
  outer = matrix(unlist(AVE_X), nrow = max(ncomp))
  AVE_outer <- as.numeric((outer %*% pjs)/sum(pjs))

  Y = shave(Y, ncomp)
  AVE_X = shave(AVE_X, ncomp)

  AVE <- list(AVE_X = AVE_X, AVE_outer = AVE_outer, AVE_inner = AVE_inner)

  if (N == 0) crit = unlist(crit)

  out <- list(Y = shave(Y, ncomp),
              a = shave(a, ncomp),
              astar = shave(astar, ncomp),
              crit = crit,
              AVE = AVE)

  class(out) <- "ngcca"
  return(out)

}
