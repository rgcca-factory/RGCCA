#' Regularized Generalized Canonical Correlation Analysis (RGCCA) is a
#' generalization of regularized canonical correlation analysis to three or more
#' sets of variables.
#' @details
#' Given \eqn{J} matrices
#' \eqn{\mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}}{X1, X2, ..., XJ} that
#' represent \eqn{J} sets of variables observed on the
#' same set of \eqn{n} individuals. The matrices
#' \eqn{\mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}}{X1, X2, ..., XJ}
#' must have the same number of rows, but may (and usually will) have different
#' numbers of columns. The aim of RGCCA is to study  the relationships between
#' these \eqn{J} blocks of variables. It constitutes a general framework for
#' many multi-block data analysis methods. It combines the power of multi-block
#' data analysis methods (maximization of well identified criteria) and the
#' flexibility of PLS path modeling (the researcher decides which blocks are
#' connected and which are not). Hence, the use of RGCCA requires the
#' construction (user specified) of a design matrix,
#' (\eqn{\mathbf{connection}}{connection}), that characterize the connections
#' between blocks. Elements of the (symmetric) design matrix
#' \eqn{\mathbf{connection} = (c_{jk})}{connection = (c_jk)} is positive;
#' but usually equal to 1
#' if block \eqn{j} and block \eqn{k} are connected, and 0 otherwise. The
#' objective is to find a stationnary point related to the RGCCA optimization
#' problem. The function rgccad() implements a globally convergent algorithm
#' (i.e. monotone convergence that hits at convergence a stationary point).
#' Moreover, depending on the dimensionality of each block
#' \eqn{\mathbf{X}_j}{Xj},
#' \eqn{j = 1, \ldots, J}{j = 1, ..., J}, the primal (when \eqn{n > p_j})
#' algorithm or
#' the dual (when \eqn{n < p_j}) algorithm is used (see Tenenhaus et al. 2015).
#' Moreover, by deflation strategy, rgccad() allows to compute several RGCCA
#' block components (specified by ncomp) for each block. Using deflation, within
#' each block, block components are guaranteed to be orthogonal. The so-called
#' symmetric deflation is considered in this implementation, i.e. each block is
#' deflated with respect to its own component. It should be noted that the
#' numbers of components per block can differ from one block to another.
#' The rgcca() function can handle missing values using a NIPALS type algorithm
#' (non-linear iterative partial least squares algorithm) as described in
#' (Tenenhaus et al, 2005).
#' @inheritParams select_analysis
#' @param verbose Logical value indicating if the progress of the
#' algorithm is reported while computing.
#' @param init Character string giving the type of initialization to use in
#' the  algorithm. It could be either by Singular Value Decompostion ("svd")
#' or by random initialisation ("random") (default: "svd").
#' @param bias A logical value for biaised (\eqn{1/n}) or unbiaised
#' (\eqn{1/(n-1)}) estimator of the var/cov (default: bias = TRUE).
#' @param tol The stopping value for the convergence of the algorithm.
#' @param na.rm If TRUE, runs rgcca only on available data.
#' @param superblock TRUE if a superblock is added, FALSE otherwise (deflation
#' strategy must be adapted when a superblock is used).
#' @return \item{Y}{A list of \eqn{J} elements. Each element of the list is a
#' matrix that contains the RGCCA block components for the corresponding block.}
#' @return \item{a}{A list of \eqn{J} elements. Each element of the list \eqn{a}
#' is a matrix of block weight vectors for the corresponding block.}
#' @return \item{astar}{A list of \eqn{J} elements. Each column of astar[[j]] is
#' a vector such that Y[[j]][, h] = blocks[[j]] \%*\% astar[[j]][, h].}
#' @return \item{tau}{Either a \eqn{1 \times J}{1 x J} vector or a
#' \eqn{\mathrm{max}(ncomp) \times J}{max(ncomp) x J} matrix containing the
#' values of the regularization parameters. tau varies from 0
#' (maximizing the correlation) to 1 (maximizing the covariance).
#' If tau = "optimal" the regularization parameters are estimated for each
#' block and each dimension using the Schafer and Strimmer (2005) analytical
#' formula. If tau is a \eqn{1 \times J}{1 x J} vector, tau[j] is identical
#' across the dimensions of block \eqn{\mathbf{X}_j}{Xj}. If tau is a matrix,
#' tau[k, j] is associated with \eqn{\mathbf{X}_{jk}}{Xjk} (\eqn{k}th residual
#' matrix for block \eqn{j}). tau can be also estimated using
#' \link{rgcca_permutation}.}
#' @return \item{crit}{A list of vector of length max(ncomp). Each vector of
#' the list is related to one specific deflation stage and reports the values
#' of the criterion for this stage across iterations.}
#' @return \item{primal_dual}{A \eqn{1 \times J}{1 x J} vector that contains the
#' formulation ("primal" or "dual") applied to each of the \eqn{J} blocks within
#' the RGCCA alogrithm.}
#' @return \item{AVE}{A list of numerical values giving the goodness of fit
#' the model based on the Average Variance Explained (AVE): AVE(for each block),
#' AVE(outer model), AVE(inner model).}
#' @references Tenenhaus M., Tenenhaus A. and Groenen P. J. (2017). Regularized
#' generalized canonical correlation analysis: a framework for sequential
#' multiblock component methods. Psychometrika, 82(3), 737-777.
#' @references Tenenhaus A., Philippe C. and Frouin, V. (2015). Kernel
#' generalized canonical correlation analysis. Computational Statistics and
#' Data Analysis, 90, 114-131.
#' @references Tenenhaus A. and Tenenhaus M., (2011). Regularized Generalized
#' Canonical Correlation Analysis, Psychometrika, Vol. 76, Nr 2, pp 257-284.
#' @references Schafer J. and Strimmer K. (2005). A shrinkage approach to
#' large-scale covariance matrix estimation and implications for functional
#' genomics. Statist. Appl. Genet. Mol. Biol. 4:32.
#' @title Regularized Generalized Canonical Correlation Analysis (RGCCA)
#' @examples
#' #############
#' # Example 1 #
#' #############
#' data(Russett)
#' X_agric <- as.matrix(Russett[, c("gini", "farm", "rent")])
#' X_ind <- as.matrix(Russett[, c("gnpr", "labo")])
#' X_polit <- as.matrix(Russett[, c("demostab", "dictator")])
#' blocks <- list(X_agric, X_ind, X_polit)
#' # Define the design matrix (output = connection)
#' connection <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
#' fit.rgcca <- rgccad(blocks, connection,
#'   tau = c(1, 1, 1),
#'   scheme = "factorial"
#' )
#' lab <- as.vector(apply(Russett[, 9:11], 1, which.max))
#' plot(fit.rgcca$Y[[1]], fit.rgcca$Y[[2]],
#'   col = "white",
#'   xlab = "Y1 (Agric. inequality)", ylab = "Y2 (Industrial Development)"
#' )
#' text(fit.rgcca$Y[[1]], fit.rgcca$Y[[2]], rownames(Russett),
#'   col = lab, cex = .7
#' )
#'
#' ############################################
#' # Example 2: RGCCA and mutliple components #
#' ############################################
#' ############################
#' # plot(y1, y2) for (RGCCA) #
#' ############################
#' fit.rgcca <- rgccad(blocks, connection,
#'   tau = rep(1, 3), ncomp = c(2, 2, 1),
#'   scheme = "factorial", verbose = TRUE
#' )
#' layout(t(1:2))
#' plot(fit.rgcca$Y[[1]][, 1], fit.rgcca$Y[[2]][, 1],
#'   col = "white",
#'   xlab = "Y1 (Agric. inequality)", ylab = "Y2 (Industrial Development)",
#'   main = "Factorial plan of RGCCA"
#' )
#' text(fit.rgcca$Y[[1]][, 1], fit.rgcca$Y[[2]][, 1], rownames(Russett),
#'   col = lab, cex = .6
#' )
#' plot(fit.rgcca$Y[[1]][, 1], fit.rgcca$Y[[1]][, 2],
#'   col = "white",
#'   xlab = "Y1 (Agric. inequality)", ylab = "Y2 (Agric. inequality)",
#'   main = "Factorial plan of RGCCA"
#' )
#' text(fit.rgcca$Y[[1]][, 1], fit.rgcca$Y[[1]][, 2], rownames(Russett),
#'   col = lab, cex = .6
#' )
#'
#' ######################################
#' # example 3: RGCCA and leave one out #
#' ######################################
#' Ytest <- matrix(0, 47, 3)
#' fit.rgcca <- rgccad(blocks, connection,
#'   tau = rep(1, 3), ncomp = rep(1, 3),
#'   scheme = "factorial", verbose = TRUE
#' )
#' for (i in 1:nrow(Russett)) {
#'   B <- lapply(blocks, function(x) x[-i, ])
#'   B <- lapply(B, scale)
#'
#'   resB <- rgccad(B, connection,
#'     tau = rep(1, 3), scheme = "factorial",
#'     verbose = FALSE
#'   )
#'   # look for potential conflicting sign among components within
#'   # the loo loop.
#'   for (k in 1:length(B)) {
#'     if (cor(fit.rgcca$a[[k]], resB$a[[k]]) >= 0) {
#'       resB$a[[k]] <- resB$a[[k]]
#'     } else {
#'       resB$a[[k]] <- -resB$a[[k]]
#'     }
#'   }
#'   Btest <- lapply(blocks, function(x) x[i, ])
#'   Btest[[1]] <- (Btest[[1]] - attr(B[[1]], "scaled:center")) /
#'     (attr(B[[1]], "scaled:scale"))
#'   Btest[[2]] <- (Btest[[2]] - attr(B[[2]], "scaled:center")) /
#'     (attr(B[[2]], "scaled:scale"))
#'   Btest[[3]] <- (Btest[[3]] - attr(B[[3]], "scaled:center")) /
#'     (attr(B[[3]], "scaled:scale"))
#'   Ytest[i, 1] <- Btest[[1]] %*% resB$a[[1]]
#'   Ytest[i, 2] <- Btest[[2]] %*% resB$a[[2]]
#'   Ytest[i, 3] <- Btest[[3]] %*% resB$a[[3]]
#' }
#' lab <- apply(Russett[, 9:11], 1, which.max)
#' plot(fit.rgcca$Y[[1]], fit.rgcca$Y[[2]],
#'   col = "white",
#'   xlab = "Y1 (Agric. inequality)", ylab = "Y2 (Ind. Development)"
#' )
#' text(fit.rgcca$Y[[1]], fit.rgcca$Y[[2]], rownames(Russett), col = lab)
#' text(Ytest[, 1], Ytest[, 2], substr(rownames(Russett), 1, 1), col = lab)
#' @export rgccad
#' @importFrom graphics abline axis close.screen grid legend lines par points
#' rect screen segments split.screen text
#' @importFrom stats binomial glm lm predict sd var weighted.mean
#' @importFrom utils read.table write.table
#' @importFrom stats as.formula qt

rgccad <- function(blocks, connection = 1 - diag(length(blocks)),
                   tau = rep(1, length(blocks)),
                   ncomp = rep(1, length(blocks)), scheme = "centroid",
                   init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,
                   na.rm = TRUE, quiet = FALSE, superblock = FALSE) {
  update_col_n <- function(x, y, n) {
    x[, n] <- y
    return(x)
  }

  if (verbose) {
    if (mode(scheme) != "function") {
      cat(
        "Computation of the RGCCA block components based on the",
        scheme, "scheme \n"
      )
    } else {
      cat("Computation of the RGCCA block components based on the g scheme \n")
    }

    if (!is.numeric(tau)) {
      cat("Optimal Shrinkage intensity parameters are estimated \n")
    } else {
      cat("Shrinkage intensity parameters are chosen manually \n")
    }
  }

  ##### Initialization #####
  # ndefl number of deflation per block
  ndefl <- ncomp - 1
  N <- max(ndefl)
  J <- length(blocks)
  pjs <- sapply(blocks, NCOL)
  nb_ind <- NROW(blocks[[1]])
  AVE_inner <- rep(NA, max(ncomp))

  crit <- list()
  R <- blocks

  a <- lapply(seq(J), function(b) matrix(NA, pjs[[b]], N + 1))
  Y <- lapply(seq(J), function(b) matrix(NA, nb_ind, N + 1))

  if (superblock) {
    astar <- matrix(NA, pjs[J], N + 1)
    P <- matrix(NA, pjs[J], N)
  } else {
    astar <- a
    P <- lapply(seq(J), function(b) matrix(NA, pjs[[b]], N))
  }

  # Whether primal or dual
  primal_dual <- rep("primal", J)
  primal_dual[which(nb_ind < pjs)] <- "dual"

  # Save computed shrinkage parameter in a new variable
  computed_tau <- tau
  if (is.vector(tau)) {
    computed_tau <- matrix(
      rep(tau, N + 1),
      nrow = N + 1, J, byrow = T
    )
  }

  ##### Computation of RGCCA components #####
  for (n in seq(N + 1)) {
    if (verbose) {
      cat(paste0(
        "Computation of the RGCCA block components #", n,
        " is under progress...\n"
      ))
    }
    gcca_result <- rgccak(R, connection,
      tau = computed_tau[n, ], scheme = scheme,
      init = init, bias = bias, tol = tol,
      verbose = verbose, na.rm = na.rm
    )

    # Store tau, AVE_inner, crit
    computed_tau[n, ] <- gcca_result$tau
    AVE_inner[n] <- gcca_result$AVE_inner
    crit[[n]] <- gcca_result$crit

    # Store Y, a, factors and weights
    Y <- lapply(seq(J), function(b) update_col_n(Y[[b]], gcca_result$Y[, b], n))
    a <- lapply(seq(J), function(b) update_col_n(a[[b]], gcca_result$a[[b]], n))

    # Compute astar
    if (superblock) {
      if (n == 1) {
        astar[, 1] <- a[[J]][, 1, drop = FALSE]
      } else {
        astar[, n] <- gcca_result$a[[J]] -
          astar[, seq(n - 1), drop = F] %*%
          drop(t(a[[J]][, n]) %*% P[, seq(n - 1), drop = F])
      }
    } else {
      if (n == 1) {
        astar <- a
      } else {
        astar <- lapply(seq(J), function(b) {
          update_col_n(
            astar[[b]],
            gcca_result$a[[b]] - astar[[b]][, seq(n - 1), drop = F] %*%
              drop(t(a[[b]][, n]) %*% P[[b]][, seq(n - 1), drop = F]),
            n
          )
        })
      }
    }

    # Deflation procedure
    if (n == N + 1) break
    if (superblock) {
      defl_result <- deflation(R[[J]], gcca_result$Y[, J])
      P[, n] <- defl_result$p
      cumsum_pjs <- cumsum(pjs)[seq_len(J - 1)]
      inf_pjs <- c(0, cumsum_pjs[seq_len(J - 2)]) + 1
      R <- lapply(seq(J - 1), function(b) {
        x <- defl_result$R[, inf_pjs[b]:cumsum_pjs[b], drop = FALSE]
        colnames(x) <- colnames(defl_result$R)[inf_pjs[b]:cumsum_pjs[b]]
        return(x)
      })
      R[[J]] <- defl_result$R
    } else {
      defl_result <- defl_select(gcca_result$Y, R,
        ndefl, n - 1, J,
        na.rm = na.rm
      )
      R <- defl_result$resdefl
      P <- lapply(seq(J), function(b) {
        update_col_n(
          P[[b]], defl_result$pdefl[[b]], n
        )
      })
    }
  }

  ##### Generation of the output #####
  AVE_X <- lapply(seq(J), function(b) {
    apply(
      cor(blocks[[b]], Y[[b]], use = "pairwise.complete.obs")^2, 2, mean
    )
  })

  outer <- matrix(unlist(AVE_X), nrow = max(ncomp))
  AVE_outer <- as.vector((outer %*% pjs) / sum(pjs))
  AVE_X <- shave(AVE_X, ncomp)
  AVE <- list(AVE_X = AVE_X, AVE_outer = AVE_outer, AVE_inner = AVE_inner)

  if (N == 0) {
    crit <- unlist(crit)
    computed_tau <- as.vector(computed_tau)
  }

  out <- list(
    Y = Y,
    a = a,
    astar = astar,
    tau = computed_tau,
    crit = crit, primal_dual = primal_dual,
    AVE = AVE
  )

  class(out) <- "rgccad"
  return(out)
}
