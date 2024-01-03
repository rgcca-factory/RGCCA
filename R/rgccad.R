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
#' @inheritParams rgcca
#' @param na.rm If TRUE, runs rgcca only on available data.
#' @param disjunction If TRUE, the response block is a one-hot encoded
#' qualitative block.
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
#' @importFrom graphics text

rgccad <- function(blocks, connection = 1 - diag(length(blocks)),
                   tau = rep(1, length(blocks)),
                   ncomp = rep(1, length(blocks)), scheme = "centroid",
                   init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,
                   na.rm = TRUE, superblock = FALSE, 
                   response = NULL, disjunction = NULL,
                   n_iter_max = 1000, comp_orth = TRUE,
                   groups = NULL, supergroup = FALSE) {
  if (verbose) {
    scheme_str <- ifelse(is(scheme, "function"), "user-defined", "scheme")
    cat(
      "Computation of the RGCCA block components based on the",
      scheme_str, "scheme \n"
    )
    tau_str <- ifelse(
      is.numeric(tau),
      "Shrinkage intensity parameters are chosen manually \n",
      "Optimal shrinkage intensity parameters are estimated \n"
    )
    cat(tau_str)
  }
  
  ##### Initialization #####
  # ndefl number of deflation per block
  ndefl <- ncomp - 1
  N <- max(ndefl)
  J <- length(blocks)
  pjs <- vapply(blocks, NCOL, FUN.VALUE = 1L)
  nb_ind <- vapply(blocks, NROW, FUN.VALUE = 1L)
  
  crit <- list()
  R <- blocks
  
  a <- lapply(seq(J), function(b) c())
  Y <- lapply(seq(J), function(b) c())
  if (!is.null(groups)) {
    L <- lapply(seq(J), function(b) c())
  }
  
  
  if ((superblock || supergroup) && comp_orth) {
    P <- c()
  } else {
    P <- lapply(seq(J), function(b) c())
  }
  
  # Whether primal or dual
  primal_dual <- rep("primal", J)
  if (is.null(groups)){
    primal_dual[which(nb_ind < pjs)] <- "dual"
  }
  
  # Save computed shrinkage parameter in a new variable
  computed_tau <- tau
  if (is.vector(tau)) {
    computed_tau <- matrix(
      rep(tau, N + 1),
      nrow = N + 1, J, byrow = TRUE
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
                          verbose = verbose, na.rm = na.rm, 
                          n_iter_max = n_iter_max, groups = groups
    )
    
    # Store tau, crit
    computed_tau[n, ] <- gcca_result$tau
    crit[[n]] <- gcca_result$crit
    
    # Store Y, a, factors and weights, 
    # as well as loadings L in the multi-group case
    a <- lapply(seq(J), function(b) cbind(a[[b]], gcca_result$a[[b]]))
    if (!is.null(groups)) {
      Y <- lapply(seq(J), function(b) cbind(Y[[b]], gcca_result$Y[[b]]))
      L <- lapply(seq(J), function(b) cbind(L[[b]], gcca_result$L[, b]))
    } else{
      Y <- lapply(seq(J), function(b) cbind(Y[[b]], gcca_result$Y[, b]))
    }
    
    # Deflation procedure
    if (n == N + 1) break
    defl_result <- deflate(gcca_result$a, gcca_result$Y, R, P, ndefl, n,
                           superblock, supergroup, comp_orth, response, na.rm)
    R <- defl_result$R
    P <- defl_result$P
  }
  
  # If there is a superblock and weight vectors are orthogonal, it is possible
  # to have non meaningful weights associated to blocks that have been set to
  # zero by the deflation
  if ((superblock || supergroup) && !comp_orth) {
    a <- lapply(a, function(x) {
      if (ncol(x) > nrow(x)) {
        x[, seq(nrow(x) + 1, ncol(x))] <- 0
      }
      return(x)
    })
  }
  
  ##### Generation of the output #####
  if (N == 0) {
    crit <- unlist(crit)
    computed_tau <- as.numeric(computed_tau)
  } else {
    computed_tau <- apply(computed_tau, 2, as.numeric)
  }
  
  astar <- compute_astar(a, P, superblock, supergroup, comp_orth, N)
  
  if (!is.null(groups)) {
    out <- list(
      Y = Y,
      a = a,
      L = L,
      astar = astar,
      tau = computed_tau,
      crit = crit, primal_dual = primal_dual)
  } else {
    out <- list(
      Y = Y,
      a = a,
      astar = astar,
      tau = computed_tau,
      crit = crit, primal_dual = primal_dual)
  }
  
  class(out) <- "rgccad"
  return(out)
}
