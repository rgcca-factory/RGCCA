#' Regularized Generalized Canonical Correlation Analysis (RGCCA) is a
#' generalization of regularized canonical correlation analysis to three or more
#' sets of variables.
#' @details
#' Given \eqn{J} matrices \eqn{\mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}}
#' that represent \eqn{J} sets of variables observed on the same set of \eqn{n}
#' individuals. The matrices \eqn{\mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}}
#' must have the same number of rows, but may (and usually will) have different
#' numbers of columns. The aim of RGCCA is to study  the relationships between
#' these \eqn{J} blocks of variables. It constitutes a general framework for
#' many multi-block data analysis methods. It combines the power of multi-block
#' data analysis methods (maximization of well identified criteria) and the
#' flexibility of PLS path modeling (the researcher decides which blocks are
#' connected and which are not). Hence, the use of RGCCA requires the
#' construction (user specified) of a design matrix, (\eqn{\mathbf{connection}}), that
#' characterize the connections between blocks. Elements of the (symmetric)
#' design matrix \eqn{\mathbf{connection} = (c_{jk})} is positive ; but usually equal to 1
#' if block \eqn{j} and block \eqn{k} are connected, and 0 otherwise. The
#' objective is to find a stationnary point related to the RGCCA optimization
#' problem. The function rgccad() implements a globally convergent algorithm
#' (i.e. monotone convergence that hits at convergence a stationary point).
#' Moreover, depending on the dimensionality of each block \eqn{\mathbf{X}_j},
#' \eqn{j = 1, \ldots, J}, the primal (when \eqn{n > p_j}) algorithm or
#' the dual (when \eqn{n < p_j}) algorithm is used (see Tenenhaus et al. 2015).
#' Moreover, by deflation strategy, rgccad() allow to compute several RGCCA
#' block components (specified by ncomp) for each block. Using deflation, within
#' each block, block components are guaranteed to be orthogonal. The so-called
#' symmetric deflation is considered in this implementation, i.e. each block is
#' deflated with respect to its own component. It should be noted that the
#' numbers of components per block can differ from one block to another.
#' The rgcca() function can handle missing values using a NIPALS type algorithm
#' (non-linear iterative partial least squares algorithm) as described in
#' (Tenenhaus et al, 2005).
#' @inheritParams select_analysis
#' @param blocks A list that contains the J blocks of variables X1, X2, ..., XJ.
#' Block Xj is a matrix of dimension n x p_j where n is the number of
#' observations and p_j the number of variables.
#' @param connection  A symmetric matrix (J*J) that describes the relationships between
#' blocks.
#' @param tau Either a 1*J vector or a max(ncomp)*J matrix containing
#' the values of the regularization parameters (default: tau = 1, for each
#' block and each dimension). The regularization parameters varies from 0
#' (maximizing the correlation) to 1 (maximizing the covariance). If
#' tau = "optimal" the regularization paramaters are estimated for each block
#' and each dimension using the Schafer and Strimmer (2005) analytical formula.
#' If tau is a 1*J vector, tau[j] is identical across the dimensions
#' of block Xj. If tau is a matrix, tau[k, j] is associated with
#' X_jk (kth residual matrix for block j). The regularization parameters can
#' also be estimated using \link{rgcca_permutation} or \link{rgcca_cv}.
#' @param verbose Logical value indicating if the progress of the
#' algorithm is reported while computing.
#' @param quiet Logical value indicating if warning messages are reported.
#' @param init Character string giving the type of initialization to use in
#' the  algorithm. It could be either by Singular Value Decompostion ("svd")
#' or by random initialisation ("random") (default: "svd").
#' @param bias A logical value for biaised (\eqn{1/n}) or unbiaised
#' (\eqn{1/(n-1)}) estimator of the var/cov (default: bias = TRUE).
#' @param tol The stopping value for the convergence of the algorithm.
#' @param na.rm If TRUE, runs rgcca only on available data.
#' @return \item{Y}{A list of \eqn{J} elements. Each element of the list is a
#' matrix that contains the RGCCA block components for the corresponding block.}
#' @return \item{a}{A list of \eqn{J} elements. Each element of the list \eqn{a}
#' is a matrix of block weight vectors for the corresponding block.}
#' @return \item{astar}{A list of \eqn{J} elements. Each element of astar is a
#' matrix defined as Y[[j]][, h] = blocks[[j]]\%*\%astar[[j]][, h].}
#' @return \item{tau}{Either a 1*J vector or a \eqn{\mathrm{max}(ncomp) \times J}
#' matrix containing the values of the regularization parameters. tau varies
#' from 0 (maximizing the correlation) to 1 (maximizing the covariance).
#' If tau = "optimal" the regularization paramaters are estimated for each
#' block and each dimension using the Schafer and Strimmer (2005) analytical
#' formula. If tau is a \eqn{1\times J} vector, tau[j] is identical across the
#' dimensions of block \eqn{\mathbf{X}_j}. If tau is a matrix, tau[k, j] is
#' associated with \eqn{\mathbf{X}_{jk}} (\eqn{k}th residual matrix for
#' block \eqn{j}). tau can be also estimated using \link{rgcca_permutation}.}
#' @return \item{crit}{A list of vector of length max(ncomp). Each vector of
#' the list is related to one specific deflation stage and reports the values
#' of the criterion for this stage across iterations.}
#' @return \item{primal_dual}{A \eqn{1 \times J} vector that contains the
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
#' X_agric =as.matrix(Russett[,c("gini","farm","rent")])
#' X_ind = as.matrix(Russett[,c("gnpr","labo")])
#' X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
#' blocks = list(X_agric, X_ind, X_polit)
#' #Define the design matrix (output = connection)
#' connection = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
#' fit.rgcca = rgccad(blocks, connection, tau = c(1, 1, 1), scheme = "factorial")
#' lab = as.vector(apply(Russett[, 9:11], 1, which.max))
#' plot(fit.rgcca$Y[[1]], fit.rgcca$Y[[2]], col = "white",
#'      xlab = "Y1 (Agric. inequality)", ylab = "Y2 (Industrial Development)")
#' text(fit.rgcca$Y[[1]], fit.rgcca$Y[[2]], rownames(Russett),
#'      col = lab, cex = .7)
#'
#' ############################################
#' # Example 2: RGCCA and mutliple components #
#' ############################################
#' ############################
#' # plot(y1, y2) for (RGCCA) #
#' ############################
#' fit.rgcca = rgccad(blocks, connection, tau = rep(1, 3), ncomp = c(2, 2, 1),
#'                      scheme = "factorial", verbose = TRUE)
#' layout(t(1:2))
#' plot(fit.rgcca$Y[[1]][, 1], fit.rgcca$Y[[2]][, 1], col = "white",
#' xlab = "Y1 (Agric. inequality)", ylab = "Y2 (Industrial Development)",
#' main = "Factorial plan of RGCCA")
#' text(fit.rgcca$Y[[1]][, 1], fit.rgcca$Y[[2]][, 1], rownames(Russett),
#' col = lab, cex = .6)
#' plot(fit.rgcca$Y[[1]][, 1], fit.rgcca$Y[[1]][, 2], col = "white",
#' xlab = "Y1 (Agric. inequality)", ylab = "Y2 (Agric. inequality)",
#' main = "Factorial plan of RGCCA")
#' text(fit.rgcca$Y[[1]][, 1], fit.rgcca$Y[[1]][, 2], rownames(Russett),
#' col = lab, cex = .6)
#'
#' ######################################
#' # example 3: RGCCA and leave one out #
#' ######################################
#' Ytest = matrix(0, 47, 3)
#' fit.rgcca = rgccad(blocks, connection, tau = rep(1, 3), ncomp = rep(1, 3),
#'                      scheme = "factorial", verbose = TRUE)
#' for (i in 1:nrow(Russett)){
#'  B = lapply(blocks, function(x) x[-i, ])
#'  B = lapply(B, scale)
#'
#'  resB = rgccad(B, connection, tau = rep(1, 3), scheme = "factorial", verbose = FALSE)
#'  #  look for potential conflicting sign among components within the loo loop.
#'  for (k in 1:length(B)){
#'    if (cor(fit.rgcca$a[[k]], resB$a[[k]]) >= 0)
#'      resB$a[[k]] = resB$a[[k]] else resB$a[[k]] = -resB$a[[k]]
#'  }
#'  Btest = lapply(blocks, function(x) x[i, ])
#'  Btest[[1]] = (Btest[[1]]-attr(B[[1]],"scaled:center"))/
#'                   (attr(B[[1]],"scaled:scale"))
#'  Btest[[2]] = (Btest[[2]]-attr(B[[2]],"scaled:center"))/
#'                   (attr(B[[2]],"scaled:scale"))
#'  Btest[[3]] = (Btest[[3]]-attr(B[[3]],"scaled:center"))/
#'                   (attr(B[[3]],"scaled:scale"))
#'  Ytest[i, 1] = Btest[[1]]%*%resB$a[[1]]
#'  Ytest[i, 2] = Btest[[2]]%*%resB$a[[2]]
#'  Ytest[i, 3] = Btest[[3]]%*%resB$a[[3]]
#' }
#' lab = apply(Russett[, 9:11], 1, which.max)
#' plot(fit.rgcca$Y[[1]], fit.rgcca$Y[[2]], col = "white",
#'      xlab = "Y1 (Agric. inequality)", ylab = "Y2 (Ind. Development)")
#' text(fit.rgcca$Y[[1]], fit.rgcca$Y[[2]], rownames(Russett), col = lab)
#' text(Ytest[, 1], Ytest[, 2], substr(rownames(Russett), 1, 1), col = lab)
#' @export rgccad
#' @param superblock TRUE if a superblock is added, FALSE ifelse (useful for deflation strategy when a superblock is used)
#' @importFrom grDevices dev.off png rainbow
#' @importFrom graphics abline axis close.screen grid legend lines par points
#' rect screen segments split.screen text
#' @importFrom stats binomial glm lm predict sd var weighted.mean
#' @importFrom utils read.table write.table
#' @importFrom stats as.formula qt
#' @importFrom grDevices graphics.off

rgccad = function(blocks, connection = 1 - diag(length(blocks)),
                  tau = rep(1, length(blocks)),
                  ncomp = rep(1, length(blocks)), scheme = "centroid",
                  init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,
                  na.rm = TRUE, quiet = FALSE, superblock = FALSE)
{

  if (mode(scheme) != "function") {
    if (verbose)
      cat("Computation of the RGCCA block components based on the",
          scheme, "scheme \n")
  }
  if (mode(scheme) == "function" & verbose)
    cat("Computation of the RGCCA block components based on the g scheme \n")

  if (!is.numeric(tau) & verbose) {
    cat("Optimal Shrinkage intensity paramaters are estimated \n")
  }
  else {
    if (is.numeric(tau) & verbose) {
      cat("Shrinkage intensity paramaters are chosen manually \n")
    }
  }

  # ndefl number of deflation per block
  ndefl <- ncomp - 1
  N <- max(ndefl)
  J <- length(blocks)
  pjs <- sapply(blocks,NCOL)
  nb_ind <- NROW(blocks[[1]])
  AVE_X <- list()
  AVE_outer <- rep(NA, max(ncomp))

  Y <- NULL
  P <- a <- astar <- list()
  crit <- list()
  AVE_inner <- rep(NA,max(ncomp))

  for (b in seq_len(J)){
    a[[b]] <- matrix(NA, pjs[[b]], N + 1)
    Y[[b]] <- matrix(NA, nb_ind, N + 1)
  }

  if(!superblock){
    for (b in seq_len(J)) astar[[b]] <- matrix(NA, pjs[b], N + 1)
  }else{
    astar <- matrix(NA, pjs[J], N + 1)
  }

  # Whether primal or dual
  primal_dual = rep("primal", J)
  primal_dual[which(nb_ind < pjs)] = "dual"

  # Save computed shrinkage parameter in a new variable
  computed_tau = tau
  if (is.vector(tau)) computed_tau = matrix(NA, nrow = N + 1, J)

  # First component block
  if (is.vector(tau))
    rgcca.result <- rgccak(blocks, connection, tau = tau, scheme = scheme,
                           init = init, bias = bias, tol = tol,
                           verbose = verbose, na.rm = na.rm)
  else
    rgcca.result <- rgccak(blocks, connection, tau = tau[1, ], scheme = scheme,
                           init = init, bias = bias, tol = tol,
                           verbose = verbose, na.rm = na.rm)
  computed_tau[1, ] = rgcca.result$tau

  for (b in seq_len(J)) Y[[b]][, 1] <- rgcca.result$Y[, b, drop = FALSE]
  for (b in seq_len(J)) a[[b]][, 1] <- rgcca.result$a[[b]]

  ifelse(!superblock,
         astar <- a,
         astar[, 1] <- a[[J]][, 1, drop = FALSE])

  AVE_inner[1] <- rgcca.result$AVE_inner
  crit[[1]] <- rgcca.result$crit

  if (N > 0) {
    R <- blocks

    if(!superblock){
      for (b in seq_len(J)) P[[b]] <- matrix(NA, pjs[b], N)
    }else{
      P <- matrix(NA, pjs[J], N)
    }

    for (n in 2:(N + 1)) {
      if (verbose)
        cat(paste0("Computation of the RGCCA block components #", n, " is under
                 progress...\n"))

      if(!superblock){
        defl.result <- defl.select(rgcca.result$Y, R,
                                   ndefl, n - 1, J,
                                   na.rm = na.rm)
        R <- defl.result$resdefl
        for (b in seq_len(J)) P[[b]][, n - 1] <- defl.result$pdefl[[b]]
      }else{
        defl.result <- deflation(R[[J]], rgcca.result$Y[, J])
        R[[J]] <- defl.result$R
        P[, n - 1] <- defl.result$p
        cumsum_pjs <- cumsum(pjs)[seq_len(J-1)]
        inf_pjs <- c(0, cumsum_pjs[seq_len(J-2)])+1
        for(j in seq_len(J-1)){
          R[[j]] = R[[J]][ , inf_pjs[j]:cumsum_pjs[j], drop=FALSE]
          rownames(R[[j]]) <- rownames(R[[j]])
          colnames(R[[j]]) <- colnames(R[[J]])[inf_pjs[j]:cumsum_pjs[j]]
        }
      }

      if (is.vector(tau))
        rgcca.result <- rgccak(R, connection, tau = tau, scheme = scheme,
                               init = init, bias = bias, tol = tol,
                               verbose = verbose, na.rm = na.rm)
      else
        rgcca.result <- rgccak(R, connection, tau = tau[n, ], scheme = scheme,
                               init = init, bias = bias, tol = tol,
                               verbose = verbose, na.rm = na.rm)

      computed_tau[n, ] <- rgcca.result$tau

      AVE_inner[n] <- rgcca.result$AVE_inner
      crit[[n]] <- rgcca.result$crit

      for (b in seq_len(J)) Y[[b]][, n] <- rgcca.result$Y[, b]
      for (b in seq_len(J)) a[[b]][, n] <- rgcca.result$a[[b]]

      if(!superblock){
        for (b in seq_len(J))
          astar[[b]][, n] <- rgcca.result$a[[b]] -
            astar[[b]][, (1:(n - 1)), drop = F] %*%
        drop( t(a[[b]][, n]) %*% P[[b]][, 1:(n - 1), drop = F])
      }else{
        astar[, n] <- rgcca.result$a[[J]] -
          astar[, (1:(n - 1)), drop = F] %*%
            drop(t(a[[J]][, n]) %*% P[, 1:(n - 1), drop = F])
      }
    }
  }

  for (j in seq_len(J)) AVE_X[[j]] = apply(
    cor(blocks[[j]], Y[[j]], use = "pairwise.complete.obs")^2, 2, mean)

  outer = matrix(unlist(AVE_X), nrow = max(ncomp))

  for (j in seq_len(max(ncomp)))
    AVE_outer[j] <- sum(pjs * outer[j,])/sum(pjs)

  AVE_X = shave(AVE_X, ncomp)

  AVE <- list(AVE_X = AVE_X, AVE_outer = AVE_outer, AVE_inner = AVE_inner)

  if (N == 0) {
    crit = unlist(crit)
    computed_tau = as.vector(computed_tau)
  }

  out <- list(Y = Y,
              a = a,
              astar = astar,
              tau = computed_tau,
              crit = crit, primal_dual = primal_dual,
              AVE = AVE)

   class(out) <- "rgccad"

  return(out)
}
