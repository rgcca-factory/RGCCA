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
#' @param scale Logical value indicating if blocks are standardized.
#' @param scale_block Logical value indicating if each block is divided by
#' the square root of its number of variables.
#' @param verbose Logical value indicating if the progress of the
#' algorithm is reported while computing.
#' @param quiet Logical value indicating if warning messages are reported.
#' @param init Character string giving the type of initialization to use in
#' the  algorithm. It could be either by Singular Value Decompostion ("svd")
#' or by random initialisation ("random") (default: "svd").
#' @param bias A logical value for biaised (\eqn{1/n}) or unbiaised
#' (\eqn{1/(n-1)}) estimator of the var/cov (default: bias = TRUE).
#' @param tol The stopping value for the convergence of the algorithm.
#' @param prescaling Logical value indicating if the scaling has been done
#' outside of the function.
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
#' @return \item{call}{Call of the function}
#' @references Tenenhaus M., Tenenhaus A. and Groenen P. J. (2017). Regularized
#' generalized canonical correlation analysis: a framework for sequential
#' multiblock component methods. Psychometrika, 82(3), 737-777.
#' @references Tenenhaus A., Philippe C and Frouin, V. (2015). Kernel
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
#' fit.rgcca = rgccad(blocks, connection, tau = c(1, 1, 1), scheme = "factorial",
#' scale = TRUE)
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
#'                      scheme = "factorial", verbose = TRUE, scale = TRUE,
#'                      scale_block = FALSE)
#' for (i in 1:nrow(Russett)){
#'  B = lapply(blocks, function(x) x[-i, ])
#'  B = lapply(B, scale)
#'
#'  resB = rgccad(B, connection, tau = rep(1, 3), scheme = "factorial",
#'  scale = TRUE, scale_block = FALSE, verbose = FALSE)
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
#' @importFrom grDevices dev.off png rainbow
#' @importFrom graphics abline axis close.screen grid legend lines par points
#' rect screen segments split.screen text
#' @importFrom stats binomial glm lm predict sd var weighted.mean
#' @importFrom utils read.table write.table
#' @importFrom stats as.formula qt
#' @importFrom grDevices graphics.off

rgccad=function (blocks, connection = 1 - diag(length(blocks)), tau = rep(1, length(blocks)),
                 ncomp = rep(1, length(blocks)), scheme = "centroid", scale = TRUE,
                 init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,
                 scale_block = TRUE, na.rm = TRUE,
                 prescaling = FALSE, quiet = FALSE)
{

  shave.matlist <- function(mat_list, nb_cols)
    mapply(function(m, nbcomp) m[, 1:nbcomp, drop = FALSE],
           mat_list, nb_cols,
           SIMPLIFY = FALSE
           )

  shave.veclist <- function(vec_list, nb_elts)
    mapply(function(m, nbcomp) m[1:nbcomp], vec_list, nb_elts, SIMPLIFY = FALSE)

  A0 = blocks
  call=list(blocks = blocks, connection = connection,  ncomp = ncomp, scheme = scheme, scale = scale,
            init = init, bias = bias, tol = tol, verbose = verbose,
            scale_block = scale_block, na.rm = na.rm)

  if (any(ncomp < 1)) {stop_rgcca("Compute at least one component per block!")}
  pjs <- sapply(blocks, NCOL)
  nb_row <- NROW(blocks[[1]])

  if (any(ncomp - pjs > 0))
    stop_rgcca("For each block, choose a number of components smaller than the
               number of variables!")

  if (mode(scheme) != "function") {
    if ((scheme != "horst") & (scheme != "factorial") & (scheme != "centroid")){
      stop_rgcca("Choose one of the three following schemes: horst, centroid,
                 factorial or design the g function")
    }
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

  if(!prescaling)
      blocks <- scaling(blocks, scale = scale, bias = bias, scale_block = scale_block)

  # Superblock option
  if(!is.matrix(connection)&& connection == "superblock")
  {
    #Construction of the superblock
    blocks = c(blocks, list(do.call(cbind, blocks)))
    #Construction of the corresponding design matrix
    connection = matrix(0, length(blocks), length(blocks))
    connection[length(blocks), 1:(length(blocks)-1)]=1
    connection = connection+t(connection)
    # Shrinkage parameters

    if(is.null(tau)){
      message("the shrinkage parameters have been
              automatically set to 1 for all blocks (incl. superblock)")
      tau = rep(1, NCOL(connection))
    }
    if(length(tau) == NCOL(connection)-1){
      message("the shrinkage parameter for the superblock has been
              automatically set to 1")
      tau=c(tau,1)
    }
    if(!((length(tau) == NCOL(connection)-1) | (length(tau) == NCOL(connection))))
      stop_rgcca("the length of the vector of shinkage parameters is not
                 appropriate.")

    # number of components per block
    if(is.null(ncomp)){
      message("the number of components per block has been
              automatically set to 1 for all blocks/superblock)")
      ncomp = rep(1, NCOL(connection))
    }
    if(length(ncomp) == NCOL(connection)-1){
      message("the number of global components has been
              automatically set to max(ncomp)")
      ncomp =c(ncomp,max(ncomp))
    }
    if(!((length(ncomp) == NCOL(connection)-1) | (length(tau) == NCOL(connection))))
      stop_rgcca("the ncomp argument has been filled inappropriately.")

    # number of variables per block
    pjs = c(pjs, sum(pjs))
  }

  AVE_X = list()
  AVE_outer <- vector()
  ndefl <- ncomp - 1
  N <- max(ndefl)
  nb_ind <- NROW(blocks[[1]])
  J <- length(blocks)

  # Whether primal or dual
  primal_dual = rep("primal", J)
  primal_dual[which(nb_row < pjs)] = "dual"

  # One component per block
  if(N == 0){
    result <- rgccak(blocks, connection, tau = tau, scheme = scheme, init = init,
                     bias = bias, tol = tol, verbose = verbose,
                     na.rm=na.rm,
                     scale_block=scale_block, scale=scale)

    Y <- NULL
    for (b in 1:J) Y[[b]] <- result$Y[, b, drop = FALSE]
    for (j in 1:J)
      AVE_X[[j]] = mean(cor(blocks[[j]], Y[[j]],use="pairwise.complete.obs")^2,
                        na.rm=TRUE)

    AVE_outer <- sum(pjs * unlist(AVE_X))/sum(pjs)
    AVE <- list(AVE_X = AVE_X, AVE_outer = AVE_outer,
                AVE_inner = result$AVE_inner)
    a <- lapply(result$a, cbind)

     for (b in 1:J) {
      rownames(a[[b]]) = colnames(blocks[[b]])
      rownames(Y[[b]]) = rownames(blocks[[b]])
      colnames(Y[[b]]) = "comp1"
    }

    tau=result$tau

    if(is.vector(tau))
      names(tau) = names(blocks)

    out <- list(Y = Y, a = a, astar = a, connection = connection,  scheme = scheme,
                ncomp = ncomp, crit = result$crit,
                primal_dual = primal_dual,
                AVE = AVE, A = A0, tau = tau, call = call)

    class(out) <- "rgccad"

    return(out)
  }

  Y <- NULL
  crit = list()
  AVE_inner <- rep(NA, max(ncomp))
  R <- blocks
  P <- a <- astar <- NULL
  if (is.numeric(tau))
  {
      tau_mat = tau
      if(is.vector(tau_mat)) names(tau_mat) = names(blocks)
      if(is.matrix(tau_mat))colnames(tau_mat) = names(blocks)
  }
  else
  {
      tau_mat = matrix(NA, max(ncomp), J)
      colnames(tau_mat) = names(blocks)
  }

  for (b in 1:J) P[[b]] <- a[[b]] <- astar[[b]] <- matrix(NA, pjs[[b]], N + 1)
  for (b in 1:J) Y[[b]] <- matrix(NA, nb_ind, N + 1)
  for (n in 1:N) {
     if (verbose)
      cat(paste0("Computation of the RGCCA block components #", n, " is under
                 progress...\n"))
    if (is.vector(tau))
      rgcca.result <- rgccak(R, connection, tau = tau, scheme = scheme,init = init,
                             bias = bias, tol = tol, verbose = verbose,
                             na.rm = na.rm)
    else rgcca.result <- rgccak(R, connection, tau = tau[n, ], scheme = scheme,
                                init = init, bias = bias, tol = tol,
                                verbose = verbose, na.rm = na.rm)

    if (!is.numeric(tau)) tau_mat[n, ] = rgcca.result$tau

    AVE_inner[n] <- rgcca.result$AVE_inner
    crit[[n]] <- rgcca.result$crit

    # deflation
    for (b in 1:J) Y[[b]][, n] <- rgcca.result$Y[, b]
    defla.result <- defl.select(rgcca.result$Y, R, ndefl , n, nbloc = J)
    R <- defla.result$resdefl
    for (b in 1:J) P[[b]][, n] <- defla.result$pdefl[[b]]
    for (b in 1:J) a[[b]][, n] <- rgcca.result$a[[b]]
    if (n == 1)
    {
      for (b in 1:J) astar[[b]][, n] <- rgcca.result$a[[b]]
    }
    else {
      for (b in 1:J)
      {
        astar[[b]][, n] <- rgcca.result$a[[b]]-astar[[b]][, (1:n-1), drop = F]%*%
          drop(t(a[[b]][, n])%*%P[[b]][, 1:(n-1), drop = F])
      }
    }
  }
  if (verbose)
    cat(paste0("Computation of the RGCCA block components #",
               N + 1, " is under progress ... \n"))
  if (is.vector(tau))
    rgcca.result <- rgccak(R, connection, tau = tau, scheme = scheme, init = init,
                           bias = bias, tol = tol, verbose = verbose)
  else rgcca.result <- rgccak(R, connection, tau = tau[N+1, ], scheme = scheme,
                              init = init, bias = bias, tol = tol,
                              verbose = verbose)

  crit[[N+1]] <- rgcca.result$crit
  if (!is.numeric(tau))
    tau_mat[N+1, ] = rgcca.result$tau
  AVE_inner[max(ncomp)] <- rgcca.result$AVE_inner
  for (b in 1:J) {
    Y[[b]][, N+1] <- rgcca.result$Y[, b]
    a[[b]][, N+1] <- rgcca.result$a[[b]]
    astar[[b]][, N+1] <- rgcca.result$a[[b]]-astar[[b]][, (1:N), drop = F]%*%
      drop(t(a[[b]][, (N+1)]) %*%P[[b]][, 1:N, drop = F])
    rownames(a[[b]]) = rownames(astar[[b]]) = colnames(blocks[[b]])
    rownames(Y[[b]]) = rownames(blocks[[b]])
    colnames(Y[[b]]) = paste0("comp", 1:max(ncomp))
  }

  for (j in 1:J)
      AVE_X[[j]] = apply(cor(blocks[[j]], Y[[j]], use="pairwise.complete.obs")^2,
                         2, mean, na.rm = TRUE)

  outer = matrix(unlist(AVE_X), nrow = max(ncomp))

  for (j in 1:max(ncomp))
    AVE_outer[j] <- sum(pjs * outer[j,],na.rm=na.rm)/sum(pjs)

  Y = shave.matlist(Y, ncomp)
  names(Y)=names(blocks)
  names(a)=names(blocks)
  AVE_X = shave.veclist(AVE_X, ncomp)

  AVE <- list(AVE_X = AVE_X, AVE_outer_model = AVE_outer,
              AVE_inner_model = AVE_inner)

  out <- list(Y = shave.matlist(Y, ncomp), a = shave.matlist(a,ncomp),
              astar = shave.matlist(astar, ncomp),
              tau = tau_mat,
              crit = crit, primal_dual = primal_dual,
              AVE = AVE, A = blocks, call = call)

   class(out) <- "rgccad"

  return(out)
}
