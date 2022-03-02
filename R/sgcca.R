#' SGCCA extends RGCCA to address the issue of variable selection. Specifically,
#' RGCCA is combined with an L1-penalty that gives rise to Sparse GCCA (SGCCA)
#' which is implemented in the function sgcca().
#' Given \eqn{J} matrices \eqn{X_1, X_2, ..., X_J}, that represent
#' \eqn{J} sets of variables observed on the same set of \eqn{n} individuals.
#' The matrices \eqn{X_1, X_2, ..., X_J} must have the same number of rows, but
#' may (and usually will) have different numbers of columns. Blocks are not
#' necessarily fully connected within the SGCCA framework. Hence the use of
#' SGCCA requires the construction (user specified) of a design matrix (\eqn{connection})
#' that characterizes the connections between blocks. Elements of the symmetric
#' design matrix \eqn{connection = (c_{jk})} are equal to 1 if block \eqn{j} and block
#' \eqn{k} are connected, and 0 otherwise. The SGCCA algorithm is very similar
#' to the RGCCA algorithm and keeps the same monotone convergence properties
#' (i.e. the bounded criteria to be maximized increases at each step of the
#' iterative procedure and hits at convergence a stationary point).
#' Moreover, using a deflation strategy, sgcca() enables the computation of
#' several SGCCA block components (specified by ncomp) for each block. Block
#' components for each block are guaranteed to be orthogonal when using this
#' deflation strategy. The so-called symmetric deflation is considered in this
#' implementation, i.e. each block is deflated with respect to its own
#' component. Moreover, we stress that the numbers of components per block
#' could differ from one block to another.
#' @inheritParams select_analysis
#' @inheritParams rgccad
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
#' require(gliomaData)
#' data(ge_cgh_locIGR)
#'
#' blocks <- ge_cgh_locIGR$multiblocks
#' Loc <- factor(ge_cgh_locIGR$y)
#' levels(Loc) <- colnames(ge_cgh_locIGR$multiblocks$y)
#' connection <-  matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
#' tau = c(1, 1, 0)
#'
#' # rgcca algorithm using the dual formulation for X1 and X2
#' # and the dual formulation for X3
#' blocks[[3]] = blocks[[3]][, -3]
#' result.rgcca = rgcca(blocks = blocks, connection = connection, tau = tau,
#'                      ncomp = c(2, 2, 1), scheme = "factorial",
#'                      verbose = TRUE, method = "rgcca")
#' # sgcca algorithm
#' result.sgcca = rgcca(blocks = blocks, connection = connection,
#'                      sparsity = c(.071,.2, 1), ncomp = c(2, 2, 1),
#'                      scheme = "centroid", verbose = TRUE, method = "sgcca")
#'
#' ############################
#' # plot(y1, y2) for (RGCCA) #
#' ############################
#' layout(t(1:2))
#' plot(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[2]][, 1],
#'      col = "white", xlab = "Y1 (GE)", ylab = "Y2 (CGH)",
#'      main = "Factorial plan of RGCCA")
#' text(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[2]][, 1],
#'      Loc, col = as.numeric(Loc), cex = .6)
#' plot(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[1]][, 2],
#'      col = "white", xlab = "Y1 (GE)",
#'      ylab = "Y2 (GE)", main = "Factorial plan of RGCCA")
#' text(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[1]][, 2],
#'      Loc, col = as.numeric(Loc), cex = .6)
#'
#' ############################
#' # plot(y1, y2) for (SGCCA) #
#' ############################
#' layout(t(1:2))
#' plot(result.sgcca$Y[[1]][, 1], result.sgcca$Y[[2]][, 1],
#'      col = "white", xlab = "Y1 (GE)",
#'      ylab = "Y2 (CGH)", main = "Factorial plan of SGCCA")
#' text(result.sgcca$Y[[1]][, 1], result.sgcca$Y[[2]][, 1],
#'      Loc, col = as.numeric(Loc), cex = .6)
#'
#' plot(result.sgcca$Y[[1]][, 1], result.sgcca$Y[[1]][, 2],
#'      col = "white", xlab = "Y1 (GE)",
#'      ylab = "Y2 (GE)", main = "Factorial plan of SGCCA")
#' text(result.sgcca$Y[[1]][, 1], result.sgcca$Y[[1]][, 2],
#'      Loc, col = as.numeric(Loc), cex = .6)
#'
#' # sgcca algorithm with multiple components and different
#' L1 penalties for each components
#' # (-> sparsity is a matrix)
#' init = "random"
#' result.sgcca = sgcca(blocks, connection,
#'                      sparsity = matrix(c(.071,.2, 1, 0.06, 0.15, 1),
#'                      nrow = 2, byrow = TRUE),
#'                      ncomp = c(2, 2, 1), scheme = "factorial", bias = TRUE,
#'                      init = init, verbose = TRUE)
#' # number of non zero elements per dimension
#' apply(result.sgcca$a[[1]], 2, function(x) sum(x!=0))
#'      #(-> 145 non zero elements for a11 and 107 non zero elements for a12)
#' apply(result.sgcca$a[[2]], 2, function(x) sum(x!=0))
#'      #(-> 85 non zero elements for a21 and 52 non zero elements for a22)
#' init = "svd"
#' result.sgcca = sgcca(blocks, connection, sparsity = matrix(c(.071,.2, 1, 0.06, 0.15, 1),
#'                      nrow = 2, byrow = TRUE),
#'                      ncomp = c(2, 2, 1), scheme = "factorial",
#'                      bias = TRUE,
#'                      init = init, verbose = TRUE)
#' }
#'@export sgcca

sgcca <- function(blocks, connection = 1 - diag(length(blocks)),
                  sparsity = rep(1, length(blocks)),
                  ncomp = rep(1, length(blocks)), scheme = "centroid",
                  init = "svd", bias = TRUE, tol = .Machine$double.eps,
                  verbose = FALSE,   quiet = FALSE, na.rm = TRUE,superblock=FALSE){

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
  if (superblock) astar <- matrix(NA, pjs[J], N + 1)

  ###################################################

  if (mode(scheme) != "function") {
    if (verbose) cat("Computation of the SGCCA block components based on the",
                     scheme, "scheme \n")
  }
  if (mode(scheme) == "function" & verbose) {
    cat("Computation of the SGCCA block components based on the g scheme \n")
  }

  ####################################
  # sgcca with 1 component per block #
  ####################################
  if (is.vector(sparsity)) {
    sgcca.result <- sgccak(blocks, connection, sparsity = sparsity,
                           scheme = scheme, init = init, bias = bias,
                           tol = tol, verbose = verbose, quiet = quiet,
                           na.rm = na.rm)
  } else {
    sgcca.result <- sgccak(blocks, connection, sparsity = sparsity[1, ],
                           scheme = scheme, init = init, bias = bias,
                           tol = tol, verbose = verbose, quiet = quiet,
                           na.rm = na.rm)
  }

  for (b in seq_len(J)) {
    Y[[b]][, 1] <- sgcca.result$Y[, b, drop = FALSE]
    a[[b]][, 1] <- sgcca.result$a[[b]]
  }
  ifelse(!superblock,
         astar <- a,
         astar[, 1] <- a[[J]][, 1, drop = FALSE])
  AVE_inner[1]               <- sgcca.result$AVE_inner
  crit[[1]]                  <- sgcca.result$crit

  ##############################################
  #               If any ncomp > 1             #
  #      Determination of SGCCA components     #
  ##############################################
  if (N > 0) {
    R <- blocks
    if(!superblock){
      for (b in seq_len(J)) P[[b]] <- matrix(NA, pjs[b], N)
    }else{
      P <- matrix(NA, pjs[J], N)
    }

    for (n in 2:(N + 1)) {
      if (verbose) message("Computation of the SGCCA block components #", n,
                              " is under progress... \n")

      # Apply deflation
      if (!superblock)
      {
        defl.result <- defl.select(sgcca.result$Y, R, ndefl, n - 1, J, na.rm = na.rm)
        R <- defl.result$resdefl
        for (b in seq_len(J))  {P[[b]][, n - 1] <- defl.result$pdefl[[b]]}
      }
      if (superblock)
      {
        defl.result = deflation(R[[J]], sgcca.result$Y[,J])
        R[[J]] = defl.result$R
        P[, n - 1] <- defl.result$p
        cumsum_pjs = cumsum(pjs)[1:(J - 1)]
        inf_pjs = c(0,cumsum_pjs[1:(J - 2)]) + 1
        for (b in seq_len(J - 1))
        {
          R[[b]] = R[[J]][, c(inf_pjs[b]:cumsum_pjs[b]), drop = FALSE]
          rownames(R[[b]]) = rownames(R[[b]])
          colnames(R[[b]]) = colnames(R[[J]])[c(inf_pjs[b]:cumsum_pjs[b])]
        }
      }

      if (is.vector(sparsity)) {
        sgcca.result <- sgccak(R, connection, sparsity = sparsity,
                               scheme = scheme, init = init, bias = bias,
                               tol = tol, verbose = verbose, quiet = quiet,
                               na.rm = na.rm)
      } else {
        sgcca.result <- sgccak(R, connection, sparsity = sparsity[n, ],
                               scheme = scheme, init = init, bias = bias,
                               tol = tol, verbose = verbose, quiet = quiet,
                               na.rm = na.rm)
      }

      AVE_inner[n] <- sgcca.result$AVE_inner
      crit[[n]] <- sgcca.result$crit


      for (b in seq_len(J))  {
        Y[[b]][, n] <- sgcca.result$Y[, b]
        a[[b]][, n] <- sgcca.result$a[[b]]
        if (!superblock)
          astar[[b]][, n] <- sgcca.result$a[[b]] -
          astar[[b]][, (1:(n - 1)), drop = F] %*%
          drop(crossprod(a[[b]][, n], P[[b]][, 1:(n - 1), drop = FALSE]))
        else
          astar[, n] <- sgcca.result$a[[J]] -
          astar[, (1:(n - 1)), drop = F] %*%
          drop(t(a[[J]][, n]) %*% P[, 1:(n - 1), drop = F])
      }

      if (!quiet){
        for (q in which(n < ndefl)) {
          if (sum(sgcca.result$a[[q]] != 0) <= 1) {
            warning(sprintf("Deflation failed because only one variable was
                            selected for block ",q,"! \n"))

          }
        }
      }
    }
  }

  for (b in seq_len(J)) {
    #Average Variance Explained (AVE) per block
    AVE_X[[b]] =  apply(cor(blocks[[b]], Y[[b]], use = "pairwise.complete.obs")^2, 2,
                      mean, na.rm = TRUE)
  }

  #AVE outer
  outer = matrix(unlist(AVE_X), nrow = max(ncomp))
  AVE_outer <- as.numeric((outer %*% pjs)/sum(pjs))

  AVE_X = shave(AVE_X, ncomp)

  AVE <- list(AVE_X = AVE_X, AVE_outer = AVE_outer, AVE_inner = AVE_inner)

  if (N == 0) crit = unlist(crit)

  out <- list(Y = Y,
              a = a,
              astar = astar,
              crit = crit,
              AVE = AVE)

  class(out) <- "sgcca"
  return(out)

}
