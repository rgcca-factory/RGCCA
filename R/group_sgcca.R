#' Group SGCCA extends SGCCA to address the issue of group variable selection.
#' Specifically, each block is either combined with an L1-penalty similarly to
#' SGCCA or with a group-LASSO penalty in order to perform group variable selection.
#' @inheritParams select_analysis
#' @inheritParams rgccad
#' @inheritParams sgcca
#' @param group_sparsity A list of length \eqn{J} (the number of blocks). Each element
#' of this list is either 'NULL', when no group-LASSO constraint is enforced to the
#' considered block. In this case, a L1 constraint is applied with the sparse
#' parameter defined in 'sparsity'. Or, each element of this list can be also a list
#' defining the groups that have to be used for the group-LASSO constraint. For now,
#' to define a group, you have to mention all the column index associated with
#' the variables contained in this group. Again, in this case, the sparse parameter
#' used for the group-LASSO constraint is the one defined in 'sparsity'. However,
#' if \eqn{G_j} is the number of groups defined for block \eqn{j}, the amount of
#' sparsity this time have to vary between \eqn{1/sqrt(G_j)} and 1 (larger values
#' of sparsity correspond to less penalization).
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
#'
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
#' group_GE  = NULL
#' group_CGH = cut(1:NCOL(blocks$CGH), breaks = 100)
#' group_CGH = lapply(levels(group_CGH), function(x) which(group_CGH %in% x))
#' group_y   = NULL
#' group_sparsity = list(group_GE  = group_GE,
#'                       group_CGH = group_CGH,
#'                       group_y   = group_y)
#'
#' # group sgcca algorithm
#' result.group_sgcca = rgcca(blocks = blocks, connection = connection,
#'                            sparsity = c(.071,0.2, 1), group_sparsity = group_sparsity, ncomp = c(2, 2, 1),
#'                            scheme = "centroid", verbose = TRUE, method = "group_sgcca")
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
#' ##################################
#' # plot(y1, y2) for (group SGCCA) #
#' ##################################
#' layout(t(1:2))
#' plot(result.group_sgcca$Y[[1]][, 1], result.group_sgcca$Y[[2]][, 1],
#'      col = "white", xlab = "Y1 (GE)",
#'      ylab = "Y2 (CGH)", main = "Factorial plan of group SGCCA")
#' text(result.group_sgcca$Y[[1]][, 1], result.group_sgcca$Y[[2]][, 1],
#'      Loc, col = as.numeric(Loc), cex = .6)
#'
#' plot(result.group_sgcca$Y[[1]][, 1], result.group_sgcca$Y[[1]][, 2],
#'      col = "white", xlab = "Y1 (GE)",
#'      ylab = "Y2 (GE)", main = "Factorial plan of group SGCCA")
#' text(result.group_sgcca$Y[[1]][, 1], result.group_sgcca$Y[[1]][, 2],
#'      Loc, col = as.numeric(Loc), cex = .6)
#' }
#'@export group_sgcca

group_sgcca <- function(blocks, connection = 1 - diag(length(blocks)),
                        sparsity = rep(1, length(blocks)), group_sparsity,
                        ncomp = rep(1, length(blocks)), scheme = "centroid",
                        init = "svd", bias = TRUE, tol = .Machine$double.eps,
                        verbose = FALSE,   quiet = FALSE, na.rm = TRUE){

  # ndefl number of deflation per block
  ndefl <- ncomp - 1
  N <- max(ndefl)
  J <- length(blocks)
  pjs <- sapply(blocks,NCOL)
  nb_ind <- NROW(blocks[[1]])
  AVE_X = list()
  AVE_outer <- rep(NA,max(ncomp))

  Y <- NULL
  a <- astar <- P <- NULL
  crit <- list()
  AVE_inner <- rep(NA,max(ncomp))

  for (b in 1:J) a[[b]] <- astar[[b]] <- matrix(NA, pjs[[b]], N + 1)
  for (b in 1:J) Y[[b]] <- matrix(NA,nb_ind, N + 1)

  active_group_sparsity = !sapply(group_sparsity, is.null)
  nb_groups             = sapply(group_sparsity[active_group_sparsity], length)
  if (is.vector(sparsity)){
    if (any(sparsity[!active_group_sparsity] < 1/sqrt(pjs[!active_group_sparsity]) | sparsity[!active_group_sparsity] > 1 ))
      stop_rgcca("L1 constraints (sparsity) must vary between 1/sqrt(p_j) and 1.")
    if (any(sparsity[active_group_sparsity] <  1/sqrt(nb_groups) | sparsity[active_group_sparsity] > 1 ))
      stop_rgcca("LG constraints (sparsity) must vary between 1/sqrt(numberOfGroups) and 1.")
  }

  if (is.matrix(sparsity)){
    if (any(apply(as.matrix(sparsity[, !active_group_sparsity]), 1, function(x) any(x < 1/sqrt(pjs[!active_group_sparsity])))))
      stop_rgcca("L1 constraints (sparsity) must vary between 1/sqrt(p_j)
                 and 1.")
    if (any(apply(as.matrix(sparsity[, active_group_sparsity]), 1, function(x) any(x < 1/sqrt(nb_groups)))))
      stop_rgcca("LG constraints (sparsity) must vary between 1/sqrt(numberOfGroups)
                 and 1.")
  }

  ###################################################

  if (mode(scheme) != "function") {
    if (verbose) cat("Computation of the group SGCCA block components based on the",
                     scheme, "scheme \n")
  }
  if (mode(scheme) == "function" & verbose) {
    cat("Computation of the group SGCCA block components based on the g scheme \n")
  }

  ####################################
  # sgcca with 1 component per block #
  ####################################
  if (is.vector(sparsity)) {
    sgcca.result <- group_sgccak(blocks, connection, sparsity = sparsity,
			   group_sparsity = group_sparsity,
                           scheme = scheme, init = init, bias = bias,
                           tol = tol, verbose = verbose, quiet = quiet,
                           na.rm = na.rm)
  } else {
    sgcca.result <- group_sgccak(blocks, connection, sparsity = sparsity[1, ],
			   group_sparsity = group_sparsity,
                           scheme = scheme, init = init, bias = bias,
                           tol = tol, verbose = verbose, quiet = quiet,
                           na.rm = na.rm)
  }

  for (b in 1:J) Y[[b]][, 1] <- sgcca.result$Y[, b, drop = FALSE]
  for (b in 1:J) a[[b]][, 1] <- sgcca.result$a[[b]]
  astar                      <- a
  AVE_inner[1]               <- sgcca.result$AVE_inner
  crit[[1]]                  <- sgcca.result$crit

  ##############################################
  #               If any ncomp > 1             #
  #      Determination of SGCCA components     #
  ##############################################
  if (N > 0) {
    R <- blocks
    for (b in 1:J) P[[b]] <- matrix(NA, pjs[[b]], N)

    for (n in 2:(N + 1)) {
      if (verbose) cat(paste0("Computation of the SGCCA block components #", n,
                              " is under progress... \n"))

      # Apply deflation
      defla.result <- defl.select(sgcca.result$Y, R, ndefl, n - 1, J, na.rm = na.rm)
      R <- defla.result$resdefl
      for (b in 1:J) P[[b]][, n - 1] <- defla.result$pdefl[[b]]

      if (is.vector(sparsity)) {
        sgcca.result <- group_sgccak(R, connection, sparsity = sparsity,
			       group_sparsity = group_sparsity,
                               scheme = scheme, init = init, bias = bias,
                               tol = tol, verbose = verbose, quiet = quiet,
                               na.rm = na.rm)
      } else {
        sgcca.result <- group_sgccak(R, connection, sparsity = sparsity[n, ],
                               group_sparsity = group_sparsity,
                               scheme = scheme, init = init, bias = bias,
                               tol = tol, verbose = verbose, quiet = quiet,
                               na.rm = na.rm)
      }

      AVE_inner[n] <- sgcca.result$AVE_inner
      crit[[n]] <- sgcca.result$crit


      for (b in 1:J) Y[[b]][, n] <- sgcca.result$Y[, b]
      for (b in 1:J) a[[b]][, n] <- sgcca.result$a[[b]]
      for (b in 1:J) astar[[b]][, n] <- sgcca.result$a[[b]] -
        astar[[b]][, (1:(n - 1)), drop = F] %*%
        drop( t(a[[b]][, n]) %*% P[[b]][, 1:(n - 1), drop = F] )

      for (q in which(n < ndefl)) if (sum(sgcca.result$a[[q]] != 0) <= 1)
      {
        if (!quiet)
        {
          warning(sprintf("Deflation failed because only one variable was
                            selected for block ",q,"! \n"))

        }
      }
    }
  }

  for (b in 1:J) {
    rownames(a[[b]]) = rownames(astar[[b]]) = colnames(blocks[[b]])
    rownames(Y[[b]]) = rownames(blocks[[b]])
    colnames(Y[[b]]) = paste0("comp", 1:max(ncomp))
  }

  #Average Variance Explained (AVE) per block
  for (j in 1:J) AVE_X[[j]] =  apply(
    cor(blocks[[j]], Y[[j]], use = "pairwise.complete.obs")^2, 2,
    function(x) {return(mean(x, is.na = TRUE))})

  #AVE outer
  outer = matrix(unlist(AVE_X), nrow = max(ncomp))
  for (j in 1:max(ncomp))
    AVE_outer[j] <- sum(pjs * outer[j, ], na.rm = T) / sum(pjs)

  Y = shave.matlist(Y, ncomp)
  AVE_X = shave.veclist(AVE_X, ncomp)

  AVE <- list(AVE_X = AVE_X, AVE_outer = AVE_outer, AVE_inner = AVE_inner)

  if (N == 0) crit = unlist(crit)

  out <- list(Y = shave.matlist(Y, ncomp),
              a = shave.matlist(a, ncomp),
              astar = shave.matlist(astar, ncomp),
              crit = crit,
              AVE = AVE)

  class(out) <- "group_sgcca"
  return(out)

}
