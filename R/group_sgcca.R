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
#' result.rgcca = rgccad(blocks, connection, tau, ncomp = c(2, 2, 1),
#' scheme = "factorial", verbose = TRUE)
#' # sgcca algorithm
#' result.sgcca = sgcca(blocks, connection, sparsity = c(.071,.2, 1), ncomp = c(2, 2, 1),
#'                      scheme = "centroid", verbose = TRUE)
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

group_sgcca <- function (blocks, connection = 1-diag(length(blocks)), sparsity = rep(1, length(blocks)), group_sparsity,
                   ncomp = rep(1, length(blocks)), scheme = "centroid", scale = TRUE,
                   init = "svd", bias = TRUE, tol = .Machine$double.eps,
                   verbose = FALSE, scale_block = TRUE, prescaling = FALSE,
                   quiet = FALSE){

  ndefl <- ncomp-1
  N <- max(ndefl)
  J <- length(blocks)
  pjs <- sapply(blocks,NCOL)
  nb_ind <- NROW(blocks[[1]])
  AVE_X = list()
  AVE_outer <- rep(NA,max(ncomp))

  if ( any(ncomp < 1) ) stop_rgcca("One must compute at least one component per
                                   block!")
  if (any(ncomp-pjs > 0)) stop_rgcca("For each block, choose a number of
                                     components smaller than the number of
                                     variables!")

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
    if ((scheme != "horst" ) & (scheme != "factorial") & (scheme != "centroid")) {
      stop_rgcca("Choose one of the three following schemes: horst, centroid,
                 factorial or design the g function")
    }
    if (verbose) cat("Computation of the group SGCCA block components based on the",
                     scheme, "scheme \n")
  }
  if (mode(scheme) == "function" & verbose) {
    cat("Computation of the group SGCCA block components based on the g scheme \n")
  }


  #-------------------------------------------------------
  if(!prescaling)
    blocks=scaling(blocks, scale = scale, bias = bias, scale_block = scale_block)

  ####################################
  # sgcca with 1 component per block #
  ####################################

  # ndefl number of deflation per block
  ndefl <- ncomp-1
  N <- max(ndefl)
  J <- length(blocks)
  pjs <- sapply(blocks,NCOL)
  nb_ind <- NROW(blocks[[1]])
  AVE_X = list()
  AVE_outer <- rep(NA,max(ncomp))
  if (N == 0) {
    result <- group_sgccak(blocks, connection, sparsity, group_sparsity = group_sparsity, scheme, init = init, bias = bias,
                           tol = tol, verbose = verbose, quiet = quiet)
    # No deflation (No residual matrices generated).
    Y <- NULL
    for (b in 1:J) Y[[b]] <- result$Y[,b, drop = FALSE]
    #Average Variance Explained (AVE) per block
    for (j in 1:J) AVE_X[[j]] =  mean(cor(blocks[[j]], Y[[j]],
                                          use="pairwise.complete.obs")^2,
                                      na.rm=TRUE)

    #AVE outer
    AVE_outer <- sum(pjs * unlist(AVE_X))/sum(pjs)

    AVE <- list(AVE_X = AVE_X,
                AVE_outer = AVE_outer,
                AVE_inner = result$AVE_inner)

    a <- lapply(result$a, cbind)

    for (b in 1:J) {
      rownames(a[[b]]) = colnames(blocks[[b]])
      rownames(Y[[b]]) = rownames(blocks[[b]])
      colnames(Y[[b]]) = "comp1"
    }

    out <- list(Y=Y, a=a, astar=a,
                connection=connection, scheme=scheme, sparsity=sparsity, ncomp=ncomp,
                crit = result$crit[length(result$crit)],
                AVE = AVE)
    class(out) <- "group_sgcca"
    return(out)
  }

  ##################
  # Initialization #
  ##################

  Y <- NULL
  R <- blocks
  P <- a <- astar <- NULL
  crit <- list()
  AVE_inner <- rep(NA,max(ncomp))

  for (b in 1:J) P[[b]] <- a[[b]] <- astar[[b]] <- matrix(NA,pjs[[b]],N+1)
  for (b in 1:J) Y[[b]] <- matrix(NA,nb_ind,N+1)

  ##############################################
  #               If any ncomp > 1             #
  #      Determination of SGCCA components     #
  ##############################################


  for (n in 1:N) {
    if (verbose) cat(paste0("Computation of the group SGCCA block components #", n,
                            " is under progress... \n"))
    if(is.vector(sparsity)){
      sgcca.result <- group_sgccak(R, connection, sparsity = sparsity , group_sparsity = group_sparsity, scheme=scheme,
                                   init = init, bias = bias, tol = tol,
                                   verbose = verbose, quiet = quiet)
    } else{
      sgcca.result <- group_sgccak(R, connection, sparsity = sparsity[n, ] , group_sparsity = group_sparsity, scheme = scheme,
                                   init = init, bias = bias, tol = tol,
                                   verbose = verbose, quiet = quiet)
    }
    AVE_inner[n] <- sgcca.result$AVE_inner
    crit[[n]] <- sgcca.result$crit


    for (b in 1:J) Y[[b]][,n] <- sgcca.result$Y[ ,b]
    for (q in which(n <= ndefl)) if(sum(sgcca.result$a[[q]]!=0) <= 1)
    {
      if(!quiet)
      {
        warning("Deflation failed because only one variable was
                selected for block ",q,"!")
      }
    }
    defla.result <- defl.select(sgcca.result$Y, R, ndefl, n, nbloc = J)
    R <- defla.result$resdefl
    for (b in 1:J) {
      P[[b]][,n] <- defla.result$pdefl[[b]]
      a[[b]][,n] <- sgcca.result$a[[b]]
    }

    if (n==1) {
      for (b in 1:J) astar[[b]][,n] <- sgcca.result$a[[b]]
    }
    else {
      for (b in 1:J) astar[[b]][,n] <- sgcca.result$a[[b]] -
          astar[[b]][,(1:n-1),drop=F] %*%
          drop( t(a[[b]][,n]) %*% P[[b]][,1:(n-1),drop=F] )
    }
  }
  if (verbose) cat(paste0("Computation of the group SGCCA block components #", N+1, "
                            is under progress...\n"))
  if(is.vector(sparsity)) {
    sgcca.result <- group_sgccak(R, connection, sparsity = sparsity, group_sparsity = group_sparsity, scheme=scheme,
                                 init = init, bias = bias, tol = tol,
                                 verbose = verbose, quiet = quiet)
  } else{
    sgcca.result <- group_sgccak(R, connection, sparsity = sparsity[N+1, ], group_sparsity = group_sparsity, scheme=scheme,
                                 init = init, bias = bias, tol = tol,
                                 verbose = verbose, quiet = quiet)
  }
  AVE_inner[max(ncomp)] <- sgcca.result$AVE_inner

  crit[[N+1]] <- sgcca.result$crit
  for (q in which(N == ndefl)) if(sum(sgcca.result$a[[q]]!=0) <= 1)
  {
    if(!quiet)
    {
      warning("Deflation failed because only one variable was
                selected for block ",q,"!")
    }
  }
  for (b in 1:J) {
    Y[[b]][,N+1]     <- sgcca.result$Y[, b]
    a[[b]][,N+1]     <- sgcca.result$a[[b]]
    astar[[b]][,N+1] <- sgcca.result$a[[b]] -  astar[[b]][,(1:N),drop=F] %*%
      drop(t(a[[b]][, N+1]) %*% P[[b]][,1:N,drop=F])
    rownames(a[[b]]) = rownames(astar[[b]]) = colnames(blocks[[b]])
    rownames(Y[[b]]) = rownames(blocks[[b]])
    colnames(Y[[b]]) = paste0("comp", 1:max(ncomp))
  }

  shave.matlist <- function(mat_list, nb_cols)
    mapply(function(m, nbcomp) m[, 1:nbcomp, drop = FALSE],
           mat_list, nb_cols, SIMPLIFY=FALSE)
  shave.veclist <- function(vec_list, nb_elts)
    mapply(function(m, nbcomp) m[1:nbcomp],
           vec_list, nb_elts, SIMPLIFY=FALSE)

  #Average Variance Explained (AVE) per block
  for (j in 1:J) AVE_X[[j]] =  apply(cor(blocks[[j]], Y[[j]],
                                         use="pairwise.complete.obs")^2, 2,
                                     function(x) {return(mean(x,is.na=TRUE))})

  #AVE outer
  outer = matrix(unlist(AVE_X), nrow = max(ncomp))
  for (j in 1:max(ncomp))
    AVE_outer[j] <- sum(pjs * outer[j, ], na.rm=T)/sum(pjs)

  Y = shave.matlist(Y, ncomp)
  AVE_X = shave.veclist(AVE_X, ncomp)

  AVE <- list(AVE_X = AVE_X,
              AVE_outer = AVE_outer,
              AVE_inner = AVE_inner)

  out <- list(Y = shave.matlist(Y, ncomp),
              a = shave.matlist(a, ncomp),
              astar = shave.matlist(astar, ncomp),
              crit = crit,
              AVE = AVE)

  class(out) <- "group_sgcca"
  return(out)

}
