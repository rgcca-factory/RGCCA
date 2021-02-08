#' MGCCA extends RGCCA to address the issue of tensor structured data. 
#' Specifically, RGCCA is combined with a Kronecker constraint that gives rise 
#' to Multiway GCCA (MGCCA) which is implemented in the function mgcca().
#' Given \eqn{J} arrays \eqn{X_1, X_2, ..., X_J}, that represent
#' \eqn{J} sets of variables observed on the same set of \eqn{n} individuals.
#' The arrays \eqn{X_1, X_2, ..., X_J} must have the same dimension on the 
#' first dimension, but may (and usually will) have different numbers of modes 
#' and mode dimensions. Blocks are not necessarily fully connected within the 
#' MGCCA framework. Hence MGCCA requires the construction (user specified) of a 
#' design matrix (\eqn{C}) that characterizes the connections between blocks. 
#' Elements of the symmetric design matrix \eqn{C = (c_{jk})} are equal to 1 if 
#' block \eqn{j} and block \eqn{k} are connected, and 0 otherwise. The MGCCA 
#' algorithm is very similar to the RGCCA algorithm and keeps the same monotone 
#' convergence properties (i.e. the bounded criteria to be maximized increases 
#' at each step of the iterative procedure and hits at convergence a stationary 
#' point).
#' Moreover, using a deflation strategy, mgcca() enables the computation of
#' several MGCCA block components (specified by ncomp) for each block. Block
#' components for each block are guaranteed to be orthogonal when using this
#' deflation strategy. The so-called symmetric deflation is considered in this
#' implementation, i.e. each block is deflated with respect to its own
#' component. Moreover, we stress that the numbers of components per block
#' could differ from one block to another.
#' @inheritParams select_analysis
#' @inheritParams rgccaNa
#' @inheritParams rgccad
#' @param M_regularisation If not NULL, a list of \eqn{J} elements. Each element
#' of \eqn{M_regularisation} is either NULL or a list of symmetric positive
#' definite regularization matrices. There must be as many matrices as modes
#' on the corresponding block and their dimensions must match the dimensions of
#' the corresponding modes. A change of variable is done at the beginning and at
#' the end of the MGCCA algorithm to apply regularization.
#' @param nstart Number of random starts to avoid finding a local optimum.
#' @param nstart_at_comp_2 Boolean to allow multiple random starts when 
#' estimating the 2nd component.
#' @param ranks A vector of \eqn{J} elements that gives the number of terms in 
#' the Kronecker constraint (hence the rank of the estimated tensor) for each
#' block.
#' @return \item{Y}{A list of \eqn{J} elements. Each element of \eqn{Y} is a
#' matrix that contains the analysis components for the corresponding block.}
#' @return \item{a}{A list of \eqn{J} elements. Each element of \eqn{a} is a
#' matrix that contains the outer weight vectors for each block.}
#' @return \item{astar}{A list of \eqn{J} elements. Each element of astar is a
#' matrix defined as Y[[j]][, h] = A[[j]]\%*\%astar[[j]][, h]}
#' @return \item{factors}{A list of \eqn{J} elements. If bloc \eqn{j} is a
#' tensor of order \eqn{d}, element \eqn{j} of \eqn{factors} is a list with
#' \eqn{d} elements and each element is a matrix that contains the outer weight
#' vectors for each block.
#' @return \item{call}{Call of the function}
#' @return \item{crit}{A vector of integer that contains for each component the
#' values of the analysis criteria across iterations.}
#' @return \item{AVE}{A list of numerical values giving the indicators of model
#' quality based on the Average Variance Explained (AVE): AVE(for each block),
#' AVE(outer model), AVE(inner model).}
#' @return \item{A}{The eventually preprocessed entry data.}
#' @references Soon.
#' @title Multiway Generalized Canonical Correlation Analysis (MGCCA)
#'@export mgcca


mgcca <- function(A, C = 1-diag(length(A)), tau = rep(1, length(A)), 
                  ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,
                  init="svd", bias = TRUE, tol = .Machine$double.eps, 
                  verbose=FALSE, scale_block = TRUE, M_regularisation = NULL, 
                  ranks = rep(1, length(A)), prescaling = FALSE, quiet = FALSE,
                  nstart = 1, nstart_at_comp_2 = FALSE) {

  call=list(A = A, C = C,  ncomp = ncomp, scheme = scheme, scale = scale,
            init = init, bias = bias, tol = tol, verbose = verbose,
            scale_block = scale_block, M_regularisation = M_regularisation,
            ranks = ranks, prescaling = prescaling)

  # Number of blocks
  J      = length(A)

  # List of 3D Tensors and 2D matrix
  DIM    <- lapply(A, dim)
  LEN    <- sapply(DIM, length)
  B_nD   <- which(LEN >= 3)   # Store which blocks are higher order tensors
  B_2D   <- which(LEN == 2)   # Store which blocks are 2D
  B_0D   <- which(LEN == 0)   # Store which blocks are 1D (stored as 0D)

  # Convert vectors to one-column matrices
  if (length(B_0D) !=0){
    for (i in B_0D){
      A[[i]]   = as.matrix(A[[i]])
      DIM[[i]] = dim(A[[i]])
    }
    B_2D = c(B_2D, B_0D)
  }

  # Dimensions of each block
  pjs <- sapply(DIM, function(x) prod(x[-1]))
  nb_ind   <- DIM[[1]][1]

  # Matricization (mode-1)
  A_m = lapply(1:J, function(x) matrix(as.vector(A[[x]]), nrow = nb_ind))

  if ( any(ncomp < 1) ) stop_rgcca("One must compute at least one component per
                                   block!")
  if (any(ncomp-pjs > 0)) stop_rgcca("For each block, choose a number of
                                     components smaller than the number of
                                     variables!")

  # Shrinkage parameters
  if(is.null(tau) & verbose){
    cat("The shrinkage parameters have been
              automatically set to 1 for all blocks \n")
    tau = rep(1, NCOL(C))
  }

  # Multiple starts are random
  if (nstart > 1) init = "random"

  ###################################################

  if (mode(scheme) != "function") {
    if ((scheme != "horst" ) & (scheme != "factorial") & (scheme != "centroid")) {
      stop_rgcca("Choose one of the three following schemes: horst, centroid,
                 factorial or design the g function")
    }
    if (verbose) cat("Computation of the SGCCA block components based on the",
                     scheme, "scheme \n")
  }
  if (mode(scheme) == "function" & verbose) {
    cat("Computation of the SGCCA block components based on the g scheme \n")
  }


  #-------------------------------------------------------
  if(!prescaling)
    A=scaling(A, scale = scale, bias = bias, scale_block = scale_block)

  ######################
  ### Initialization ###
  ######################
  AVE_outer        = vector()
  ndefl            = ncomp-1
  N                = max(ndefl)
  AVE_inner        = rep(NA,max(ncomp))
  R                = A
  R_m              = A_m

  AVE_X <- crit <- factors <- list()
  Y     <- P    <- a       <- astar <- NULL

  for (d in 1:J) P[[d]]  <- a[[d]] <- astar[[d]] <- matrix(NA,pjs[[d]],N+1)

  for (d in B_nD) {
    factors[[d]] = list()
    for (f in 1:(LEN[[d]] - 1)) {
      factors[[d]][[f]] = matrix(NA, DIM[[d]][[f + 1]], (N+1) * ranks[[d]])
    }
  }

  for (d in 1:J)  Y[[d]] = matrix(NA, nb_ind, N+1)

  if (is.numeric(tau))
  {
    tau_mat = tau
    if(is.vector(tau_mat)) {
      tau_mat = matrix(tau_mat, nrow = max(ncomp), ncol = J, byrow = TRUE)
    }
    if(is.matrix(tau_mat)) colnames(tau_mat) = names(A)
  }
  else
  {
    tau_mat = matrix(NA, max(ncomp), J)
    colnames(tau_mat) = names(A)
  }

  #########################################
  ### Determination of MGCCA components ###
  #########################################
  for (n in 1:(N+1)) {
    if (verbose) cat(paste0("Computation of the MGCCA block components #", n,
                            " is under progress...\n"))
    # n_random_starts
    for (start in 1:nstart){
      # MGCCA algorithm
      tmp_mgcca.result = mgccak(
        A          = R,
        A_m        = R_m,
        C          = C,
        tau        = tau_mat[n, ],
        scheme     = scheme,
        init       = init,
        bias       = bias,
        tol        = tol,
        verbose    = verbose,
        M_regularisation = M_regularisation,
        ranks            = ranks
      )
      tmp_crit = tmp_mgcca.result$crit[length(tmp_mgcca.result$crit)]
      if(start == 1){
        mgcca.result = tmp_mgcca.result
        best_crit    = tmp_crit
      }else{
        if(tmp_crit > best_crit){
          mgcca.result = tmp_mgcca.result
          best_crit    = tmp_crit
        }
      }
    }
    if (!nstart_at_comp_2) nstart = 1
    nstart_at_comp_2 = FALSE
    
    # Store tau, AVE_inner, crit
    if (!is.numeric(tau)) tau_mat[n, ] = mgcca.result$tau
    AVE_inner[n] = mgcca.result$AVE_inner
    crit[[n]]    = mgcca.result$crit

    # Store Y, a, and factors
    for (d in 1:J) Y[[d]][,n] = mgcca.result$Y[ , d]
    for (d in 1:J) a[[d]][,n] = mgcca.result$a[[d]]
    for (d in B_nD) {
      for (f in 1:(LEN[[d]] - 1)) {
        idx                      = seq((n - 1) * ranks[[d]] + 1, n * ranks[[d]])
        factors[[d]][[f]][, idx] = mgcca.result$factors[[d]][[f]]
      }
    }

    # Deflation procedure
    defla.result = defl.select(mgcca.result$Y, R, ndefl, n, nbloc = J)
    R            = defla.result$resdefl
    R_m          = NULL # Let the inner function take care of the matricization

    # Store projection matrices for deflation
    for (d in 1:J) P[[d]][,n] = defla.result$pdefl[[d]]

    # Compute astar
    for (d in 1:J){
      if (n == 1){
        astar[[d]][,n] = mgcca.result$a[[d]]
      }else{
        astar[[d]][,n] = mgcca.result$a[[d]] - astar[[d]][,(1:n-1), drop=F] %*%
          drop( t(a[[d]][,n]) %*% P[[d]][,1:(n-1),drop=F] )
      }
    }
  }

  #############
  ### Names ###
  #############
  for (d in B_nD) {
    for (f in 1:(LEN[[d]] - 1)) {
      rownames(factors[[d]][[f]]) = dimnames(A[[d]])[[f + 1]]
    }
  }
  for (d in 1:J){
    rownames(a[[d]]) = rownames(astar[[d]]) = colnames(A[[d]])
    rownames(Y[[d]]) = rownames(A[[d]])
    colnames(Y[[d]]) = paste0("comp", 1:max(ncomp))
  }

  # Average Variance Explained (AVE) per block
  for (j in 1:J)
    AVE_X[[j]] = apply(cor(A_m[[j]], Y[[j]], use="pairwise.complete.obs")^2,
                       2, mean, na.rm = TRUE)
  # AVE outer
  if (N == 0) {
    AVE_outer = sum(pjs * unlist(AVE_X)) / sum(pjs)
  }else{
    outer = matrix(unlist(AVE_X), nrow = max(ncomp))
    for (j in 1:max(ncomp)) AVE_outer[j] = sum(pjs * outer[j, ])/sum(pjs)
  }

  # Remove unused components
  for (d in B_nD) {
    factors[[d]] = shave.matlist(factors[[d]], ncomp * ranks[d])
  }

  # AVE
  AVE = list(AVE_X           = AVE_X,
             AVE_outer_model = AVE_outer,
             AVE_inner_model = AVE_inner)
  
  call$tau = tau_mat

  # output
  out = list(Y       = shave.matlist(Y, ncomp),
             a       = shave.matlist(a, ncomp),
             factors = factors,
             astar   = shave.matlist(astar, ncomp),
             crit    = crit,
             AVE     = AVE,
             A       = A,
             call    = call)

  class(out) = "mgcca"
  return(out)
}
