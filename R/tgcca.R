tgcca <- function(A, C = 1-diag(length(A)), tau = rep(1, length(A)),
                  ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,
                  init="svd", bias = TRUE, tol = 1e-8, verbose=FALSE,
                  scale_block = TRUE, kronecker_covariance = FALSE,
                  ranks = rep(1, length(A)), prescaling = FALSE, quiet = FALSE,
                  n_run = 1, n_cores = 1) {

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
  pjs    <- sapply(DIM, function(x) prod(x[-1]))
  nb_ind <- DIM[[1]][1]

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

  ###################################################

  if (mode(scheme) != "function") {
    if ((scheme != "horst" ) & (scheme != "factorial") & (scheme != "centroid")) {
      stop_rgcca("Choose one of the three following schemes: horst, centroid,
                 factorial or design the g function")
    }
    if (verbose) cat("Computation of the TGCCA block components based on the",
                     scheme, "scheme \n")
  }
  if (mode(scheme) == "function" & verbose) {
    cat("Computation of the TGCCA block components based on the g scheme \n")
  }


  #-------------------------------------------------------
  if(!prescaling)
    A=scaling(A, scale = scale, bias = bias, scale_block = scale_block)

  # Matricization (mode-1)
  A_m = lapply(1:J, function(x) {
    m = matrix(as.vector(A[[x]]), nrow = nb_ind)
    rownames(m) = rownames(A[[x]])
    if (!is.null(dimnames(A[[x]]))) {
      grid        = do.call(expand.grid, dimnames(A[[x]])[-1])
      colnames(m) = do.call(paste, c(grid, sep = " x "))
    }
    return(m)
  })

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
  Y     <- P    <- a       <- astar   <- list()
  weights <- replicate(J, c(), simplify = FALSE)


  ranks = matrix(ranks, nrow = max(ncomp), ncol = J, byrow = T)
  colnames(ranks) = names(A)
  for (d in B_2D) ranks[, d] <- 1
  for (d in 1:J) P[[d]]  <- a[[d]] <- astar[[d]] <- matrix(NA,pjs[[d]],N+1)

  for (d in B_nD) {
    factors[[d]] = list()
    for (f in 1:(LEN[[d]] - 1)) {
      factors[[d]][[f]] = matrix(NA, DIM[[d]][[f + 1]], sum(ranks[, d]))
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
  ### Determination of TGCCA components ###
  #########################################
  for (n in 1:(N+1)) {
    if (verbose) cat(paste0("Computation of the TGCCA block components #", n,
                            " is under progress...\n"))
    # n_random_starts
    mgcca.result = tgccak(
      A                       = R,
      A_m                     = R_m,
      C                       = C,
      tau                     = tau_mat[n, ],
      scheme                  = scheme,
      init                    = init,
      bias                    = bias,
      tol                     = tol,
      verbose                 = verbose,
      kronecker_covariance    = kronecker_covariance,
      ranks                   = ranks[n, ],
      n_run                   = n_run,
      n_cores                 = n_cores
    )

    # Store tau, AVE_inner, crit
    if (!is.numeric(tau)) tau_mat[n, ] = mgcca.result$tau
    AVE_inner[n] = mgcca.result$AVE_inner
    crit[[n]]    = mgcca.result$crit

    # Store Y, a, factors and weights
    for (d in 1:J) Y[[d]][,n] = mgcca.result$Y[ , d]
    for (d in 1:J) a[[d]][,n] = mgcca.result$a[[d]]
    for (d in B_nD) {
      idx = seq(c(0, cumsum(ranks[, d]))[n] + 1, c(0, cumsum(ranks[, d]))[n + 1])
      for (f in 1:(LEN[[d]] - 1)) {
        factors[[d]][[f]][, idx] = mgcca.result$factors[[d]][[f]]
      }
      weights[[d]] = c(weights[[d]], mgcca.result$weights[[d]])
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
        astar[[d]][,n] = mgcca.result$a[[d]] - astar[[d]][,1:(n-1), drop=F] %*%
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
    rownames(a[[d]]) = rownames(astar[[d]]) = colnames(A_m[[d]])
    rownames(Y[[d]]) = rownames(A[[d]])
    colnames(Y[[d]]) = paste0("comp", 1:max(ncomp))
  }

  # Average Variance Explained (AVE) per block
  for (j in 1:J)
    AVE_X[[j]] = apply(cor2(A[[j]], Y[[j]], use="pairwise.complete.obs")^2,
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
    factors[[d]] = shave.matlist(factors[[d]], sum(ranks[, d]))
  }

  # AVE
  AVE = list(AVE_X           = shave.veclist(AVE_X, ncomp),
             AVE_outer_model = AVE_outer,
             AVE_inner_model = AVE_inner)

  # output
  out = list(blocks  = A,
             Y       = shave.matlist(Y, ncomp),
             a       = shave.matlist(a, ncomp),
             factors = factors,
             weights = weights,
             astar   = shave.matlist(astar, ncomp),
             crit    = crit,
             AVE     = AVE,
             tau     = tau_mat)

  class(out) = "mgcca"
  return(out)
}
