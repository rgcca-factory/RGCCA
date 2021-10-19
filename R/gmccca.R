gmgcca <- function(A, C = 1-diag(length(A)), tau = rep(1, length(A)),
                  ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,
                  init="svd", bias = TRUE, tol = 1e-8, verbose=FALSE,
                  scale_block = TRUE, regularisation_matrices = NULL,
                  ranks = rep(1, length(A)), prescaling = FALSE, quiet = FALSE) {

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
    if (verbose) cat("Computation of the gMGCCA block components based on the",
                     scheme, "scheme \n")
  }
  if (mode(scheme) == "function" & verbose) {
    cat("Computation of the gMGCCA block components based on the g scheme \n")
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

  AVE_X <- crit <- factors <- weights <- list()
  Y     <- P    <- a       <- astar   <- list()

  for (d in B_2D) ranks[d] <- 1
  for (d in 1:J) P[[d]]  <- a[[d]] <- astar[[d]] <- matrix(NA,pjs[[d]],N+1)

  for (d in B_nD) {
    factors[[d]] = list()
    for (f in 1:(LEN[[d]] - 1)) {
      factors[[d]][[f]] = matrix(NA, DIM[[d]][[f + 1]], (N+1) * ranks[[d]])
    }
  }

  for (d in 1:J)  Y[[d]] = matrix(NA, nb_ind, N+1)

  ##########################################
  ### Determination of gMGCCA components ###
  ##########################################
  gmgcca.result = gmgccak(
    A                       = R,
    A_m                     = R_m,
    C                       = C,
    tau                     = tau,
    scheme                  = scheme,
    init                    = init,
    bias                    = bias,
    tol                     = tol,
    verbose                 = verbose,
    regularisation_matrices = regularisation_matrices,
    ranks                   = ranks,
    ncomp                   = ncomp
  )

  # Store tau, AVE_inner, crit
  AVE_inner = gmgcca.result$AVE_inner
  crit      = gmgcca.result$crit

  # Store Y, a, factors and weights
  a       = gmgcca.result$a
  astar   = gmgcca.result$a
  factors = gmgcca.result$factors
  weights = gmgcca.result$weights
  Y       = gmgcca.result$Y

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
    factors[[d]] = shave.matlist(factors[[d]], ncomp[d] * ranks[d])
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
             tau     = tau)

  class(out) = "gmgcca"
  return(out)
}
