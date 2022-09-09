mgccak <- function (A, A_m = NULL, C, tau = rep(1, length(A)), scheme = "centroid",
                    verbose = FALSE, init="svd", bias = TRUE, tol = 1e-8,
                    regularisation_matrices, ranks= rep(1, length(A)),
                    n_run = 1, n_cores = 1) {

  call = list(A = A, A_m = A_m, C = C, scheme = scheme, verbose = verbose, init = init,
              bias = bias, tol = tol, ranks = ranks)

  if(mode(scheme) != "function")
  {
    if(!scheme %in% c("horst", "factorial", "centroid"))
    {stop_rgcca("Please choose scheme as 'horst', 'factorial', 'centroid'")}
    if(scheme == "horst"){ g <- function(x) x}
    if(scheme == "factorial"){ g <- function(x)  x^2}
    if(scheme == "centroid"){g <- function(x) abs(x)}
  }
  else g <- scheme

  ######################
  ### Initialization ###
  ######################

  J      <- length(A)

  # List of 2D matrix and higher order tensors
  DIM    <- lapply(A, dim)
  LEN    <- unlist(lapply(DIM, length))
  B_nD   <- which(LEN >= 4)   # Store which blocks are higher order tensors
  B_3D   <- which(LEN == 3)   # Store which blocks are 3D
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
  n   <- DIM[[1]][1]
  n_iter_max <- 1000L

  # Matricization (mode-1)
  if(is.null(A_m)){
    A_m = lapply(1:J, function(x) matrix(as.vector(A[[x]]), nrow = n))
  }

  P <- M_inv_sqrt <- list()

  # Determination of the regularization matrix
  for (j in 1:J){
    if (j > length(regularisation_matrices)) {
      reg_matrices = NULL
    } else {
      reg_matrices = regularisation_matrices[[j]]
    }
    reg_matrices = parse_regularisation_matrices(
      reg_matrices = reg_matrices,
      tau          = tau[j],
      A            = A[[j]],
      DIM          = DIM[[j]],
      j            = j,
      bias         = bias
    )
    P[[j]]          = reg_matrices$P
    M_inv_sqrt[[j]] = reg_matrices$M_inv_sqrt
    tau[j]          = reg_matrices$tau
  }

  models <- pbapply::pblapply(seq(n_run), function(run_number) {
    init <- ifelse(run_number == 1, "svd", "random")

    core_mgcca(A, P, DIM, LEN, B_2D, B_3D, B_nD, init, g, verbose, C,
               tol, n_iter_max, bias, ranks)

  }, cl = n_cores)

  best_model <- which.max(lapply(models, function(m) m$crit[length(m$crit)]))
  a <- models[[best_model]]$a
  Y <- models[[best_model]]$Y
  factors <- models[[best_model]]$factors
  weights <- models[[best_model]]$weights
  crit <- models[[best_model]]$crit
  iter <- length(unlist(crit))

  # Inverse change of variables if needed
  if (length(M_inv_sqrt) > 0)  { # If no regularization matrix, list is empty
    for (j in 1:J) {
      if (j <= length(M_inv_sqrt) && !is.null(M_inv_sqrt[[j]])) {
        if (j %in% B_2D) {
          a[[j]] = M_inv_sqrt[[j]] %*% a[[j]]
        } else {
          for (d in 1:(LEN[[j]] - 1)) {
            factors[[j]][[d]] = M_inv_sqrt[[j]][[d]] %*% factors[[j]][[d]]
          }
          a[[j]] = weighted_kron_sum(factors[[j]], weights[[j]])
        }
      }
    }
  }

  # Final messages
  if (iter > 1000) {
    warning("The MGCCA algorithm did not converge after 1000 iterations.")
  }
  if (iter < 1000 & verbose) {
    cat("The MGCCA algorithm converged to a stationary point after ",
        iter-1, " iterations \n")
  }
  if (verbose) {
    plot(crit[1:iter], xlab = "iteration", ylab = "criteria")
  }

  AVEinner <- sum(C * cor(Y)^2/2)/(sum(C)/2)

  result <- list(Y         = Y,
                 a         = a,
                 factors   = factors,
                 weights   = weights,
                 crit      = crit,
                 AVE_inner = AVEinner,
                 call      = call,
                 tau       = tau)

  return(result)
}
