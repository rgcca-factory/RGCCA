ns_mgccak <- function (A, A_m = NULL, C, tau = rep(1, length(A)), scheme = "centroid",
                    verbose = FALSE, init="svd", bias = TRUE, tol = 1e-8,
                    regularisation_matrices, ranks= rep(1, length(A)),
                    kronecker_covariance = F, n_run = 1, n_cores = 1) {

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

  ### Get useful constants
  J      <- length(A)
  DIM    <- lapply(A, dim)
  LEN    <- unlist(lapply(DIM, length))
  B_2D   <- which(LEN <= 2)   # Store which blocks are 2D
  B_0D   <- which(LEN == 0)   # Store which blocks are 1D (stored as 0D)

  # Dimensions of each block
  pjs <- sapply(DIM, function(x) prod(x[-1]))
  n   <- DIM[[1]][1]

  # Convert vectors to one-column matrices
  if (length(B_0D) != 0) {
    for (i in B_0D) {
      A[[i]]   = as.matrix(A[[i]])
    }
  }

  # Dimensions of each block
  n   <- nrow(A[[1]])

  # Matricization (mode-1)
  if(is.null(A_m)){
    A_m = lapply(1:J, function(x) matrix(as.vector(A[[x]]), nrow = n))
  }

  if (kronecker_covariance) {
    XtX = lapply(1:J, function(j) {
      if (j %in% B_2D) {
        return((1 - tau[j]) * crossprod(A_m[[j]]) / (n - 1 + bias) + tau[j] * diag(pjs[j]))
      }
      fac = estimate_kronecker_covariance(A[[j]])
      to  = tau[j] ^ (1 / length(fac))
      fac = lapply(fac, function(x) {
        (1 - to) * x + to * diag(nrow(x))
      })
      return(Reduce("%x%", rev(fac)))
    })
  } else {
    XtX = lapply(1:J, function(j) {
      if (!(j %in% B_2D)) return(NULL)
      (1 - tau[j]) * crossprod(A_m[[j]]) / (n - 1 + bias) + tau[j] * diag(pjs[j])
    })
  }

  XtX_sing <- lapply(seq(J), function(j) {
    if (j %in% B_2D) return(NULL)
    if (pjs[j] > n) {
      return(
        eigen(
          tcrossprod(A_m[[j]]), symmetric = TRUE, only.values = TRUE
        )$values[1] * (1 - tau[j]) / (n - 1 + bias) + tau[j]
      )
    }
    return(
      eigen(
        crossprod(A_m[[j]]), symmetric = TRUE, only.values = TRUE
      )$values[1] * (1 - tau[j]) / (n - 1 + bias) + tau[j]
    )
  })

  models <- pbapply::pblapply(seq(n_run), function(run_number) {
    init <- ifelse(run_number == 1, "svd", "random")

    core_ns_mgcca(A, A_m, XtX, XtX_sing, DIM, LEN, B_2D, B_3D, B_nD, init, g,
                  verbose, C, tol, n_iter_max, bias, ranks, tau)

  }, cl = n_cores)

  best_model <- which.max(lapply(models, function(m) m$crit[length(m$crit)]))
  a <- models[[best_model]]$a
  Y <- models[[best_model]]$Y
  factors <- models[[best_model]]$factors
  weights <- models[[best_model]]$weights
  crit <- models[[best_model]]$crit
  iter <- length(unlist(crit))

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
