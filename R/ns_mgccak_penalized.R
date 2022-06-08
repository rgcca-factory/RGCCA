ns_mgccak_penalized <- function(A, A_m = NULL, C, tau = rep(1, length(A)), scheme = "centroid",
                       verbose = FALSE, init="svd", bias = TRUE, tol = 1e-8,
                       regularisation_matrices, ranks= rep(1, length(A)),
                       kronecker_covariance = F, penalty_coef = 0) {

  call = list(A = A, A_m = A_m, C = C, scheme = scheme, verbose = verbose, init = init,
              bias = bias, tol = tol, ranks = ranks)

  criter = function() {
    cur_crit = c()
    for (i in 1:J) {
      for (j in 1:J) {
        cur_crit = c(cur_crit, C[i, j] * sum(diag(g(crossprod(Y[[i]], Y[[j]])))))
      }
    }
    return(sum(cur_crit))
  }

  criterion = function() {
    return(criter() - penalty_coef * penalty)
  }

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
  B_0D   <- which(LEN == 0)   # Store which blocks are 1D (stored as 0D)

  # Convert vectors to one-column matrices
  if (length(B_0D) != 0) {
    for (i in B_0D) {
      A[[i]]   = as.matrix(A[[i]])
    }
  }

  # Dimensions of each block
  n   <- nrow(A[[1]])
  Y   <- matrix(0, n, J)

  # Matricization (mode-1)
  if(is.null(A_m)){
    A_m = lapply(1:J, function(x) matrix(as.vector(A[[x]]), nrow = n))
  }

  # Initialization of vector a (weight vector)
  res_init = ns_mgcca_penalized_init(A, A_m, tau = tau, ranks = ranks, init = init,
                           bias = bias, kronecker_covariance = kronecker_covariance)
  a = res_init$a; factors = res_init$factors; weights = res_init$weights; XtX = res_init$XtX; XtX_sing = res_init$XtX_sing; penalty = res_init$penalty

  # Initialization of vector Y
  for (j in 1:J) Y[, j] <- A_m[[j]] %*% a[[j]]

  # Initialize other parameters
  crit_old = criterion()
  iter     = 1
  crit     = crit_first_part = crit_second_part = numeric()
  a_old    = a

  dg = Deriv::Deriv(g)

  # MGCCA algorithm
  repeat {
    res_update = ns_mgcca_penalized_update(A, A_m, a, factors, weights, XtX, XtX_sing, Y, g, dg, C, ranks = ranks, bias = bias, penalty_coef = penalty_coef)
    a = res_update$a; factors = res_update$factors; weights = res_update$weights; Y = res_update$Y; penalty = res_update$penalty

    crit[iter] <- criterion()
    crit_first_part[iter] <- criter()
    crit_second_part[iter] <- penalty_coef * penalty

    if (verbose)
    {
      cat(" Iter: ", formatC(iter, width = 3, format = "d"),
          " Fit:", formatC(crit[iter], digits = 8,
                           width = 10, format = "f"),
          " Dif: ", formatC(crit[iter] - crit_old, digits = 8,
                            width = 10, format = "f"), "\n")
    }

    stopping_criteria = c(drop(crossprod(unlist(a, F, F) - unlist(a_old, F, F))),
                          abs(crit[iter] - crit_old))
    # Criterion must increase
    if ( crit[iter] - crit_old < -tol)
    {stop_rgcca("Convergence error: criterion did not increase monotonously")}
    if (any(stopping_criteria < tol) | (iter > 1000)) break

    crit_old = crit[iter]
    a_old <- a
    iter <- iter + 1
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
    plot(crit, xlab = "iteration", ylab = "criteria")
    plot(crit_first_part[1:iter], xlab = "iteration", ylab = "Criterion (correlations)")
    plot(crit_second_part[1:iter], xlab = "iteration", ylab = "Criterion (penalty)")
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
