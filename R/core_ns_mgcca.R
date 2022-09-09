core_ns_mgcca <- function(A, A_m, XtX, XtX_sing, DIM, LEN, B_2D, B_3D, B_nD,
                          init, g, verbose, C, tol, n_iter_max, bias, ranks,
                          tau) {
  J      <- length(A)
  n <- NROW(A[[1]])

  # Initialization of vector a (weight vector)
  res_init = ns_mgcca_init(A, A_m, XtX, XtX_sing, tau = tau, ranks = ranks,
                           init = init, bias = bias)
  a = res_init$a; factors = res_init$factors; weights = res_init$weights

  # Initialization of vector Y
  Y   <- matrix(0, n, J)
  for (j in 1:J) Y[, j] <- A_m[[j]] %*% a[[j]]

  # Initialize other parameters
  crit_old = sum(C * g(cov2(Y, bias = bias)))
  iter     = 1
  crit     = numeric()
  a_old    = a

  dg = Deriv::Deriv(g)

  # MGCCA algorithm
  repeat {
    res_update = ns_mgcca_update(A, A_m, a, factors, weights, XtX, XtX_sing, Y, g, dg, C, ranks = ranks, bias = bias, tau = tau)
    a = res_update$a; factors = res_update$factors; weights = res_update$weights; Y = res_update$Y

    crit[iter] <- sum(C*g(cov2(Y, bias = bias)))

    if (verbose)
    {
      cat(" Iter: ", formatC(iter, width = 3, format = "d"),
          " Fit:", formatC(crit[iter], digits = 8,
                           width = 10, format = "f"),
          " Dif: ", formatC(crit[iter] - crit_old, digits = 8,
                            width = 10, format = "f"), "\n")
    }

    stopping_criteria = c(
      drop(crossprod(unlist(a, F, F) - unlist(a_old, F, F))),
      crit[iter] - crit_old
    )

    # Criterion must increase
    if ( crit[iter] - crit_old < -tol)
    {stop_rgcca("Convergence error: criterion did not increase monotonously")}
    if (any(stopping_criteria < tol) | (iter > 1000)) break

    crit_old = crit[iter]
    a_old <- a
    iter <- iter + 1
  }

  return(list(a = a, Y = Y, weights = weights, factors = factors, crit = crit))
}
