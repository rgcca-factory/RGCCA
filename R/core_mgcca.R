core_mgcca <- function(A, P, DIM, LEN, B_2D, B_3D, B_nD, init, g, verbose, C,
                       tol, n_iter_max, bias, ranks) {
  J <- length(P)
  pjs <- sapply(DIM, function(x) prod(x[-1]))
  n   <- DIM[[1]][1]

  a <- factors <- weights <- list()
  for (j in 1:J) {
    factors[[j]] <- list()
  }

  # Initialization of vector a (weight vector)
  for (j in 1:J) {
    if (is.list(init)) {
      if (j %in% B_nD) {
        factors[[j]] <- init$factors[[j]]
        weights[[j]] <- init$weights[[j]]
        a[[j]]       <- weighted_kron_sum(factors[[j]], weights[[j]])
      } else {
        a[[j]] <- init$a[[j]][[1]]
      }
    } else if (init=="svd") {
      # SVD Initialization of a_j
      if (j %in% B_2D) {
        a[[j]] <- initsvd(A[[j]], dual = FALSE)
      } else {
        for (d in 1:(LEN[[j]] - 1)) {
          factors[[j]][[d]] <- svd(apply(A[[j]], d+1, c), nu=0, nv=ranks[[j]])$v
        }
        weights[[j]] <- rep(1 / sqrt(ranks[j]), ranks[j])
        a[[j]]       <- weighted_kron_sum(factors[[j]], weights[[j]])
      }
    } else if (init == "random") {
      # Random Initialisation of a_j
      A_random <- array(rnorm(n = pjs[[j]], mean = 0, sd = 1), dim = DIM[[j]][-1])
      if (j %in% B_2D) {
        a[[j]] <- matrix(A_random / sqrt(drop(crossprod(A_random))))
      } else {
        for (d in 1:(LEN[[j]] - 1)) {
          factors[[j]][[d]] <- svd(apply(A_random, d, c), nu=0, nv=ranks[[j]])$v
        }
        weights[[j]] <- rep(1 / sqrt(ranks[j]), ranks[j])
        a[[j]]       <- weighted_kron_sum(factors[[j]], weights[[j]])
      }
    } else {
      stop_rgcca("init should be either random or by SVD.")
    }
  }
  # Initialization of vector Y
  Y <- matrix(0, n, J)
  for (j in 1:J) Y[, j] <- P[[j]] %*% a[[j]]

  # Initialize other parameters
  crit_old = sum(C * g(cov2(Y, bias = bias)))
  iter     = 1
  crit     = numeric()
  Z        = matrix(0, n, J)
  a_old    = a

  dg = Deriv::Deriv(g)

  # MGCCA algorithm
  repeat {
    for (j in 1:J){
      # Apply the derivative on the current variables
      dgx    = dg(cov2(Y[, j], Y, bias = bias))
      dgx    = matrix(rep(dgx, n), n, J, byrow = TRUE)
      Z[, j] = rowSums(matrix(rep(C[j, ], n), n, J, byrow = TRUE) * dgx * Y)

      if (j %in% B_3D) { # 3D Tensors
        Q                 = matrix(t(Z[, j]) %*% P[[j]], nrow = DIM[[j]][3],
                                   ncol = DIM[[j]][2], byrow = T)
        SVD               = svd(x = Q, nu = ranks[[j]], nv = ranks[[j]])
        factors[[j]][[1]] = SVD$v
        factors[[j]][[2]] = SVD$u
        weights[[j]]      = SVD$d[1:ranks[j]] / sqrt(sum(SVD$d[1:ranks[j]] ^ 2))
        a[[j]]            = weighted_kron_sum(factors[[j]], weights[[j]])
        Y[, j]            = P[[j]] %*% a[[j]]

      } else if (j %in% B_nD) { # higher order Tensors
        for (d in 1:(LEN[[j]] - 1)) {
          Q                 = array(t(P[[j]]) %*% Z[, j], dim = DIM[[j]][-1])
          Q                 = unfold(Q, mode = d)
          other_factors     = list_khatri_rao(factors[[j]][-d]) %*% diag(weights[[j]])

          if (d == 1) {
            SVD               = svd(x = Q %*% other_factors, nu = ranks[j],
                                    nv = ranks[j])
            factors[[j]][[d]] = SVD$u %*% t(SVD$v)
          } else {
            factors[[j]][[d]] = apply(Q %*% other_factors, 2, function(x) x / norm(x, type = "2"))
          }

          a[[j]]            = weighted_kron_sum(factors[[j]], weights[[j]])
          Y[, j]            = P[[j]] %*% a[[j]]

          dgx              = dg(cov2(Y[, j], Y, bias = bias))
          dgx              = matrix(rep(dgx, n), n, J, byrow = TRUE)
          Z[, j]           = rowSums(
            matrix(rep(C[j, ], n), n, J, byrow = TRUE) * dgx * Y)
        }

        tmp          = t(Z[, j]) %*% P[[j]] %*% list_khatri_rao(factors[[j]])
        weights[[j]] = drop(tmp) / norm(drop(tmp), type = "2")

        a[[j]]       = weighted_kron_sum(factors[[j]], weights[[j]])
        Y[, j]       = P[[j]] %*% a[[j]]

      } else { # Matrices
        Q      = t(P[[j]]) %*% Z[,j]
        a[[j]] = Q / norm(Q, type = "2")
        Y[, j] = P[[j]] %*% a[[j]]
      }
    }
    # Store previous criterion
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
      abs(crit[iter] - crit_old) / crit[iter]
    )
    # Criterion must increase
    if ( crit[iter] - crit_old < -tol)
    {stop_rgcca("Convergence error: criterion did not increase monotonously")}
    if (any(stopping_criteria < tol) | (iter > n_iter_max)) break

    crit_old = crit[iter]
    a_old <- a
    iter <- iter + 1
  }

  return(list(a = a, factors = factors, weights = weights, Y = Y, crit = crit))
}
