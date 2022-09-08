core_nn_gmgcca <- function(A, A_m, ncomp, init, g, penalty_coef, A_sing,
                           verbose, C, orth_Y, tol, n_iter_max) {
  ### Utility functions
  penalty = function() {
    cur_penalty = c()
    for (j in 1:J) {
      if (orth_Y) {
        tmp = crossprod(Y[[j]])^2
      } else {
        tmp = crossprod(a[[j]])^2 * (g(A_sing[[j]]^2) / ncomp[j])
      }
      diag(tmp) = 0
      cur_penalty = c(cur_penalty, sum(tmp))
    }
    return(sum(cur_penalty))
  }

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
    return(criter() - penalty_coef * penalty())
  }

  alpha_i_k = function(i, k) {
    eigen(crossprod(a[[i]][, -k]), symmetric = TRUE, only.values = TRUE)$values[1]
  }

  compute_dgx = function(dg, j, k) {
    sapply(1:J, function(i)
      C[j, i] * Y[[i]][, k] * dg(drop(crossprod(Y[[j]][, k], Y[[i]][, k])))
    )
  }

  # Initialization
  J <- length(A) # number of blocks
  pjs <- sapply(A_m, NCOL) # number of variables per block
  Z <- matrix(0, NROW(A[[1]]), J)

  DIM    <- lapply(A, dim)
  LEN    <- unlist(lapply(DIM, length))
  B_2D   <- which(LEN == 2)   # Store which blocks are 2D

  a <- factors <- weights <- list()
  for (j in 1:J) {
    factors[[j]] <- list()
  }

  # Initialization of vector a (weight vector)
  for (j in 1:J) {
    if (init == "svd") {
      # SVD Initialization of a_j
      if (j %in% B_2D) {
        SVD = svd(A[[j]], nu = 0, nv = ncomp[j])
        a[[j]] <- abs(SVD$v)
        a[[j]] <- apply(a[[j]], 2, function(x) x / norm(x, type = "2"))
      } else {
        for (d in 1:(LEN[[j]] - 1)) {
          SVD = svd(apply(A[[j]], d + 1, c), nu = 0, nv = ncomp[j])
          factors[[j]][[d]] <- abs(SVD$v)
          factors[[j]][[d]] <- apply(factors[[j]][[d]], 2, function(x) x / norm(x, type = "2"))
        }
        a[[j]]       <- list_khatri_rao(factors[[j]])
      }
    } else if (init == "random") {
      # Random Initialisation of a_j
      if (j %in% B_2D) {
        A_random <- matrix(rnorm(n = pjs[[j]], mean = 0, sd = 1), nrow = pjs[j], ncol = ncomp[j])
        a[[j]] <- abs(A_random)
        a[[j]] <- apply(a[[j]], 2, function(x) x / norm(x, type = "2"))
      } else {
        A_random <- array(rnorm(n = pjs[[j]], mean = 0, sd = 1), dim = DIM[[j]][-1])
        for (d in 1:(LEN[[j]] - 1)) {
          SVD = svd(apply(A_random, d, c), nu = 0, nv = ncomp[j])
          factors[[j]][[d]] <- abs(SVD$v)
          factors[[j]][[d]] <- apply(factors[[j]][[d]], 2, function(x) x / norm(x, type = "2"))
        }
        a[[j]]       <- list_khatri_rao(factors[[j]])
      }
    } else {
      stop_rgcca("init should be either random or by SVD.")
    }
  }
  # Initialization of vector Y
  Y = lapply(1:J, function(j) A_m[[j]] %*% a[[j]])

  # Initialize other parameters
  iter <- 1
  n_iter_max <- 1000L
  crit <- numeric(n_iter_max)
  crit_first_part = numeric(n_iter_max)
  crit_second_part = numeric(n_iter_max)
  crit_old <- criterion()
  a_old = a
  dg = Deriv::Deriv(g)

  # MGCCA algorithm
  repeat {
    for (j in 1:J) {
      for (k in 1:ncomp[j]) {
        dgx      = compute_dgx(dg, j, k)
        Z[, j]   = apply(dgx, 1, sum)

        if (j %in% B_2D) {
          if (orth_Y) {
            Q <- t(A[[j]]) %*% Z[, j]
            alpha <- A_sing[[j]]^4 * eigen(crossprod(a[[j]][, -k]), symmetric = TRUE, only.values = TRUE)$values[1]
            a[[j]][, k] <- Q - 2 * penalty_coef * (t(A[[j]]) %*% (Y[[j]][, -k] %*% (t(Y[[j]][, -k]) %*% Y[[j]][, k])) - alpha * a[[j]][, k])
            a[[j]][, k] <- pmax(0, a[[j]][, k])
            a[[j]][, k] <- a[[j]][, k] / norm(a[[j]][, k], type = "2")
            Y[[j]][, k] = A[[j]] %*% a[[j]][, k]
          } else {
            Q <- t(A[[j]]) %*% Z[, j]
            alpha <- alpha_i_k(j, k)
            a[[j]][, k] <- Q - 2 * penalty_coef * (a[[j]][, -k] %*% crossprod(a[[j]][, -k], a[[j]][, k]) - alpha * a[[j]][, k]) * (g(A_sing[[j]]^2) / ncomp[j])
            a[[j]][, k] <- pmax(0, a[[j]][, k])
            a[[j]][, k] <- a[[j]][, k] / norm(a[[j]][, k], type = "2")
            Y[[j]][, k] = A[[j]] %*% a[[j]][, k]
          }
        } else {
          # Tensor case
          Q <- t(A_m[[j]]) %*% Z[, j]
          if (orth_Y) {
            alpha <- A_sing[[j]]^4 * eigen(crossprod(a[[j]][, -k]), symmetric = TRUE, only.values = TRUE)$values[1]
            Q <- Q - 2 * penalty_coef * (t(A_m[[j]]) %*% (Y[[j]][, -k] %*% (t(Y[[j]][, -k]) %*% Y[[j]][, k])) - alpha * a[[j]][, k])
          } else {
            alpha <- alpha_i_k(j, k)
            Q <- Q - 2 * penalty_coef * (a[[j]][, -k] %*% crossprod(a[[j]][, -k], a[[j]][, k]) - alpha * a[[j]][, k]) * (g(A_sing[[j]]^2) / ncomp[j])
          }
          Q <- matrix(Q, nrow = DIM[[j]][2])

          # Update first factor, while second fixed
          factors[[j]][[1]][, k] <- pmax(0, Q %*% factors[[j]][[2]][, k])
          factors[[j]][[1]][, k] <- factors[[j]][[1]][, k] / norm(factors[[j]][[1]][, k], type = "2")

          # Update second factor, while first fixed
          factors[[j]][[2]][, k] <- pmax(0, t(Q) %*% factors[[j]][[1]][, k])
          factors[[j]][[2]][, k] <- factors[[j]][[2]][, k] / norm(factors[[j]][[2]][, k], type = "2")

          a[[j]] <- list_khatri_rao(factors[[j]])
          Y[[j]][, k] = A_m[[j]] %*% a[[j]][, k]
        }
      }
    }

    # Store previous criterion
    crit[iter] <- criterion()
    crit_first_part[iter] <- criter()
    crit_second_part[iter] <- penalty_coef * penalty()

    if (verbose & (iter%%1) == 0){
      cat(" Iter: ", formatC(iter, width = 3, format = "d"),
          " Fit:", formatC(crit[iter], digits = 8,
                           width = 10, format = "f"),
          " Dif: ", formatC(crit[iter] - crit_old, digits = 8,
                            width = 10, format = "f"), "\n")
    }

    stopping_criteria = c(drop(crossprod(unlist(a, F, F) - unlist(a_old, F, F))),
                          abs((crit[iter] - crit_old) / crit[iter])
    )

    if (crit[iter] - crit_old < -tol) stop_rgcca("Convergence issue")
    if (any(stopping_criteria < tol) | (iter > n_iter_max)) {
      break
    }
    crit_old = crit[iter]
    a_old <- a
    iter <- iter + 1
  }

  return(list(a = a, crit = crit, Y = Y, factors = factors))
}
