gmgccak_penalized <- function(A, A_m = NULL, C, tau = rep(1, length(A)), scheme = "centroid",
                     verbose = FALSE, init="svd", bias = TRUE, tol = 1e-8,
                     ncomp = rep(1, length(A)), penalty_coef = 1, na.rm = TRUE,
                     orth_Y = TRUE) {

  call = list()

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

  if (mode(scheme) != "function") {
    if (scheme == "horst") {g <- function(x) x ; ctrl = FALSE}
    if (scheme == "factorial") {g <- function(x)  x^2 ; ctrl = TRUE}
    if (scheme == "centroid") {g <- function(x) abs(x) ; ctrl = TRUE}
  } else {
    # check for parity of g
    g <- scheme ; ctrl = !any(g(-5:5) != g(5:-5))
  }

  dg = Deriv::Deriv(g)

  J <- length(A) # number of blocks
  pjs <- sapply(A_m, NCOL) # number of variables per block
  Z <- matrix(0,NROW(A[[1]]),J)

  A_sing <- lapply(A_m, function(x) svd(x)$d[1])

  # List of 2D matrix and higher order tensors
  DIM    <- lapply(A, dim)
  LEN    <- unlist(lapply(DIM, length))
  B_2D   <- which(LEN == 2)   # Store which blocks are 2D
  B_0D   <- which(LEN == 0)   # Store which blocks are 1D (stored as 0D)

  # Convert vectors to one-column matrices
  if (length(B_0D) != 0) {
    for (i in B_0D) {
      A[[i]]   = as.matrix(A[[i]])
      DIM[[i]] = dim(A[[i]])
    }
    B_2D = c(B_2D, B_0D)
  }

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
        a[[j]] <- SVD$v
      } else {
        for (d in 1:(LEN[[j]] - 1)) {
          SVD = svd(apply(A[[j]], d + 1, c), nu = 0, nv = ncomp[j])
          factors[[j]][[d]] <- SVD$v
        }
        a[[j]]       <- list_khatri_rao(factors[[j]])
      }
    } else if (init == "random") {
      # Random Initialisation of a_j
      A_random <- array(rnorm(n = pjs[[j]], mean = 0, sd = 1), dim = DIM[[j]][-1])
      if (j %in% B_2D) {
        SVD = svd(A_random, nu = 0, nv = ncomp[j])
        a[[j]] <- SVD$v
      } else {
        for (d in 1:(LEN[[j]] - 1)) {
          SVD = svd(apply(A_random, d, c), nu = 0, nv = ncomp[j])
          factors[[j]][[d]] <- SVD$v
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
            # a[[j]][, k] <- pmax(0, a[[j]][, k])
            a[[j]][, k] <- a[[j]][, k] / norm(a[[j]][, k], type = "2")
            Y[[j]][, k] = A[[j]] %*% a[[j]][, k]
          } else {
            Q <- t(A[[j]]) %*% Z[, j]
            alpha <- alpha_i_k(j, k)
            a[[j]][, k] <- Q - 2 * penalty_coef * (a[[j]][, -k] %*% crossprod(a[[j]][, -k], a[[j]][, k]) - alpha * a[[j]][, k]) * (g(A_sing[[j]]^2) / ncomp[j])
            # a[[j]][, k] <- pmax(0, a[[j]][, k])
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
          SVD <- svd(Q, nu = 1, nv = 1)
          factors[[j]][[1]][, k] <- SVD$u
          factors[[j]][[2]][, k] <- SVD$v
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

  # AVEinner <- sum(C * cor(Y)^2/2)/(sum(C)/2)
  AVEinner = NULL

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
