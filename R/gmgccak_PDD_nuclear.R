# In this version, we add a penalty on the nuclear norm of the component
# To make it work, we write Y as Y = PD with D diagonal and we maximize
# the nuclear norm of P
gmgccak_PDD_nuclear <- function(A, C, tau = rep(1, length(A)), scheme = "centroid",
                                   verbose = FALSE, init = "svd", bias = TRUE,
                                   tol = 1e-08, na.rm = TRUE,
                                   ncomp = rep(1, length(A)),
                                   eta_decay = 0.9, rho_decay = 0.9,
                                   tol_in_decay = 0.9, rho = 2,
                                   penalty_coef = 1, A_m = NULL) {
  ### Utility functions
  criterion <- function() {
    cur_crit <- c()
    for (i in 1:J) {
      for (j in 1:J) {
        cur_crit <- c(cur_crit, C[i, j] * sum(diag(g(crossprod(Y[[i]], Y[[j]])))))
      }
      cur_crit <- c(cur_crit, penalty_coef * sum(svd(P[[i]])$d))
    }
    return(sum(cur_crit))
  }

  # Returns a n x ncomp matrix
  compute_dgx <- function(dg, j) {
    res <- lapply(seq(J), function(k) {
      2 * C[j, k] * Y[[k]] %*% diag(dg(diag(crossprod(Y[[j]], Y[[k]]))), nrow = ncol(Y[[k]]))
    })
    Reduce("+", res)
  }

  # A vector of differences between the supposed to be equal variables
  equality_constraints <- function() {
    vapply(seq(J), function(j) {
      norm(Y[[j]] - Z[[j]], type = "F")^2
    }, FUN.VALUE = 1.)
  }

  crit_lagrangian <- function() {
    criterion() - 1 / (2 * rho) * Reduce("+", vapply(seq(J), function(j) {
      norm(Y[[j]] - Z[[j]] + rho * Mu[[j]], type = "F")^2
    }, FUN.VALUE = 1.))
  }

  lagrange_search <- function(d, x, target, max_iter = 1000) {
    min_val <- max(0, sqrt(drop(crossprod(x)) / target) - max(d))
    max_val <- max(0, sqrt(drop(crossprod(x)) / target) - min(d))
    it <- 1
    repeat {
      val <- (min_val + max_val) / 2
      value <- drop(t(x) %*% diag(1 / (d + val)^2, nrow = length(d)) %*% x)
      if (abs(value - target) < tol | it > max_iter) {
        break
      }
      if (value < target) {
        max_val <- val
      } else {
        min_val <- val
      }
      it <- it + 1
    }
    return(val)
  }

  a_update <- function(j) {
    res <- list()
    Q <- grad + t(A_m[[j]]) %*% (Z[[j]] / rho - Mu[[j]])
    Q <- a[[j]] + (Q - t(A_m[[j]]) %*% Y[[j]] / rho) * rho / lambda[j]

    if (j %in% B_2D) {
      res$a <- apply(Q, 2, function(x) {
        x / norm(x, type = "2")
      })
    } else {
      res$factors <- factors[[j]]
      for (d in 1:(LEN[j] - 1)) {
        fac <- Reduce(khatri_rao, rev(res$factors[-d]))
        for (h in seq(ncomp[j])) {
          q <- unfold(array(Q[, h], dim = DIM[[j]][-1]), mode = d)
          res$factors[[d]][, h] <- q %*% fac[, h]
        }
        res$factors[[d]] <- apply(
          res$factors[[d]], 2, function(x) x / norm(x, type = "2")
        )
      }
      res$a <- list_khatri_rao(res$factors)
    }

    return(res)
  }

  # Write Z as a close to orthogonal matrix times a diagonal one and alternate
  # between the two
  Z_update <- function(j) {
    Q <- Y[[j]] + rho * Mu[[j]]
    SVD <- svd(P[[j]])
    U <- SVD$u
    s <- SVD$d
    V <- SVD$v

    # Update V - we minimize an upper bound
    E <- (D[[j]] %*% t(Q) %*% U %*% diag(s, nrow = ncol(Q)) -
            D[[j]]^2 %*% V %*% diag(s^2, nrow = ncol(Q))
    ) / (s[1]^2 * max(D[[j]]^2))
    SVD <- svd(E + V)
    V <- SVD$u %*% t(SVD$v)

    # Update U
    SVD <- svd(Q %*% D[[j]] %*% V %*% diag(s, nrow = ncol(Q)))
    U <- SVD$u %*% t(SVD$v)

    # Update s
    ds <- diag(t(V) %*% D[[j]]^2 %*% V)
    x <- pmax(
      0,
      rho * penalty_coef * rep(1, ncol(Q)) + diag(t(U) %*% Q %*% D[[j]] %*% V)
    )
    # Test norm of s, if it's below the number of components, the lagrange
    # multiplier associated to the norm constraint is null
    if (drop(t(x) %*% diag(1 / ds^2, nrow = ncol(Q)) %*% x) <= ncol(Q)) {
      s <- drop(diag(1 / ds, nrow = ncol(Q)) %*% x)
    } else {
      val <- lagrange_search(ds, x, ncol(Q))
      s <- drop(diag(1 / (ds + val), nrow = ncol(Q)) %*% x)
    }

    # Update D
    new_D <- drop(diag(
      1 / diag(V %*% diag(s^2, nrow = ncol(Q)) %*% t(V)),
      nrow = ncol(Q)
    ) %*% diag(V %*% diag(s, nrow = ncol(Q)) %*% t(U) %*% Q))

    return(list(
      D = diag(new_D, nrow = ncol(Q)),
      P = U %*% diag(s, nrow = ncol(Q)) %*% t(V)
    ))
  }

  if (mode(scheme) != "function") {
    if (scheme == "horst") {
      g <- function(x) x
      ctrl <- FALSE
    }
    if (scheme == "factorial") {
      g <- function(x) x^2
      ctrl <- TRUE
    }
    if (scheme == "centroid") {
      g <- function(x) abs(x)
      ctrl <- TRUE
    }
  } else {
    # check for parity of g
    g <- scheme
    ctrl <- !any(g(-5:5) != g(5:-5))
  }

  dg <- Deriv::Deriv(g)

  J <- length(A) # number of blocks
  n <- NROW(A[[1]]) # number of individuals

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

  pjs <- sapply(A_m, NCOL) # number of variables per block
  a <- factors <- P <- D <- Z <- list()
  for (j in 1:J) {
    factors[[j]] <- list()
  }

  for (j in seq(J)) {
    if (init == "svd") {
      # Initialization by SVD
      if (j %in% B_2D) {
        a[[j]] <- svd(A[[j]], nu = 0, nv = ncomp[j])$v
      } else {
        for (d in 1:(LEN[[j]] - 1)) {
          factors[[j]][[d]] <- svd(apply(A[[j]], d + 1, c), nu = 0, nv = ncomp[j])$v
        }
        a[[j]] <- list_khatri_rao(factors[[j]])
      }
    } else if (init == "random") {
      # Random Initialisation of a_j
      if (j %in% B_2D) {
        a[[j]] <- apply(
          matrix(rnorm(pjs[j] * ncomp[j]), pjs[j]), 2,
          function(x) x / norm(x, type = "2")
        )
      } else {
        for (d in 1:(LEN[[j]] - 1)) {
          factors[[j]][[d]] <- apply(
            matrix(rnorm(DIM[[j]][d + 1] * ncomp[j]), DIM[[j]][d + 1]), 2,
            function(x) x / norm(x, type = "2")
          )
        }
        a[[j]] <- list_khatri_rao(factors[[j]])
      }
    } else {
      stop_rgcca("init should be either random or by SVD.")
    }
  }


  Y <- lapply(seq(J), function(j) {
    pm(A_m[[j]], a[[j]], na.rm = na.rm)
  })

  # Find the spectral norm of t(A) %*% A
  lambda <- vapply(A_m, function(x) {
    if (nrow(x) > ncol(x)) {
      eigen(crossprod(x), symmetric = TRUE, only.values = TRUE)$values[1]
    } else {
      eigen(tcrossprod(x), symmetric = TRUE, only.values = TRUE)$values[1]
    }
  }, FUN.VALUE = 1.)

  for (j in seq(J)) {
    SVD <- svd(Y[[j]])
    P[[j]] <- SVD$u %*% t(SVD$v)
    D[[j]] <- diag(SVD$d, nrow = ncol(Y[[j]]))
    Z[[j]] <- P[[j]] %*% D[[j]]
  }
  Mu <- lapply(seq(J), function(j) matrix(0, nrow = n, ncol = ncomp[j]))

  iter <- 1
  n_iter_max <- 1000L
  iter_in_max <- 1000L
  tol_in <- 1e-3
  crit <- numeric(n_iter_max)
  crit_old <- criterion()
  eta <- eta_decay * max(equality_constraints())

  ### PDD algorithm
  repeat {
    iter_in <- 1
    a_old <- a
    Z_old <- Z
    crit_lag_old <- crit_lagrangian()
    c_lag_old <- crit_lag_old

    # Update of the primal variables
    repeat {
      for (j in 1:J) {
        dgx <- compute_dgx(dg, j)
        grad <- pm(t(A_m[[j]]), dgx, na.rm = na.rm)

        res <- a_update(j)
        a[[j]] <- res$a
        if (!(j %in% B_2D)) {
          factors[[j]] <- res$factors
        }
        Y[[j]] <- pm(A_m[[j]], a[[j]], na.rm = na.rm)

        c_lag <- crit_lagrangian()
        if (((c_lag - c_lag_old)  < 0)) {
          print(paste0("issue with update of a[[", j, "]]"))
          # browser()
        }
        c_lag_old <- c_lag
      }

      for (j in seq(J)) {
        res <- Z_update(j)
        P[[j]] <- res$P
        D[[j]] <- res$D
        Z[[j]] <- P[[j]] %*% D[[j]]

        c_lag <- crit_lagrangian()
        # if ((c_lag - c_lag_old)  < 0) {
        #   print(paste0("issue with update of Z[[", j, "]]: ", (c_lag - c_lag_old) / abs(c_lag_old)))
        #   # browser()
        # }
        c_lag_old <- c_lag
      }

      crit_lag <- crit_lagrangian()

      # stopping_crit <- max(
      #   max(unlist(a, FALSE, FALSE) - unlist(a_old, FALSE, FALSE)),
      #   max(unlist(Z, FALSE, FALSE) - unlist(Z_old, FALSE, FALSE))
      # )

      stopping_crit <- (crit_lag - crit_lag_old) / abs(crit_lag_old)

      if (stopping_crit < -tol_in) {
        # stop("Convergence issue")
        print(stopping_crit)
      }

      if (stopping_crit < tol_in | iter_in > iter_in_max) {
        # print(iter_in)
        break
      }

      a_old <- a
      Z_old <- Z
      crit_lag_old <- crit_lag
      iter_in <- iter_in + 1
    }


    # Choice between Augmented Lagrangian or penalty method
    if (sum(equality_constraints() < eta)) {
      # Update of dual variable
      Mu <- lapply(seq(J), function(j) {
        Mu[[j]] + 1 / rho * (Z[[j]] - Y[[j]])
      })
    } else {
      # Update of rho
      rho <- rho_decay * rho
    }

    crit[iter] <- criterion()

    if (verbose & (iter %% 1) == 0) {
      cat(
        " Iter: ", formatC(iter, width = 3, format = "d"),
        " Fit:", formatC(crit[iter],
                         digits = 8,
                         width = 10, format = "f"
        ),
        " Dif: ", formatC(crit[iter] - crit_old,
                          digits = 8,
                          width = 10, format = "f"
        ), "\n"
      )
    }

    stopping_criterion <- max(equality_constraints())

    # if (crit[iter] - crit_old < -tol) {
    #   # stop_rgcca("Convergence issue")
    #   print("Convergence issue")
    # }
    if ((stopping_criterion < tol) | (iter > n_iter_max)) {
      break
    }
    crit_old <- crit[iter]
    iter <- iter + 1
    eta <- eta_decay * min(eta, max(equality_constraints()))
    tol_in <- max(tol, tol_in_decay * tol_in)
  }

  for (j in seq_len(J)) {
    if (ctrl & a[[j]][1] < 0) {
      a[[j]] <- -a[[j]]
      Y[[j]] <- pm(A_m[[j]], a[[j]], na.rm = na.rm)
    }
  }

  crit <- crit[which(crit != 0)]

  if (iter > n_iter_max) {
    warning(
      "The RGCCA algorithm did not converge after", n_iter_max,
      " iterations."
    )
  }
  if (iter < n_iter_max & verbose) {
    cat(
      "The RGCCA algorithm converged to a stationary point after",
      iter - 1, "iterations \n"
    )
  }
  if (verbose) {
    plot(crit, xlab = "iteration", ylab = "criteria")
  }

  # AVEinner <- sum(C * cor(Y)^2/2)/(sum(C)/2)
  AVEinner <- NULL

  result <- list(
    Y = Y, a = a, crit = crit,
    AVE_inner = AVEinner, tau = tau
  )
  return(result)
}
