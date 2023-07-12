# Try PDD to look for mixed solutions of the optimization problem where
# on one side we impose Kronecker constraints, and orthogonality constraints
# on the other side
# We add upper bound norm constraints on Y
gmgccak_PDD5 <- function(A, C, tau = rep(1, length(A)), scheme = "centroid",
                         verbose = FALSE, init = "svd", bias = TRUE,
                         tol = 1e-08, na.rm = TRUE,
                         ncomp = rep(1, length(A)),
                         eta_decay = 0.9, rho_decay = 0.9,
                         tol_in_decay = 0.9, rho = 2, A_m = NULL,
                         penalty_coef = 0, orth_Y = FALSE) {
  orth_Y = TRUE
  rho = 10

  ### Utility functions
  criterion <- function() {
    cur_crit <- c()
    for (i in seq(J)) {
      for (j in seq(J)) {
        if (orth_Y) {
          cur_crit <- c(
            cur_crit,
            C[i, j] * sum(diag(g(crossprod(Y[[i]], Z[[j]]))))
          )
        } else {
          cur_crit <- c(
            cur_crit,
            C[i, j] * sum(diag(g(crossprod(Y[[i]], A_m[[j]] %*% Z[[j]]))))
          )
        }
      }
    }
    return(sum(cur_crit))
  }

  # Returns a n x ncomp matrix
  compute_dgx <- function(dg, j) {
    if (orth_Y) {
      res <- lapply(seq(J), function(k) {
        C[j, k] * Z[[k]] %*% diag(dg(diag(crossprod(Y[[j]], Z[[k]]))), nrow = ncol(Z[[k]]))
      })
    } else {
      res <- lapply(seq(J), function(k) {
        C[j, k] * A_m[[k]] %*% Z[[k]] %*% diag(dg(diag(crossprod(Y[[j]], A_m[[k]] %*% Z[[k]]))), nrow = ncol(Z[[k]]))
      })
    }
    Reduce("+", res)
  }

  compute_dgx_Z <- function(dg, j) {
    if (orth_Y) {
      res <- lapply(seq(J), function(k) {
        C[j, k] * Y[[k]] %*% diag(dg(diag(crossprod(Z[[j]], Y[[k]]))), nrow = ncol(Y[[k]]))
      })
    } else {
      res <- lapply(seq(J), function(k) {
        C[j, k] * Y[[k]] %*% diag(dg(diag(crossprod(A_m[[j]] %*% Z[[j]], Y[[k]]))), nrow = ncol(Y[[k]]))
      })
    }
    Reduce("+", res)
  }

  # A vector of differences between the supposed to be equal variables
  equality_constraints <- function() {
    vapply(seq(J), function(j) {
      if (orth_Y) {
        norm(Y[[j]] - Z[[j]], type = "F")^2
      } else {
        norm(a[[j]] - Z[[j]], type = "F")^2
      }
    }, FUN.VALUE = 1.)
  }

  crit_lagrangian <- function() {
    criterion() - 1 / (2 * rho) * Reduce("+", vapply(seq(J), function(j) {
      if (orth_Y) {
        norm(Y[[j]] - Z[[j]] + rho * Mu[[j]], type = "F")^2
      } else {
        norm(a[[j]] - Z[[j]] + rho * Mu[[j]], type = "F")^2
      }
    }, FUN.VALUE = 1.))
  }

  a_update <- function(j) {
    grad <- compute_dgx(dg, j)

    res <- list()
    if (orth_Y) {
      Q <- 1 / lambda[j] * (
        t(A_m[[j]]) %*% (Z[[j]] - Y[[j]] + rho * (grad - Mu[[j]]))
      ) + a[[j]]
    } else {
      Q <- t(A_m[[j]]) %*% grad + (Z[[j]] / rho - Mu[[j]])
    }


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

  Z_update <- function(j) {
    grad <- compute_dgx_Z(dg, j)

    if (orth_Y) {
      # Write Z as an orthogonal matrix times
      # a diagonal one with a norm constraint
      Q <- Y[[j]] + rho * (grad + Mu[[j]])
      delta <- sqrt(diag(crossprod(Z[[j]])))
      SVD <- svd(Q %*% diag(delta, nrow = ncol(Q)))
      delta <- diag(SVD$v %*% t(SVD$u) %*% Q)
      if (norm(delta, type = "2")^2 <= (ncol(Q) * lambda[j])) {
        delta <- delta
      } else {
        delta <- delta * sqrt(lambda[j] * ncol(Q)) / norm(delta, type = "2")
      }
      SVD$u %*% t(SVD$v) %*% diag(delta, nrow = ncol(Q))
    } else {
      grad <- pm(t(A_m[[j]]), grad, na.rm = na.rm)
      Q <- grad + (a[[j]] / rho + Mu[[j]])
      SVD <- svd(Q)
      SVD$u %*% t(SVD$v)
    }
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
  a <- factors <- Z <- list()
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

  if (orth_Y) {
    # Z <- lapply(Y, function(x) {
    #   SVD <- svd(x)
    #   SVD$u %*% t(SVD$v) %*% diag(SVD$d, nrow = ncol(x))
    # })
    Z <- lapply(Y, function(x) {
      SVD <- svd(matrix(rnorm(prod(dim(x))), nrow = nrow(x)))
      SVD$u %*% t(SVD$v) %*% diag(sqrt(diag(crossprod(x))), nrow = ncol(x))
    })
  } else {
    Z <- lapply(a, function(x) {
      SVD <- svd(x)
      SVD$u %*% t(SVD$v)
    })
  }

  if (orth_Y) {
    # Mu <- lapply(seq(J), function(j) matrix(0, nrow = n, ncol = ncomp[j]))
    # # Standard deviation such that E[Mu[[j]][, h]] = lambda[j]
    # Mu <- lapply(seq(J), function(j) matrix(
    #   rnorm(n * ncomp[j], sd = sqrt(lambda[j] / n)), nrow = n, ncol = ncomp[j]
    # ))
    # Standard deviation such that E[Mu[[j]][, h]] = lambda[j]
    Mu <- lapply(seq(J), function(j) vapply(seq(ncomp[j]), function(h) {
      rnorm(n, sd = norm(Y[[j]][, h], type = "2") / sqrt(n))
    }, FUN.VALUE = double(n)))
  } else {
    Mu <- lapply(seq(J), function(j) matrix(0, nrow = pjs[j], ncol = ncomp[j]))
  }

  iter <- 1
  n_iter_max <- 1000L
  iter_in_max <- 1000L
  tol_in <- 1e-3
  crit <- numeric(n_iter_max)
  crit_old <- criterion()
  eta <- eta_decay * max(equality_constraints())

  a0 <- a
  Y0 <- Y
  Z0 <- Z

  counts <- c(0, 0)

  ### PDD algorithm
  repeat {
    iter_in <- 1
    a_old <- a
    Z_old <- Z
    crit_lag_old <- crit_lagrangian()

    # Update of the primal variables
    repeat {
      res <- lapply(seq(J), a_update)
      a <- lapply(res, '[[', "a")

      for (j in c(B_3D, B_nD)) {
        factors[[j]] <- res[[j]]$factors
      }

      Y <- lapply(seq(J), function(j) {
        Y[[j]] <- pm(A_m[[j]], a[[j]], na.rm = na.rm)
      })

      Z <- lapply(seq(J), Z_update)

      crit_lag <- crit_lagrangian()

      # stopping_crit <- max(
      #   max(unlist(a, FALSE, FALSE) - unlist(a_old, FALSE, FALSE)),
      #   max(unlist(Z, FALSE, FALSE) - unlist(Z_old, FALSE, FALSE))
      # )

      stopping_crits <- c((crit_lag - crit_lag_old) / abs(crit_lag_old), max(
        max(abs(unlist(a) - unlist(a_old))), max(abs(unlist(Z) - unlist(Z_old)))
      ))
      # print(paste0(iter_in, ": ", stopping_crit))

      if (stopping_crits[1] < -tol_in) {
        # stop("Convergence issue")
        print(stopping_crit)
      }

      if (any(stopping_crits < tol_in) | iter_in > iter_in_max) {
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
      counts[1] <- counts[1] + 1
      # Update of dual variable
      if (orth_Y) {
        Mu <- lapply(seq(J), function(j) {
          Mu[[j]] + 1 / rho * (Z[[j]] - Y[[j]])
        })
      } else {
        Mu <- lapply(seq(J), function(j) {
          Mu[[j]] + 1 / rho * (Z[[j]] - a[[j]])
        })
      }
    } else {
      counts[2] <- counts[2] + 1
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

  print(paste0("Number of dual updates: ", counts[1]))
  print(paste0("Number of rho updates: ", counts[2]))

  result <- list(
    Y = Y, a = a, crit = crit,
    AVE_inner = AVEinner, tau = tau
  )
  return(result)
}
