tgccak <- function(A, A_m = NULL, C, tau = rep(1, length(A)),
                   scheme = "centroid", verbose = FALSE,
                   init = "svd", bias = TRUE, tol = 1e-8,
                   regularisation_matrices, ranks = rep(1, length(A)),
                   n_run = 1, n_cores = 1, orth_modes = 1) {
  call <- list(
    A = A, A_m = A_m, C = C, scheme = scheme,
    verbose = verbose, init = init,
    bias = bias, tol = tol, ranks = ranks
  )

  if (mode(scheme) != "function") {
    if (!scheme %in% c("horst", "factorial", "centroid")) {
      stop_rgcca("Please choose scheme as 'horst', 'factorial', 'centroid'")
    }
    if (scheme == "horst") {
      g <- function(x) x
    }
    if (scheme == "factorial") {
      g <- function(x) x^2
    }
    if (scheme == "centroid") {
      g <- function(x) abs(x)
    }
  } else {
    g <- scheme
  }

  dg <- Deriv::Deriv(g)

  criterion <- function() {
    cur_crit <- 0
    for (i in seq(J)) {
      for (j in seq(J)) {
        cur_crit <- cur_crit + C[i, j] * g(crossprod(Y[[i]], Y[[j]]))
      }
    }
    return(drop(cur_crit) / (n - 1 + bias))
  }

  # Compute the gradient with respect to a[[j]] : dim = pjs[j]
  compute_dgx <- function(j) {
    res <- lapply(seq(J), function(i) {
      2 * C[j, i] * Y[[i]] %*% dg(crossprod(Y[[i]], Y[[j]]))
    })
    drop(t(A_m[[j]]) %*% Reduce("+", res))
  }

  ######################
  ### Initialization ###
  ######################

  J <- length(A)

  # List of 2D matrix and higher order tensors
  DIM <- lapply(A, dim)
  LEN <- unlist(lapply(DIM, length))
  B_nD <- which(LEN >= 4) # Store which blocks are higher order tensors
  B_3D <- which(LEN == 3) # Store which blocks are 3D
  B_2D <- which(LEN == 2) # Store which blocks are 2D
  B_0D <- which(LEN == 0) # Store which blocks are 1D (stored as 0D)

  # Convert vectors to one-column matrices
  if (length(B_0D) != 0) {
    for (i in B_0D) {
      A[[i]] <- as.matrix(A[[i]])
      DIM[[i]] <- dim(A[[i]])
    }
    B_2D <- c(B_2D, B_0D)
  }

  # Dimensions of each block
  n <- DIM[[1]][1]
  pjs <- vapply(DIM, function(x) prod(x[-1]), FUN.VALUE = 1.)
  n_iter_max <- 1000L

  # Matricization (mode-1)
  if (is.null(A_m)) {
    A_m <- lapply(1:J, function(x) matrix(as.vector(A[[x]]), nrow = n))
  }

  # Initialization of factors
  factors <- lapply(seq(J), function(j) {
    if (j %in% B_2D) {
      return(NULL)
    }
    if (init == "svd") {
      lapply(seq(2, LEN[[j]]), function(d) {
        svd(matrix(aperm(
          A[[j]], perm = c(d, seq_along(DIM[[j]])[-d])
        ), nrow = DIM[[j]][d]), nu = ranks[j])$u
      })
    } else {
      lapply(seq(2, LEN[[j]]), function(d) {
        svd(matrix(
          rnorm(DIM[[j]][d] * ranks[j], mean = 0, sd = 1), nrow = DIM[[j]][d]
        ), nu = ranks[j])$u
      })
    }
  })

  # Initialization of the regularization matrices
  M <- lapply(seq(J), function(j) {
    if (tau[j] == 1) {
      return(NULL)
    }
    # tau[j] * diag(pjs[j]) + (1 - tau[j]) * crossprod(A_m[[j]]) / (n - 1 + bias)
    tau[j] * diag(pjs[j]) + (1 - tau[j]) * crossprod(A_m[[j]])
  })

  Minv <- lapply(seq(J), function(j) {
    if (is.null(M[[j]])) {
      return(NULL)
    }
    if (!(j %in% B_2D)) {
      return(NULL)
    }
    return(solve(M[[j]]))
  })

  # Initialization of weights vectors a
  a <- lapply(seq(J), function(j) {
    if (j %in% B_2D) {
      if (init == "svd") {
        return(initsvd(A[[j]], dual = FALSE))
      } else {
        return(rnorm(n = pjs[[j]], mean = 0, sd = 1))
      }
    }
    Reduce(khatri_rao, rev(factors[[j]])) %*% rep(1, ranks[j])
  })

  # Normalize a and factors
  for (j in seq(J)) {
    if (is.null(M[[j]])) {
      a_norm <- sqrt(drop(crossprod(a[[j]])))
    } else {
      a_norm <- sqrt(drop(t(a[[j]]) %*% M[[j]] %*% a[[j]]))
    }
    a[[j]] <- a[[j]] / a_norm
    if (!(j %in% B_2D)) {
      factors[[j]][[1]] <- factors[[j]][[1]] / a_norm
    }
  }

  # Initialization of vector Y
  Y <- lapply(seq(J), function(j) {
    A_m[[j]] %*% a[[j]]
  })

  # Initialization of other parameters
  iter <- 1
  crit <- numeric()
  a_old <- a
  crit_old <- criterion()

  # TGCCA algorithm
  repeat {
    for (j in seq(J)) {
      grad_j <- compute_dgx(j)

      if (j %in% B_2D) {
        if (is.null(M[[j]])) {
          a[[j]] <- grad_j / sqrt(drop(crossprod(grad_j)))
        } else {
          grad_j <- Minv[[j]] %*% grad_j
          a[[j]] <- grad_j / sqrt(drop(t(grad_j) %*% M[[j]] %*% grad_j))
        }
        Y[[j]] <- A_m[[j]] %*% a[[j]]

      } else if (j %in% B_3D) {
        # Solve for first factor
        grad <- matrix(grad_j, nrow = DIM[[j]][2], ncol = DIM[[j]][3])

        # if (!is.null(M[[j]])) {
        #   Mx <- Minv[[j]] %*% grad_j
        #   Pi_Mx <- tcrossprod(Mx) / drop(crossprod(Mx))
        #   Pi_B <- factors[[j]][[2]] %*% solve(crossprod(factors[[j]][[2]])) %*% t(factors[[j]][[2]])
        #   # U <- t(Pi_B) %x% diag(DIM[[j]][2]) + Pi_Mx
        #   U <- (Pi_B %x% diag(DIM[[j]][2])) %*% (diag(DIM[[j]][2] * DIM[[j]][3]) - Pi_Mx)
        #
        #   tgcca_partial_update(M, grad_j, Pi_B, Pi_Mx, DIM, j)
        #
        #   y <- ginv(U) %*% (t(Pi_B) %x% diag(DIM[[j]][2])) %*% grad_j
        #   y <- y * sqrt(drop((t(grad_j) %*% Minv[[j]] %*% grad_j) / (t(y) %*% Minv[[j]] %*% y)))
        #
        #   grad <- matrix(Minv[[j]] %*% (as.vector(grad) + y), nrow = DIM[[j]][2])
        # }
        #
        # factors[[j]][[1]] <- t(solve(crossprod(factors[[j]][[2]]), t(grad %*% factors[[j]][[2]])))

        factors[[j]][[1]] <- solution(grad, factors[[j]][[2]], M[[j]])

        a[[j]] <- khatri_rao(factors[[j]][[2]], factors[[j]][[1]]) %*% rep(1, ranks[j])
        # if (is.null(M[[j]])) {
        #   a_norm <- sqrt(drop(crossprod(a[[j]])))
        # } else {
        #   a_norm <- sqrt(drop(t(a[[j]]) %*% M[[j]] %*% a[[j]]))
        # }
        # a[[j]] <- a[[j]] / a_norm
        # factors[[j]][[1]] <- factors[[j]][[1]] / a_norm
        Y[[j]] <- A_m[[j]] %*% a[[j]]

        # Solve for second factor
        idx1to2 <- as.vector(matrix(seq_along(grad_j), nrow = DIM[[j]][3], byrow = TRUE))
        # idx2to1 <- as.vector(matrix(seq_along(grad_j), nrow = DIM[[j]][2], byrow = TRUE))

        grad <- matrix(grad_j, nrow = DIM[[j]][2], ncol = DIM[[j]][3])

        # if (!is.null(M[[j]])) {
        #   # U <- svd(factors[[j]][[1]])$u
        #   # grad <- tcrossprod(U) %*% grad
        #
        #   # grad <- matrix(
        #   #   Minv[[j]][idx1to2, idx1to2] %*% as.vector(t(grad)), nrow = DIM[[j]][2]
        #   # )[idx2to1, idx2to1]
        #   grad <- matrix(
        #     Minv[[j]][idx1to2, idx1to2] %*% as.vector(t(grad)),
        #     nrow = DIM[[j]][2], byrow = TRUE
        #   )
        # }
        #
        # factors[[j]][[2]] <- t(solve(crossprod(factors[[j]][[1]]), t(t(grad) %*% factors[[j]][[1]])))

        factors[[j]][[2]] <- solution(t(grad), factors[[j]][[1]], M[[j]][idx1to2, idx1to2])

        a[[j]] <- khatri_rao(factors[[j]][[2]], factors[[j]][[1]]) %*% rep(1, ranks[j])
        # if (is.null(M[[j]])) {
        #   a_norm <- sqrt(drop(crossprod(a[[j]])))
        # } else {
        #   a_norm <- sqrt(drop(t(a[[j]]) %*% M[[j]] %*% a[[j]]))
        # }
        # a[[j]] <- a[[j]] / a_norm
        # factors[[j]][[2]] <- factors[[j]][[2]] / a_norm
        Y[[j]] <- A_m[[j]] %*% a[[j]]

      } else {
        for (d in seq(1, LEN[[j]] - 1)) {
          grad <- array(grad_j, dim = DIM[[j]][-1])
          grad <- matrix(aperm(
            grad, perm = c(d, seq_along(DIM[[j]][-1])[-d])
          ), nrow = DIM[[j]][d + 1])
          Z <- Reduce(khatri_rao, rev(factors[[j]][-d]))

          if (!is.null(M[[j]])) {
            grad <- Minv[[j]] %*% as.vector(grad)
            grad <- matrix(grad, nrow = DIM[[j]][d + 1])
          }

          factors[[j]][[d]] <- t(solve(crossprod(Z), t(grad %*% Z)))

          # Update and normalize all quantities
          a[[j]] <- Reduce(khatri_rao, rev(factors[[j]])) %*% rep(1, ranks[j])
          if (is.null(M[[j]])) {
            a_norm <- sqrt(drop(crossprod(a[[j]])))
          } else {
            a_norm <- sqrt(drop(t(a[[j]]) %*% M[[j]] %*% a[[j]]))
          }
          a[[j]] <- a[[j]] / a_norm
          factors[[j]][[d]] <- factors[[j]][[d]] / a_norm
          Y[[j]] <- A_m[[j]] %*% a[[j]]
          print(paste0("block ", j, ", mode ", d, ": ", criterion()))
        }
      }
    }

    crit[iter] <- criterion()
    if (verbose)
    {
      cat(" Iter: ", formatC(iter, width = 3, format = "d"),
          " Fit:", formatC(crit[iter], digits = 8,
                           width = 10, format = "f"),
          " Dif: ", formatC(crit[iter] - crit_old, digits = 8,
                            width = 10, format = "f"), "\n")
    }

    stopping_criteria <- c(
      drop(crossprod(unlist(a, FALSE, FALSE) - unlist(a_old, FALSE, FALSE))),
      abs(crit[iter] - crit_old) / crit[iter]
    )
    # Criterion must increase
    if (crit[iter] - crit_old < -tol) {
      stop_rgcca("Convergence error: criterion did not increase monotonously")
    }
    if (any(stopping_criteria < tol) | (iter > n_iter_max)) break

    crit_old <- crit[iter]
    a_old <- a
    iter <- iter + 1
  }

  # Final messages
  if (iter > 1000) {
    warning("The MGCCA algorithm did not converge after 1000 iterations.")
  }
  if (iter < 1000 & verbose) {
    cat(
      "The MGCCA algorithm converged to a stationary point after ",
      iter - 1, " iterations \n"
    )
  }
  if (verbose) {
    plot(crit[1:iter], xlab = "iteration", ylab = "criteria")
  }

  Y <- do.call(cbind, Y)
  AVEinner <- sum(C * cor(Y)^2 / 2) / (sum(C) / 2)

  weights <- lapply(seq(J), function(j) rep(1, ranks[j]))

  result <- list(
    Y = Y,
    a = a,
    factors = factors,
    weights = weights,
    crit = crit,
    AVE_inner = AVEinner,
    call = call,
    tau = tau
  )

  return(result)
}
