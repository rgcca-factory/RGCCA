gmgccak_aux_Y <- function(
  A, A_m, C, tau, scheme, init, bias, tol, verbose, ranks, ncomp, na.rm = TRUE
) {
  criterion <- function() {
    cur_crit <- 0
    for (i in seq(J)) {
      cur_crit <- cur_crit + sum(diag(g(crossprod(U[[i]], Y[[i]]))))
      for (j in seq(J)) {
        cur_crit <- cur_crit + C[i, j] * sum(diag(g(crossprod(U[[i]], U[[j]]))))
      }
    }
    return(cur_crit / (n - 1 + bias))
  }

  # Returns a n x ncomp matrix
  compute_dgx <- function(j) {
    U[[j]] %*% diag(dg(diag(crossprod(Y[[j]], U[[j]]))))
  }

  compute_dgU <- function(j) {
    res <- lapply(seq_along(U), function(i) {
      C[i, j] * U[[i]] %*% diag(dg(diag(crossprod(U[[i]], U[[j]]))), nrow = ncol(U[[i]]))
    })
    Reduce("+", res) + Y[[j]] %*% diag(dg(diag(crossprod(Y[[j]], U[[j]]))))
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
  nc <- max(ncomp)
  pjs <- sapply(A_m, NCOL) # number of variables per block

  DIM    <- lapply(A, dim)
  LEN    <- unlist(lapply(DIM, length))
  B_nD   <- which(LEN >= 4)   # Store which blocks are higher order tensors
  B_3D   <- which(LEN == 3)   # Store which blocks are 3D
  B_2D   <- which(LEN < 3)    # Store which blocks are 2D

  if (length(B_nD) > 0) {
    stop_rgcca("not implemented.")
  }

  # Compute projection matrices
  Pi <- lapply(A_m, function(X) {
    tcrossprod(qr.Q(qr(X)))
  })

  ### Initialization ###
  weights <- lapply(seq(J), function(j) {
    rep(1 / ranks[j]^2, ranks[j] * ncomp[j])
  })

  if (init == "svd") {
    factors <- lapply(seq(J), function(j) {
      if (j %in% c(B_3D, B_nD)) {
        return(lapply(seq_len(LEN[[j]] - 1), function(d) {
          svd(apply(A[[j]], d + 1, c), nu = 0, nv = ranks[j] * ncomp[j])$v
        }))
      } else {
        return(NULL)
      }
    })

    a <- lapply(seq(J), function(j) {
      if (j %in% c(B_3D, B_nD)) {
        vapply(seq(ncomp[j]), function(k) {
          idx <- seq(1 + (k - 1) * ranks[j], k * ranks[j])
          fac <- lapply(factors[[j]], "[", , idx, drop = FALSE)
          weighted_kron_sum(fac, weights[[j]][idx])
        }, FUN.VALUE = double(pjs[j]))
      } else {
        svd(A_m[[j]], nu = 0, nv = ncomp[j])$v
      }
    })
  } else {
    stop_rgcca("not implemented.")
  }

  Y <- lapply(1:J, function(j) pm(A_m[[j]], a[[j]], na.rm = na.rm))

  U <- lapply(seq_along(A), function(j) {
    svd(Pi[[j]] %*% matrix(rnorm(n * ncomp[j]), n), nu = ncomp[j])$u
  })

  iter <- 1
  n_iter_max <- 1000L
  crit <- numeric(n_iter_max)
  crit_old <- criterion()
  a_old <- a

  repeat {
    # Update U
    for (j in seq_along(U)) {
      grad <- compute_dgU(j)
      SVD <- svd(Pi[[j]] %*% grad, nv = nc, nu = nc)
      U[[j]] <- SVD$u %*% t(SVD$v)
    }

    # Update a, factors, and weights
    for (j in B_2D) {
      grad <- pm(t(A_m[[j]]), compute_dgx(j), na.rm = na.rm)
      a[[j]] <- apply(grad, 2, function(x) {
        x / norm(x, type = "2")
      })
    }

    for (j in B_3D) {
      grad <- pm(t(A_m[[j]]), compute_dgx(j), na.rm = na.rm)
      for (k in seq(ncomp[j])) {
        grad_k <- matrix(
          grad[, k], nrow = DIM[[j]][3], ncol = DIM[[j]][2], byrow = TRUE
        )
        idx <- seq(1 + (k - 1) * ranks[j], k * ranks[j])
        SVD <- svd(x = grad_k, nu = ranks[[j]], nv = ranks[[j]])
        factors[[j]][[1]][, idx] <- SVD$v
        factors[[j]][[2]][, idx] <- SVD$u
        weights[[j]][idx] <- SVD$d[seq(ranks[j])] / sqrt(sum(SVD$d[seq(ranks[j])] ^ 2))
        a[[j]][, k] <- weighted_kron_sum(
          lapply(factors[[j]], "[", , idx, drop = FALSE), weights[[j]][idx]
        )
      }
    }

    # Update Y
    Y <- lapply(seq_along(a), function(j) {
      pm(A_m[[j]], a[[j]], na.rm = na.rm)
    })

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

    stopping_criteria <- c(
      drop(crossprod(unlist(a, F, F) - unlist(a_old, F, F))),
      abs(crit[iter] - crit_old)
    )

    if (crit[iter] - crit_old < -tol) stop_rgcca("Convergence issue")
    if (any(stopping_criteria < tol) | (iter > 1000)) {
      break
    }
    crit_old <- crit[iter]
    a_old <- a
    iter <- iter + 1
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
    Y = Y, a = a, factors = factors, weights = weights, crit = crit,
    AVE_inner = AVEinner, tau = tau, V = U
  )
  return(result)
}
