core_tgcca <- function(A, P, DIM, LEN, B_2D, B_3D, B_nD, init, g, verbose, C,
                       tol, n_iter_max, bias, ranks, M, Minv) {

  criterion <- function() {
    cur_crit <- 0
    for (i in seq(J)) {
      for (j in seq(J)) {
        cur_crit <- cur_crit + C[i, j] * g(crossprod(Y[[i]], Y[[j]]))
      }
    }
    return(drop(cur_crit) / g(n - 1 + bias))
  }

  # Compute the gradient with respect to a[[j]] : dim = pjs[j]
  compute_dgx <- function(j) {
    res <- lapply(seq(J), function(i) {
      2 * C[j, i] * Y[[i]] %*% dg(crossprod(Y[[i]], Y[[j]]))
    })
    drop(t(P[[j]]) %*% Reduce("+", res))
  }

  # Compute the update
  solution <- function(X, B, M = NULL, Md = NULL) {
    p_a <- nrow(X)
    R <- ncol(B)

    # Construct P and P^{-1/2} for change of variable
    if (is.null(M)) {
      SVD <- svd(crossprod(B))
      Pm12 <- (SVD$u %*% diag(1 / sqrt(SVD$d), nrow = R) %*% t(SVD$v)) %x%
        diag(p_a)
    } else {
      if (!is.null(Md)) {
        SVD <- svd(Md)
        Pm12 <- (SVD$u %*% diag(1 / sqrt(SVD$d), nrow = R) %*% t(SVD$v)) %x% M
      } else {
        P <- (t(B) %x% diag(p_a)) %*% M %*% (B %x% diag(p_a))
        SVD <- svd(P)
        Pm12 <- SVD$u %*% diag(1 / sqrt(SVD$d), nrow = R * p_a) %*% t(SVD$v)
      }
    }

    # Solve for u = P^{1/2}Vec(A)
    u <- Pm12 %*% as.vector(X %*% B)
    u <- u / norm(u, type = "2")

    # Invert change of variable
    A <- matrix(Pm12 %*% u, nrow = p_a, ncol = R)

    return(A)
  }

  ### Initialization
  J <- length(P)
  pjs <- sapply(DIM, function(x) prod(x[-1]))
  n   <- DIM[[1]][1]

  # Initialization of factors
  factors <- lapply(seq(J), function(j) {
    if (j %in% B_2D) {
      return(NULL)
    }
    if (init == "svd") {
      lapply(seq(2, LEN[[j]]), function(d) {
        x <- svd(matrix(aperm(
          array(P[[j]], dim = DIM[[j]]), perm = c(d, seq_along(DIM[[j]])[-d])
        ), nrow = DIM[[j]][d]), nu = ranks[j])$u
        x / norm(x, type = "F")
      })
    } else {
      lapply(seq(2, LEN[[j]]), function(d) {
        x <- svd(matrix(
          rnorm(DIM[[j]][d] * ranks[j], mean = 0, sd = 1), nrow = DIM[[j]][d]
        ), nu = ranks[j])$u
        x / norm(x, type = "F")
      })
    }
  })

  # Initialization of weights vectors a
  a <- lapply(seq(J), function(j) {
    if (j %in% B_2D) {
      if (init == "svd") {
        return(initsvd(P[[j]], dual = FALSE))
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
      if (is.list(M[[j]])) {
        a_norm <- sqrt(drop(t(rep(1, ranks[j])) %*%
          Reduce("*", lapply(seq(LEN[j] - 1), function(d) {
            t(factors[[j]][[d]]) %*% M[[j]][[d]] %*% factors[[j]][[d]]
          })) %*% rep(1, ranks[j])))
      } else {
        a_norm <- sqrt(drop(t(a[[j]]) %*% M[[j]] %*% a[[j]]))
      }
    }
    a[[j]] <- a[[j]] / a_norm
    if (!(j %in% B_2D)) {
      factors[[j]][[1]] <- factors[[j]][[1]] / a_norm
    }
  }

  # Initialization of vector Y
  Y <- lapply(seq(J), function(j) {
    P[[j]] %*% a[[j]]
  })

  # Initialize other parameters
  iter <- 1
  crit <- numeric()
  a_old <- a
  crit_old <- criterion()

  dg <- Deriv::Deriv(g)

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
        Y[[j]] <- P[[j]] %*% a[[j]]

      } else {
        for (d in seq(1, LEN[j] - 1)) {
          grad <- array(grad_j, dim = DIM[[j]][-1])
          grad <- matrix(aperm(
            grad, perm = c(d, seq_along(DIM[[j]][-1])[-d])
          ), nrow = DIM[[j]][d + 1])

          Z <- Reduce(khatri_rao, rev(factors[[j]][-d]))
          if (is.list(M[[j]])) {
            Md <- Reduce("*", lapply(seq(1, LEN[j] - 1)[-d], function(m) {
              t(factors[[j]][[m]]) %*% M[[j]][[m]] %*% factors[[j]][[m]]
            }))
            factors[[j]][[d]] <- solution(grad, Z, Minv[[j]][[d]], Md)
          } else {
            Z <- Reduce(khatri_rao, rev(factors[[j]][-d]))
            idx <- array(seq_along(grad_j), dim = DIM[[j]][-1])
            idx <- aperm(idx, perm = c(d, seq_along(DIM[[j]][-1])[-d]))
            factors[[j]][[d]] <- solution(grad, Z, M[[j]][idx, idx])
          }

          if (d > 1) {
            f_norm <- norm(factors[[j]][[d]], type = "F")
            factors[[j]][[d]] <- factors[[j]][[d]] / f_norm
            factors[[j]][[1]] <- factors[[j]][[1]] * f_norm
          }

          a[[j]] <- Reduce(khatri_rao, rev(factors[[j]])) %*% rep(1, ranks[j])
          Y[[j]] <- P[[j]] %*% a[[j]]
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

  Y <- do.call(cbind, Y)
  weights <- lapply(seq(J), function(j) rep(1, ranks[j]))

  return(list(a = a, factors = factors, weights = weights, Y = Y, crit = crit))
}
