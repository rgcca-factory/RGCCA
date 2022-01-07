gmgccak <- function (A, A_m = NULL, C, tau = rep(1, length(A)), scheme = "centroid",
                    verbose = FALSE, init="svd", bias = TRUE, tol = 1e-8,
                    regularisation_matrices, ranks= rep(1, length(A)),
                    ncomp = rep(1, length(A)), penalty_coef = 0) {

  call = list(A = A, A_m = A_m, C = C, scheme = scheme, verbose = verbose, init = init,
              bias = bias, tol = tol, ranks = ranks, ncomp = ncomp)

  penalty = function() {
    cur_penalty = c()
    for (j in 1:J) {
      cur_penalty = c(cur_penalty, sum(crossprod(a[[j]])^2))
    }
    return(sum(cur_penalty))
  }

  criterion = function() {
    cur_crit = c()
    for (i in 1:J) {
      for (j in 1:J) {
        cur_crit = c(cur_crit, C[i, j] * sum(diag(g(crossprod(Y[[i]], Y[[j]])))))
      }
    }
    return(sum(cur_crit) - penalty_coef * penalty())
  }

  M_i_k = function(i, k) {
    tcrossprod(a[[i]][, -k])
  }

  factors_comp = function(i, k) {
    lapply(factors[[i]], "[", seq((k - 1) * ranks[i] + 1, k * ranks[i]))
  }

  # Returns a N x J matrix
  compute_dgx = function(dg, j, k) {
    sapply(1:J, function(i)
      C[j, i] * Y[[i]][, k] * dg(drop(crossprod(Y[[j]][, k], Y[[i]][, k])))
    )
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

  # Dimensions of each block
  pjs <- sapply(DIM, function(x) prod(x[-1]))
  n   <- DIM[[1]][1]

  Y   <- matrix(0, n, J)

  # Matricization (mode-1)
  if(is.null(A_m)){
    A_m = lapply(1:J, function(x) matrix(as.vector(A[[x]]), nrow = n))
  }

  a <- factors <- weights <- M_inv_sqrt <- P <- P_m <- list()
  for (j in 1:J) {
    factors[[j]] <- list()
  }

  # Initialization of vector a (weight vector)
  for (j in 1:J) {
    if (is.list(init)) {
      if (j %in% B_nD) {
        factors[[j]] <- init$factors[[j]]
        weights[[j]] <- init$weights[[j]]
        a[[j]]       <- list_khatri_rao(factors[[j]])
      } else {
        a[[j]] <- init$a[[j]][[1]]
      }
    } else if (init=="svd") {
      # SVD Initialization of a_j
      if (j %in% B_2D) {
        SVD    = svd(A[[j]], nu=0, nv=ncomp[j])
        a[[j]] = SVD$v

      } else {
        for (d in 1:(LEN[[j]] - 1)) {
          SVD = svd(apply(A[[j]], d+1, c), nu=0, nv=ncomp[j] * ranks[j])
          factors[[j]][[d]] = SVD$v
        }
        weights[[j]] = matrix(1 / sqrt(ranks[j]), nrow = ncomp[j], ncol = ranks[j])
        a[[j]]       = mix_weighted_kron_sum(factors[[j]], weights[[j]], ncomp[j], ranks[j])
      }
    } else if (init == "random") {
      # Random Initialisation of a_j
      A_random <- array(rnorm(n = pjs[[j]], mean = 0, sd = 1), dim = DIM[[j]][-1])
      if (j %in% B_2D) {
        SVD    = svd(A_random, nu=0, nv=ncomp[j])
        a[[j]] = SVD$v
      } else {
        for (d in 1:(LEN[[j]] - 1)) {
          SVD = svd(apply(A_random, d, c), nu=0, nv=ncomp[j] * ranks[j])
          factors[[j]][[d]] = SVD$v
        }
        weights[[j]] = matrix(1 / sqrt(ranks[j]), nrow = ncomp[j], ncol = ranks[j])
        a[[j]]       = mix_weighted_kron_sum(factors[[j]], weights[[j]], ncomp[j], ranks[j])
      }
    } else {
      stop_rgcca("init should be either random or by SVD.")
    }
  }
  # Initialization of vector Y
  Y = lapply(1:J, function(j) A_m[[j]] %*% a[[j]])

  # Determination of the regularization matrix
  for (j in 1:J){
    if (j > length(regularisation_matrices)) {
      reg_matrices = NULL
    } else {
      reg_matrices = regularisation_matrices[[j]]
    }
    reg_matrices = parse_regularisation_matrices(
      reg_matrices = reg_matrices,
      tau          = tau[j],
      A            = A[[j]],
      DIM          = DIM[[j]],
      j            = j,
      bias         = bias
    )
    P[[j]]          = reg_matrices$P_t
    P_m[[j]]        = reg_matrices$P
    M_inv_sqrt[[j]] = reg_matrices$M_inv_sqrt
    tau[j]          = reg_matrices$tau
  }

  # Initialize other parameters
  crit_old = criterion()
  iter     = 1
  crit     = numeric()
  Z        = matrix(0, n, J)
  a_old    = a

  dg = Deriv::Deriv(g, env = parent.frame())

  # MGCCA algorithm
  repeat {
    for (k in 1:max(ncomp)) {
      for (j in 1:J){
        if (k > ncomp[j]) next
        # Apply the derivative on the current variables
        dgx      = compute_dgx(dg, j, k)
        Z[, j]   = apply(dgx, 1, sum)

        if (j %in% B_2D) { # Matrices
          Q       = t(P[[j]]) %*% Z[, j]
          M       = penalty_coef * M_i_k(j, k)
          alpha   = svd(M, nu = 0, nv = 0)$d[1]

          a[[j]][, k]  = Q - 2 * crossprod(M - alpha * diag(nrow(M)), a[[j]][, k])
          a[[j]][, k]  = a[[j]][, k] / norm(a[[j]][, k], type = "2")
          Y[[j]][, k]  = P[[j]] %*% a[[j]][, k]
        } else {
          for (d in 1:(LEN[j] - 1)) {
            fac   = factors_comp(j, k)
            W     = list_khatri_rao(fac[-d]) %*% diag(weights[[j]][k, ])

            Q     = array(t(P[[j]]) %*% Z[, j], dim = DIM[[j]][-1])
            Q     = unfold(Q, mode = d)
            Q     = crossprod(W, Q)

            M     = t(W) %*% (penalty_coef * M_i_k(j, k)) %*% W
            alpha = svd(M, nu = 0, nv = 0)$d[1]

            SVD = svd(Q - 2 * (M + alpha * diag(nrow(M))) %*% fac[[d]])

            fac[[d]]    = SVD$u %*% t(SVD$v)
            factors[[j]][[d]][seq((k - 1) * ranks[j] + 1, k * ranks[j])] = fac[[d]]
            a[[j]][, k] = weighted_kron_sum(fac, weights[[j]][k, ])
            Y[, j][, k] = P[[j]] %*% a[[j]][, k]

            dgx      = compute_dgx(dg, j, k)
            Z[, j]   = apply(dgx, 1, sum)
          }
          fac   = factors_comp(j, k)
          W     = list_khatri_rao(fac)

          Q       = t(W) %*% t(P[[j]]) %*% Z[, j]
          M       = t(W) %*% (penalty_coef * M_i_k(j, k)) %*% W
          alpha   = svd(M, nu = 0, nv = 0)$d[1]

          weights[[j]][k, ]  = Q - 2 * crossprod(M - alpha * diag(nrow(M)), weights[[j]][k, ])
          weights[[j]][k, ]  = weights[[j]][k, ] / norm(weights[[j]][k, ], type = "2")
          a[[j]][, k]        = weighted_kron_sum(fac, weights[[j]][k, ])
          Y[[j]][, k]        = P[[j]] %*% a[[j]][, k]
        }
      }
    }

    # Store previous criterion
    crit[iter] <- criterion()

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

  # Inverse change of variables if needed
  if (length(M_inv_sqrt) > 0)  { # If no regularization matrix, list is empty
    stop_rgcca("NOT IMPLEMENTED")
    for (j in 1:J) {
      if (j <= length(M_inv_sqrt) && !is.null(M_inv_sqrt[[j]])) {
        if (j %in% B_2D) {
          a[[j]] = M_inv_sqrt[[j]] %*% a[[j]]
        } else {
          for (d in 1:(LEN[[j]] - 1)) {
            factors[[j]][[d]] = M_inv_sqrt[[j]][[d]] %*% factors[[j]][[d]]
          }
          a[[j]] = list_khatri_rao(factors[[j]])
        }
      }
    }
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
