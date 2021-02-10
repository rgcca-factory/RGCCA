mgccak <- function (A, A_m = NULL, C, tau = rep(1, length(A)), scheme = "centroid",
                    verbose = FALSE, init="svd", bias = TRUE, tol = 1e-8,
                    regularisation_matrices, ranks= rep(1, length(A))) {
  
  kron_sum <- function(factors) {
    apply(Reduce("krprod", rev(factors)), 1, sum)
  }

  call = list(A = A, A_m = A_m, C = C, scheme = scheme, verbose = verbose, init = init,
              bias = bias, tol = tol, ranks = ranks)

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

  a <- factors <- reg_matrices <- M_inv <- P <- list()
  for (j in 1:J) {
    factors[[j]] <- list()
  }

  # Initialization of vector a (weight vector)
  for (j in 1:J) {
    if (is.list(init)) {
      if (j %in% B_nD) {
        factors[[j]] <- init$factors[[j]]
        a[[j]]       <- kron_sum(factors[[j]])
      } else {
        a[[j]] <- init$a[[j]][[1]]
      }
    } else if (init=="svd") {
      # SVD Initialization of a_j
      if (j %in% B_2D) {
        a[[j]] <- initsvd(A[[j]], dual)
      } else {
        SVD               <- svd(apply(A[[j]], 2, c), nu=0, nv=ranks[[j]])
        factors[[j]][[1]] <- SVD$v %*%
          diag(SVD$d[1:ranks[[j]]]) / sqrt(sum(SVD$d[1:ranks[[j]]] ^ 2))
        for (d in 2:(LEN[[j]] - 1)) {
          factors[[j]][[d]] <- svd(apply(A[[j]], d+1, c), nu=0, nv=ranks[[j]])$v
        }
        a[[j]] <- kron_sum(factors[[j]])
      }
    } else if (init=="random") {
      # Random Initialisation of a_j
      A_random <- array(mvrnorm(n = pjs[j], mu = 0, Sigma = 1), dim = DIM[[j]])
      if (j %in% B_2D) {
        a[[j]] <- initsvd(A_random)
      } else {
        SVD               <- svd(apply(A_random, 2, c), nu=0, nv=ranks[[j]])
        factors[[j]][[1]] <- SVD$v %*%
          diag(SVD$d[1:ranks[[j]]]) / sqrt(sum(SVD$d[1:ranks[[j]]] ^ 2))
        for (d in 2:(LEN[[j]] - 1)) {
          factors[[j]][[d]] <- svd(apply(A[[j]], d+1, c), nu=0, nv=ranks[[j]])$v
        }
        a[[j]] <- kron_sum(factors[[j]])
      }
    } else {
      stop_rgcca("init should be either random or by SVD.")
    }
  }
  # Initialization of vector Y
  for (j in 1:J) Y[, j] <- (n^(-1/2)) * A_m[[j]] %*% a[[j]]

  # Determination of the M regularization matrix
  for (j in 1:J){
    reg_matrices[[j]] = parse_regularisation_matrices(
      reg_matrices = regularisation_matrices[[j]],
      tau          = tau[j],
      A            = A[[j]],
      DIM          = DIM[[j]],
      bias         = bias
    )
    P[[j]]     = reg_matrices[[j]]$P
    M_inv[[j]] = reg_matrices[[j]]$M_inv
    tau[j]     = reg_matrices[[j]]$tau
  }

  # Initialize other parameters
  crit_old = sum(C * g(cov2(Y, bias = bias)))
  iter     = 1
  crit     = numeric()
  Z        = matrix(0, n, J)
  a_old    = a

  dg = Deriv::Deriv(g, env = parent.frame())

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
        factors[[j]][[1]] = SVD$v %*%
          diag(SVD$d[1:ranks[[j]]]) / sqrt(sum(SVD$d[1:ranks[[j]]] ^ 2))
        factors[[j]][[2]] = SVD$u
        a[[j]]            = kron_sum(factors[[j]])
        Y[, j]            = P[[j]] %*% a[[j]]

      } else if (j %in% B_nD) { # higher order Tensors
        col_idx           = 1:(LEN[[j]] - 1)
        Q                 = array(t(P[[j]]) %*% Z[, j], dim = DIM[[j]][-1])
        Q                 = unfold(Q, mode = 1)
        other_factors     = kron_sum(factors[[j]][-1])
        D                 = diag(sqrt(diag(crossprod(factors[[j]][[1]]))))
        # Tandem iteration
        SVD               = svd(x = Q %*% other_factors %*% D, nu = ranks[[j]],
                                nv = ranks[[j]])
        factors[[j]][[1]] = SVD$u %*% t(SVD$v)
        weights           = diag(t(Q %*% other_factors) %*% factors[[j]][[1]])
        D                 = diag(weights / sqrt(sum(weights ^ 2)))
        factors[[j]][[1]] = factors[[j]][[1]] %*% D
        a[[j]]            = kron_sum(factors[[j]])
        Y[, j]            = P[[j]] %*% a[[j]]

        for (d in 2:(LEN[[j]] - 1)) {
          dgx              = update_dgx(scheme, Y, dg, n, J, j)
          Z[, j]           = rowSums(
            matrix(rep(C[j, ], n), n, J, byrow = TRUE) * dgx * Y)

          Q                 = array(t(P[[j]]) %*% Z[, j], dim = DIM[[j]][-1])
          Q                 = unfold(Q, mode = d)
          other_factors     = kron_sum(factors[[j]][-d])
          SVD               = svd(x = Q %*% other_factors, nu = ranks[[j]],
                                  nv = ranks[[j]])
          factors[[j]][[d]] = SVD$u %*% t(SVD$v)
          a[[j]]            = kron_sum(factors[[j]])
          Y[, j]            = P[[j]] %*% a[[j]]
        }

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
      drop(crossprod(Reduce("c", mapply("-", a, a_old)))),
      crit[iter] - crit_old
    )
    # Criterion must increase
    if ( crit[iter] - crit_old < -tol)
      {stop_rgcca("Convergence error: criterion did not increase monotonously")}
    if (any(stopping_criteria < tol) | (iter > 1000))
      {break}
    crit_old = crit[iter]
    a_old <- a
    iter <- iter + 1
  }

  # Inverse change of variables if needed
  if (length(M_inv) > 0)  { # If no regularization matrix, list is empty
    for (j in 1:J) {
      if (!is.null(M_inv[[j]])) {
        if (j %in% B_2D) {
          a[[j]] = M_inv[[j]]$Minv_sqrt %*% a[[j]]
        } else {
          for (d in 1:(LEN[[j]] - 1)) {
            factors[[j]][[d]] = M_inv[[j]][[d]]$Minv_sqrt %*% factors[[j]][[d]]
          }
          a[[j]] = kron_sum(factors[[j]])
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

  AVEinner <- sum(C * cor(Y)^2/2)/(sum(C)/2)

  result <- list(Y         = Y,
                 a         = a,
                 factors   = factors,
                 crit      = crit,
                 AVE_inner = AVEinner,
                 call      = call,
                 tau       = tau)

  return(result)
}
