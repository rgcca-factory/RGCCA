tgccak <- function(A, A_m = NULL, C, tau = rep(1, length(A)),
                   scheme = "centroid", verbose = FALSE,
                   init = "svd", bias = TRUE, tol = 1e-8,
                   ranks = rep(1, length(A)),
                   n_run = 1, n_cores = 1,
                   kronecker_covariance = 0) {
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

  # Initialization of the regularization matrices
  M <- lapply(seq(J), function(j) {
    if (tau[j] == 1) {
      return(NULL)
    }
    if (j %in% B_2D) {
      return(tau[j] * diag(pjs[j]) +
        (1 - tau[j]) * crossprod(A[[j]]) / (n - 1 + bias))
    } else {
      if (kronecker_covariance == 0) {
        return(tau[j] * diag(pjs[j]) +
          (1 - tau[j]) * crossprod(A_m[[j]]) / (n - 1 + bias))
      } else if (kronecker_covariance == 1) {
        fac <- estimate_kronecker_covariance(A[[j]])
        return(tau[j] * diag(pjs[j]) +
          (1 - tau[j]) * Reduce("%x%", rev(fac)))
      } else {
        fac <- estimate_kronecker_covariance(A[[j]])
        to <- tau[j] ^ (1 / length(fac))
        fac <- lapply(fac, function(x) {
          (1 - to) * x + to * diag(nrow(x))
        })
        return(fac)
      }
    }
  })

  Minv <- lapply(seq(J), function(j) {
    if (is.null(M[[j]])) {
      return(NULL)
    }
    if (!(j %in% B_2D)) {
      if (kronecker_covariance == 2) {
        return(lapply(M[[j]], solve))
      } else {
        return(NULL)
      }
    }
    return(solve(M[[j]]))
  })

  # Matricization (mode-1)
  if (is.null(A_m)) {
    A_m <- lapply(1:J, function(x) matrix(as.vector(A[[x]]), nrow = n))
  }

  myCluster <- makeCluster(n_cores, type = "PSOCK")
  doParallel::registerDoParallel(myCluster)
  models <- foreach(run_number = seq(n_run)) %dopar% {
    init <- ifelse(run_number == 1, init, "random")

    core_tgcca(A, A_m, DIM, LEN, B_2D, B_3D, B_nD, init, g, verbose, C,
               tol, n_iter_max, bias, ranks, M, Minv)
  }
  stopCluster(myCluster)

  best_model <- which.max(lapply(models, function(m) m$crit[length(m$crit)]))
  a <- models[[best_model]]$a
  Y <- models[[best_model]]$Y
  factors <- models[[best_model]]$factors
  weights <- models[[best_model]]$weights
  crit <- models[[best_model]]$crit
  iter <- length(unlist(crit))

  # Final messages
  if (iter > 1000) {
    warning("The TGCCA algorithm did not converge after 1000 iterations.")
  }
  if (iter < 1000 & verbose) {
    cat(
      "The TGCCA algorithm converged to a stationary point after ",
      iter - 1, " iterations \n"
    )
  }
  if (verbose) {
    plot(crit[1:iter], xlab = "iteration", ylab = "criteria")
  }

  AVEinner <- sum(C * cor(Y)^2 / 2) / (sum(C) / 2)

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
