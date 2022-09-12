core_rgccad <- function(A, C, tau, scheme, init, bias, tol,
                        verbose, na.rm, n_run, ncomp, scale, scale_block) {
  ndefl <- ncomp - 1
  N <- max(ndefl)
  nb_ind <- NROW(A[[1]])
  J <- length(A)
  pjs <- sapply(A, NCOL)

  # One component per block
  if(N == 0){
    result <- rgccak(A, C, tau = tau, scheme = scheme, init = init,
                     bias = bias, tol = tol, verbose = verbose,
                     na.rm=na.rm,
                     scale_block=scale_block, scale=scale)

    Y <- NULL
    for (b in 1:J) Y[[b]] <- result$Y[, b, drop = FALSE]
    a <- lapply(result$a, cbind)

    for (b in 1:J) {
      rownames(a[[b]]) = colnames(A[[b]])
      rownames(Y[[b]]) = rownames(A[[b]])
      colnames(Y[[b]]) = "comp1"
    }

    tau=result$tau

    if(is.vector(tau))
      names(tau) = names(A)

    return(list(a = a, Y = Y, astar = a, AVE_inner = result$AVE_inner, crit = result$crit, tau = tau))
  }

  Y <- NULL
  crit = list()
  AVE_inner <- rep(NA, max(ncomp))
  R <- A
  P <- a <- astar <- NULL
  if (is.numeric(tau))
  {
    tau_mat = tau
    if(is.vector(tau_mat)) names(tau_mat) = names(A)
    if(is.matrix(tau_mat))colnames(tau_mat) = names(A)
  }
  else
  {
    tau_mat = matrix(NA, max(ncomp), J)
    colnames(tau_mat) = names(A)
  }

  for (b in 1:J) P[[b]] <- a[[b]] <- astar[[b]] <- matrix(NA, pjs[[b]], N + 1)
  for (b in 1:J) Y[[b]] <- matrix(NA, nb_ind, N + 1)
  for (n in 1:N) {
    if (verbose)
      cat(paste0("Computation of the RGCCA block components #", n, " is under
                 progress...\n"))
    if (is.vector(tau))
      rgcca.result <- rgccak(R, C, tau = tau, scheme = scheme,init = init,
                             bias = bias, tol = tol, verbose = verbose,
                             na.rm = na.rm)
    else rgcca.result <- rgccak(R, C, tau = tau[n, ], scheme = scheme,
                                init = init, bias = bias, tol = tol,
                                verbose = verbose, na.rm = na.rm)

    if (!is.numeric(tau)) tau_mat[n, ] = rgcca.result$tau

    AVE_inner[n] <- rgcca.result$AVE_inner
    crit[[n]] <- rgcca.result$crit

    # deflation
    for (b in 1:J) Y[[b]][, n] <- rgcca.result$Y[, b]
    defla.result <- defl.select(rgcca.result$Y, R, ndefl , n, nbloc = J)
    R <- defla.result$resdefl
    for (b in 1:J) P[[b]][, n] <- defla.result$pdefl[[b]]
    for (b in 1:J) a[[b]][, n] <- rgcca.result$a[[b]]
    if (n == 1)
    {
      for (b in 1:J) astar[[b]][, n] <- rgcca.result$a[[b]]
    }
    else {
      for (b in 1:J)
      {
        astar[[b]][, n] <- rgcca.result$a[[b]]-astar[[b]][, (1:n-1), drop = F]%*%
          drop(t(a[[b]][, n])%*%P[[b]][, 1:(n-1), drop = F])
      }
    }
  }
  if (verbose)
    cat(paste0("Computation of the RGCCA block components #",
               N + 1, " is under progress ... \n"))
  if (is.vector(tau))
    rgcca.result <- rgccak(R, C, tau = tau, scheme = scheme, init = init,
                           bias = bias, tol = tol, verbose = verbose)
  else rgcca.result <- rgccak(R, C, tau = tau[N+1, ], scheme = scheme,
                              init = init, bias = bias, tol = tol,
                              verbose = verbose)

  crit[[N+1]] <- rgcca.result$crit
  if (!is.numeric(tau))
    tau_mat[N+1, ] = rgcca.result$tau
  AVE_inner[max(ncomp)] <- rgcca.result$AVE_inner
  for (b in 1:J) {
    Y[[b]][, N+1] <- rgcca.result$Y[, b]
    a[[b]][, N+1] <- rgcca.result$a[[b]]
    astar[[b]][, N+1] <- rgcca.result$a[[b]]-astar[[b]][, (1:N), drop = F]%*%
      drop(t(a[[b]][, (N+1)]) %*%P[[b]][, 1:N, drop = F])
    rownames(a[[b]]) = rownames(astar[[b]]) = colnames(A[[b]])
    rownames(Y[[b]]) = rownames(A[[b]])
    colnames(Y[[b]]) = paste0("comp", 1:max(ncomp))
  }

  return(list(a = a, Y = Y, astar = astar, AVE_inner = AVE_inner, crit = crit, tau = tau_mat))
}
