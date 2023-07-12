gmgcca_pdd = function(blocks, connection = 1 - diag(length(blocks)),
                              tau = rep(1, length(blocks)),
                              ncomp = rep(1, length(blocks)), scheme = "centroid",
                              init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,
                              na.rm = TRUE, quiet = FALSE, n_run = 1, n_cores = 1,
                              scale = TRUE, scale_block = TRUE, prescaling = FALSE,
                              eta_decay = 0.9, rho_decay = 0.9, rho = 2,
                              penalty_coef = 1)
{

  if (mode(scheme) != "function") {
    if (verbose)
      cat("Computation of the gMGCCA_PDD block components based on the",
          scheme, "scheme \n")
  }
  if (mode(scheme) == "function" & verbose)
    cat("Computation of the gMGCCA_PDD block components based on the g scheme \n")

  if (!is.numeric(tau) & verbose) {
    cat("Optimal Shrinkage intensity paramaters are estimated \n")
  }
  else {
    if (is.numeric(tau) & verbose) {
      cat("Shrinkage intensity paramaters are chosen manually \n")
    }
  }

  if (!prescaling) {
    blocks <- scaling(
      blocks, scale = scale, bias = bias, scale_block = scale_block
    )
  }

  J <- length(blocks)
  pjs <- sapply(blocks,NCOL)
  n <- NROW(blocks[[1]])
  AVE_X = list()
  AVE_outer <- rep(NA,max(ncomp))

  block_names <- names(blocks)

  # Change of variables if tau != 1
  M <- lapply(seq(J), function(j) {
    if (tau[j] == 1) return(NULL)
    x <- (1 - tau[j]) * crossprod(blocks[[j]]) / (n - 1 + bias) + tau[j] * diag(pjs[j])
    eig <- eigen(x, symmetric = TRUE)
    return(eig$vectors %*% diag(1 / sqrt(eig$values)) %*% t(eig$vectors))
  })

  blocks <- lapply(seq(J), function(j) {
    if (tau[j] == 1) return(blocks[[j]])
    blocks[[j]] %*% M[[j]]
  })
  names(blocks) <- block_names

  # Matricization (mode-1)
  A_m = lapply(1:J, function(x) {
    m = matrix(as.vector(blocks[[x]]), nrow = n)
    rownames(m) = rownames(blocks[[x]])
    if (!is.null(dimnames(blocks[[x]]))) {
      grid        = do.call(expand.grid, dimnames(blocks[[x]])[-1])
      colnames(m) = do.call(paste, c(grid, sep = " x "))
    }
    return(m)
  })

  # Call gMGCCA_PDD_nuclear
  # models <- pbapply::pblapply(seq(n_run), function(run_number) {
  #   init <- ifelse(run_number == 1, init, "random")
  #
  #   gmgccak_PDD_nuclear(blocks, connection, tau = tau, scheme = scheme,
  #                       init = init, bias = bias, tol = tol,
  #                       verbose = verbose, na.rm = na.rm, ncomp = ncomp,
  #                       eta_decay = eta_decay, rho_decay = rho_decay, rho = rho,
  #                       penalty_coef = penalty_coef, A_m = A_m)
  #
  # }, cl = n_cores)

  models <- lapply(seq(n_run), function(run_number) {
    init <- ifelse(run_number == 1, init, "random")

    gmgccak_PDD5(blocks, connection, tau = tau, scheme = scheme,
                         init = init, bias = bias, tol = tol,
                         verbose = verbose, na.rm = na.rm, ncomp = ncomp,
                         eta_decay = eta_decay, rho_decay = rho_decay, rho = rho,
                         penalty_coef = penalty_coef, A_m = A_m)

  })

  best_model <- which.max(
    lapply(models, function(m) m$crit[length(m$crit)])
  )
  a <- models[[best_model]]$a
  astar <- models[[best_model]]$a
  Y <- models[[best_model]]$Y
  crit <- models[[best_model]]$crit
  AVE_inner <- models[[best_model]]$AVE_inner

  # Invert change of variables
  a <- lapply(seq(J), function(j) {
    if (tau[j] == 1) return(a[[j]])
    M[[j]] %*% a[[j]]
  })

  for (b in 1:J) {
    rownames(a[[b]]) = rownames(astar[[b]]) = colnames(A_m[[b]])
    rownames(Y[[b]]) = rownames(blocks[[b]])
    colnames(Y[[b]]) = paste0("comp", 1:max(ncomp))
  }

  for (j in 1:J) AVE_X[[j]] = apply(
    cor(A_m[[j]], Y[[j]], use = "pairwise.complete.obs")^2, 2, mean)

  outer = matrix(unlist(AVE_X), nrow = max(ncomp))

  for (j in 1:max(ncomp))
    AVE_outer[j] <- sum(pjs * outer[j,])/sum(pjs)

  Y = shave.matlist(Y, ncomp)
  AVE_X = shave.veclist(AVE_X, ncomp)

  AVE <- list(
    AVE_X = AVE_X,
    AVE_outer_model = AVE_outer,
    AVE_inner_model = AVE_inner
  )

  out <- list(
    Y = Y,
    a = shave.matlist(a, ncomp),
    astar = shave.matlist(astar, ncomp),
    tau = tau,
    crit = crit,
    AVE = AVE
  )

  class(out) <- "gmgcca_pdd_nuclear"

  return(out)
}
