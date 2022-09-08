gmgcca_penalized = function(blocks, connection = 1 - diag(length(blocks)),
                            tau = rep(1, length(blocks)),
                            ncomp = rep(1, length(blocks)), scheme = "centroid",
                            init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,
                            quiet = FALSE, penalty_coef = 0)
{

  if (mode(scheme) != "function") {
    if (verbose)
      cat("Computation of the gRGCCA_penalized block components based on the",
          scheme, "scheme \n")
  }
  if (mode(scheme) == "function" & verbose)
    cat("Computation of the gRGCCA_penalized block components based on the g scheme \n")

  if (!is.numeric(tau) & verbose) {
    cat("Optimal Shrinkage intensity paramaters are estimated \n")
  }
  else {
    if (is.numeric(tau) & verbose) {
      cat("Shrinkage intensity paramaters are chosen manually \n")
    }
  }

  J <- length(blocks)
  pjs <- sapply(blocks,NCOL)
  AVE_X = list()
  AVE_outer <- rep(NA,max(ncomp))

  n <- NROW(blocks[[1]])
  # Matricization (mode-1)
  A_m = lapply(seq(J), function(x) {
    m = matrix(as.vector(blocks[[x]]), nrow = n)
    rownames(m) = rownames(blocks[[x]])
    if (!is.null(dimnames(blocks[[x]]))) {
      grid        = do.call(expand.grid, dimnames(blocks[[x]])[-1])
      colnames(m) = do.call(paste, c(grid, sep = " x "))
    }
    return(m)
  })

  Y <- NULL
  a <- astar <- NULL
  crit <- list()
  AVE_inner <- rep(NA,max(ncomp))

  # Save computed shrinkage parameter in a new variable
  tau = rep(1, length(blocks))

  # First component block
  gcca.result <- nn_gmgccak_penalized(blocks, A_m, connection, tau = tau, scheme = scheme,
                                   init = init, bias = bias, tol = tol,
                                   verbose = verbose, na.rm = na.rm, ncomp = ncomp,
                                   penalty_coef = penalty_coef)

  Y         <- gcca.result$Y
  a         <- gcca.result$a
  astar     <- gcca.result$a
  AVE_inner <- gcca.result$AVE_inner
  crit      <- gcca.result$crit
  factors <- gcca.result$factors

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

  AVE <- list(AVE_X = AVE_X, AVE_outer = AVE_outer, AVE_inner = AVE_inner)

  out <- list(Y = Y, a = a,
              astar = astar,
              tau = tau,
              crit = crit,
              AVE = AVE,
              factors = factors)

  class(out) <- "gmgcca_penalized"

  return(out)
}
