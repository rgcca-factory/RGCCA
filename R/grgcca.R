grgcca = function(blocks, connection = 1 - diag(length(blocks)),
                  tau = rep(1, length(blocks)),
                  ncomp = rep(1, length(blocks)), scheme = "centroid",
                  init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,
                  na.rm = TRUE, quiet = FALSE)
{

  if (mode(scheme) != "function") {
    if (verbose)
      cat("Computation of the gRGCCA block components based on the",
          scheme, "scheme \n")
  }
  if (mode(scheme) == "function" & verbose)
    cat("Computation of the gRGCCA block components based on the g scheme \n")

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
  nb_ind <- NROW(blocks[[1]])
  AVE_X = list()
  AVE_outer <- rep(NA,max(ncomp))

  Y <- NULL
  a <- astar <- NULL
  crit <- list()
  AVE_inner <- rep(NA,max(ncomp))

  # Save computed shrinkage parameter in a new variable
  tau = rep(1, length(blocks))

  # First component block
  gcca.result <- grgccak(blocks, connection, tau = tau, scheme = scheme,
                         init = init, bias = bias, tol = tol,
                         verbose = verbose, na.rm = na.rm, ncomp = ncomp)

  Y         <- gcca.result$Y
  a         <- gcca.result$a
  astar     <- gcca.result$a
  AVE_inner <- gcca.result$AVE_inner
  crit      <- gcca.result$crit

  for (b in 1:J) {
    rownames(a[[b]]) = rownames(astar[[b]]) = colnames(blocks[[b]])
    rownames(Y[[b]]) = rownames(blocks[[b]])
    colnames(Y[[b]]) = paste0("comp", 1:max(ncomp))
  }

  for (j in 1:J) AVE_X[[j]] = apply(
    cor(blocks[[j]], Y[[j]], use = "pairwise.complete.obs")^2, 2, mean)

  outer = matrix(unlist(AVE_X), nrow = max(ncomp))

  for (j in 1:max(ncomp))
    AVE_outer[j] <- sum(pjs * outer[j,])/sum(pjs)

  Y = shave(Y, ncomp)
  AVE_X = shave(AVE_X, ncomp)

  AVE <- list(AVE_X = AVE_X, AVE_outer = AVE_outer, AVE_inner = AVE_inner)

  out <- list(Y = shave(Y, ncomp), a = shave(a,ncomp),
              astar = shave(astar, ncomp),
              tau = tau,
              crit = crit,
              AVE = AVE)

   class(out) <- "grgcca"

  return(out)
}
