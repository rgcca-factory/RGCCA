nn_gmgccak_penalized <- function(A, A_m = NULL, C, tau = rep(1, length(A)), scheme = "centroid",
                              verbose = FALSE, init="svd", bias = TRUE, tol = 1e-8,
                              ncomp = rep(1, length(A)), penalty_coef = 1, na.rm = TRUE,
                              orth_Y = FALSE) {

  call = list()

  if (mode(scheme) != "function") {
    if (scheme == "horst") {g <- function(x) x ; ctrl = FALSE}
    if (scheme == "factorial") {g <- function(x)  x^2 ; ctrl = TRUE}
    if (scheme == "centroid") {g <- function(x) abs(x) ; ctrl = TRUE}
  } else {
    # check for parity of g
    g <- scheme ; ctrl = !any(g(-5:5) != g(5:-5))
  }

  A_sing <- lapply(A_m, function(x) svd(x)$d[1])

  # List of 2D matrix and higher order tensors
  DIM    <- lapply(A, dim)
  LEN    <- unlist(lapply(DIM, length))
  B_2D   <- which(LEN == 2)   # Store which blocks are 2D
  B_0D   <- which(LEN == 0)   # Store which blocks are 1D (stored as 0D)

  # Convert vectors to one-column matrices
  if (length(B_0D) != 0) {
    for (i in B_0D) {
      A[[i]]   = as.matrix(A[[i]])
      DIM[[i]] = dim(A[[i]])
    }
    B_2D = c(B_2D, B_0D)
  }

  n_iter_max <- 1000L

  models <- pbapply::pblapply(seq(100), function(run_number) {
    init <- ifelse(run_number == 1, "svd", "random")

    core_nn_gmgcca(A, A_m, ncomp, init, g, penalty_coef, A_sing, verbose, C,
                   orth_Y, tol, n_iter_max)

  }, cl = 6)

  best_model <- which.max(lapply(models, function(m) m$crit[length(m$crit)]))
  a <- models[[best_model]]$a
  Y <- models[[best_model]]$Y
  factors <- models[[best_model]]$factors
  crit <- models[[best_model]]$crit
  iter <- length(crit)

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
