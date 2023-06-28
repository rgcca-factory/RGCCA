grgccak = function(A, C, tau = rep(1, length(A)), scheme = "centroid",
                  verbose = FALSE, init = "svd", bias = TRUE,
                  tol = 1e-08, na.rm = TRUE, ncomp = rep(1, length(A))) {

  criterion = function() {
    cur_crit = c()
    for (i in 1:J) {
      for (j in 1:J) {
        cur_crit = c(cur_crit, C[i, j] * sum(diag(g(crossprod(Y[[i]], Y[[j]])))))
      }
    }
    return(sum(cur_crit) / (n - 1 + bias))
  }

  # Returns a n x ncomp x J array
  compute_dgx = function(dg, j) {
    sapply(1:J, function(k)
      C[j, k] * Y[[k]] %*% diag(dg(diag(crossprod(Y[[j]], Y[[k]]))), nrow = ncol(Y[[k]])),
      simplify = "array"
    )
  }

  if (mode(scheme) != "function") {
    if (scheme == "horst") {g <- function(x) x ; ctrl = FALSE}
    if (scheme == "factorial") {g <- function(x)  x^2 ; ctrl = TRUE}
    if (scheme == "centroid") {g <- function(x) abs(x) ; ctrl = TRUE}
  } else {
    # check for parity of g
    g <- scheme ; ctrl = !any(g(-5:5) != g(5:-5))
  }

  dg = Deriv::Deriv(g)

  J <- length(A) # number of blocks
  n <- NROW(A[[1]]) # number of individuals
  pjs <- sapply(A, NCOL) # number of variables per block
  Z   <- array(0, dim = c(n, max(ncomp), J))

  # Initialisation by SVD
  if (init == "svd") {
    a <- lapply(1:J, function(j) svd(A[[j]], nu = 0, nv = ncomp[j])$v)
  }

  else if (init == "random") {
    a <- lapply(1:J, function(j) matrix(rnorm(pjs[j] * ncomp[j]), pjs[j]))
    a <- lapply(a, function(x) apply(x, 2, function(y) y / norm(y, type = "2")))
  }

  else {
    stop_rgcca("init should be either random or by SVD.")
  }

  Y = lapply(1:J, function(j) pm(A[[j]], a[[j]], na.rm = na.rm))

  iter <- 1
  n_iter_max <- 1000L
  crit <- numeric(n_iter_max)
  crit_old <- criterion()
  a_old = a

  repeat {
    for (j in 1:J) {
      dgx      = compute_dgx(dg, j)
      Z[, , j] = apply(dgx, c(1, 2), sum)

      Az     = pm(t(A[[j]]), Z[, 1:ncomp[j], j], na.rm = na.rm)
      SVD    = svd(Az, nv = ncomp[j], nu = ncomp[j])
      a[[j]] = SVD$u %*% t(SVD$v)
      Y[[j]] = pm(A[[j]], a[[j]], na.rm = na.rm)
    }

  crit[iter] <- criterion()
  if (verbose & (iter%%1) == 0){
    cat(" Iter: ", formatC(iter, width = 3, format = "d"),
        " Fit:", formatC(crit[iter], digits = 8,
                         width = 10, format = "f"),
        " Dif: ", formatC(crit[iter] - crit_old, digits = 8,
                          width = 10, format = "f"), "\n")
  }

  stopping_criteria = c(drop(crossprod(unlist(a, F, F) - unlist(a_old, F, F))),
                        abs(crit[iter] - crit_old))

  if (crit[iter] - crit_old < -tol) stop_rgcca("Convergence issue")
  if (any(stopping_criteria < tol) | (iter > 1000)) {break}
  crit_old = crit[iter]
  a_old <- a
  iter <- iter + 1
  }

  for (j in seq_len(J)) {
    if (ctrl & a[[j]][1] < 0) {
      a[[j]] = -a[[j]]
      Y[[j]] = pm(A[[j]] , a[[j]], na.rm = na.rm)
    }
  }

  crit <- crit[which(crit != 0)]

  if (iter > n_iter_max)
    warning("The RGCCA algorithm did not converge after", n_iter_max,
            " iterations.")
  if (iter < n_iter_max & verbose)
    cat("The RGCCA algorithm converged to a stationary point after",
            iter - 1, "iterations \n")
  if (verbose){
    plot(crit, xlab = "iteration", ylab = "criteria")
  }

  # AVEinner <- sum(C * cor(Y)^2/2)/(sum(C)/2)
  AVEinner = NULL

  result <- list(Y = Y, a = a, crit = crit,
                 AVE_inner = AVEinner, tau = tau)
  return(result)
}
