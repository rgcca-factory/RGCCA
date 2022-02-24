gsgccak_penalized = function(A, C, sparsity = rep(1, length(A)), scheme = "centroid",
                             verbose = FALSE, init = "svd", bias = TRUE,
                             tol = 1e-08, na.rm = TRUE,
                             ncomp = rep(1, length(A)), penalty_coef = 0) {

  ### Utility functions
  penalty = function() {
    cur_penalty = c()
    for (j in 1:J) {
      tmp = crossprod(a[[j]])^2
      # tmp = crossprod(Y[[j]])^2
      diag(tmp) = 0
      cur_penalty = c(cur_penalty, sum(tmp))
    }
    return(sum(cur_penalty))
  }

  criter = function() {
    cur_crit = c()
    for (i in 1:J) {
      for (j in 1:J) {
        cur_crit = c(cur_crit, C[i, j] * sum(diag(g(crossprod(Y[[i]], Y[[j]])))))
      }
    }
    return(sum(cur_crit))
  }

  criterion = function() {
    return(criter() - penalty_coef * penalty())
  }

  alpha_i_k = function(i, k) {
    svd(crossprod(a[[i]][, -k]))$d[1]
    # svd(
    #   t(a[[i]][, -k]) %*% (t(A[[i]]) %*% (A[[i]] %*% (t(A[[i]]) %*% Y[[i]][, -k])))
    # )$d[1]
  }

  compute_dgx = function(dg, j, k) {
    sapply(1:J, function(i)
      C[j, i] * Y[[i]][, k] * dg(drop(crossprod(Y[[j]][, k], Y[[i]][, k])))
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

  dg = Deriv::Deriv(g, env = parent.frame())

  J <- length(A) # number of blocks
  n <- NROW(A[[1]]) # number of individuals
  pjs <- sapply(A, NCOL) # number of variables per block
  Z <- matrix(0,NROW(A[[1]]),J)

  if (any( sparsity < 1/sqrt(pjs) | sparsity > 1 ))
    stop_rgcca("L1 constraints must vary between 1/sqrt(p_j) and 1.")

  const <- sparsity * sqrt(pjs)

  A <- lapply(A, as.matrix)
  a <- list()

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

  for (j in 1:J) {
    a[[j]] <- apply(a[[j]], 2, function(x) soft.threshold(x, const[j]))
  }
  Y = lapply(1:J, function(j) pm(A[[j]], a[[j]], na.rm = na.rm))

  iter <- 1
  # n_iter_max <- 1000L
  n_iter_max <- 100000L
  crit <- numeric(n_iter_max)
  crit_first_part = numeric(n_iter_max)
  crit_second_part = numeric(n_iter_max)
  crit_old <- criterion()
  a_old = a

  repeat {
    for (j in 1:J) {
      for (k in 1:ncomp[j]) {
        dgx      = compute_dgx(dg, j, k)
        Z[, j]   = apply(dgx, 1, sum)

        Q       = t(A[[j]]) %*% Z[, j]
        alpha   = alpha_i_k(j, k)

        a[[j]][, k]  = Q - 2 * penalty_coef * (a[[j]][, -k] %*% crossprod(a[[j]][, -k], a[[j]][, k]) - alpha * a[[j]][, k])
        a[[j]][, k] = soft.threshold(a[[j]][, k], const[j])
        Y[[j]][, k] = A[[j]] %*% a[[j]][, k]
      }
    }

    crit[iter] <- criterion()
    crit_first_part[iter] <- criter()
    crit_second_part[iter] <- penalty_coef * penalty()

    if (verbose & (iter%%1) == 0){
      cat(" Iter: ", formatC(iter, width = 3, format = "d"),
          " Fit:", formatC(crit[iter], digits = 8,
                           width = 10, format = "f"),
          " Dif: ", formatC(crit[iter] - crit_old, digits = 8,
                            width = 10, format = "f"), "\n")
    }

    stopping_criteria = c(drop(crossprod(unlist(a, F, F) - unlist(a_old, F, F))),
                          abs(crit[iter] - crit_old))
    if (any(stopping_criteria < tol) | (iter > n_iter_max)) {break}
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
    plot(crit_first_part[1:iter], xlab = "iteration", ylab = "Criterion (correlations)")
    plot(crit_second_part[1:iter], xlab = "iteration", ylab = "Criterion (penalty)")
  }

  # AVEinner <- sum(C * cor(Y)^2/2)/(sum(C)/2)
  AVEinner = NULL

  result <- list(Y = Y, a = a, crit = crit,
                 AVE_inner = AVEinner)
  return(result)
}
