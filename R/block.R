### Create classes
new_block <- function(x, j, na.rm = TRUE, bias = TRUE,
                      ..., class = character()) {
  n <- NROW(x)
  p <- NCOL(x)
  N <- ifelse(bias, n, n - 1)

  x <- list(
    x = x,
    j = j,
    n = n,
    p = p,
    N = N,
    na.rm = na.rm,
    a = NULL,
    Y = NULL,
    ...
  )
  class(x) <- c(class, "block")
  return(x)
}

new_dual_block <- function(x, j, na.rm = TRUE, ..., class = character()) {
  K <- pm(x, t(x), na.rm = na.rm)
  new_block(
    x, j, na.rm, alpha = NULL, K = K, ..., class = c(class, "dual_block")
  )
}

new_primal_regularized_block <- function(x, j, tau, ...) {
  new_block(x, j, tau = tau, M = NULL, ..., class = "primal_regularized_block")
}

new_dual_regularized_block <- function(x, j, tau, ...) {
  new_dual_block(
    x, j, tau = tau, M = NULL, ..., class = "dual_regularized_block"
  )
}

new_sparse_block <- function(x, j, sparsity, tol = 1e-08, ...) {
  const <- sqrt(NCOL(x)) * sparsity
  new_block(
    x, j, sparsity = sparsity, const = const,
    tol = tol, ..., class = "sparse_block"
  )
}

new_ac_block <- function(x, j, tau, confounders, penalty_coef, algo, ...) {
  new_block(x, j, tau = tau, M = NULL, M_inv = NULL, B = NULL, 
            f = NULL, sqrt_M = NULL, sqrt_M_inv = NULL, 
            mu = NULL, f_left = NULL, f_right = NULL, 
            algo = algo, d = NULL, Q = NULL, 
            e_QM = NULL, a_MQ = NULL, e = NULL, h_tilde_QMX = NULL,
            h = NULL, h_tilde = NULL, eigen_val_Bplus1 = NULL, eigen_vec_Bplus1 = NULL,
            confounders = confounders, penalty_coef = penalty_coef, ..., 
            class = "ac_block")
}

### Utility method to choose the adequate class
create_block <- function(x, j, bias, na.rm, tau, sparsity, tol, confounders, penalty_coef, algo) {
  if (sparsity < 1) {
    res <- new_sparse_block(x, j, sparsity, tol, bias = bias, na.rm = na.rm)
  } else if (!is.null(confounders) && (penalty_coef != 0)) {
    res <- new_ac_block(x, j, tau, confounders = confounders, penalty_coef = penalty_coef, algo = algo)
  } else if (NROW(x) > NCOL(x)) {
    if (tau < 1) {
      res <- new_primal_regularized_block(x, j, tau, bias = bias, na.rm = na.rm)
    } else {
      res <- new_block(x, j, bias = bias, na.rm = na.rm)
    }
  } else {
    if (tau < 1) {
      res <- new_dual_regularized_block(x, j, tau, bias = bias, na.rm = na.rm)
    } else {
      res <- new_dual_block(x, j, bias = bias, na.rm = na.rm)
    }
  }
  return(res)
}
