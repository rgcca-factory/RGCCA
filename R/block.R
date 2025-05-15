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

new_sparse_block <- function(x, j, sparsity, tol = 1e-08, ...,
  class=character()) {
  const <- sqrt(NCOL(x)) * sparsity
  new_block(
    x, j, sparsity = sparsity, const = const,
    tol = tol, ..., class = c(class, "sparse_block")
  )
}

new_graphnet_block <- function(x, j, sparsity, lambda, graph_laplacians, 
  tol = 1e-08, ...) {
  new_sparse_block(x, j, sparsity = sparsity, tol = tol,
    lambda = lambda, graph_laplacians = graph_laplacians, ...,
    class = "graphnet_block")
}

### Utility method to choose the adequate class
create_block <- function(x, j, bias, na.rm, tau, sparsity, lambda,
  graph_laplacian, tol) {
  if (!is.null(graph_laplacian)) {
    res <- new_graphnet_block(x, j, sparsity, lambda, graph_laplacian, tol,
      bias = bias, na.rm = na.rm)
  } else if (sparsity < 1) {
    res <- new_sparse_block(x, j, sparsity, tol, bias = bias, na.rm = na.rm)
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
