set.seed(0)

g <- function(x) x^2
dg <- Deriv::Deriv(g, env = parent.frame())

verify_norm_constraint <- function(a, const, tol) {
  for (j in seq_along(a)) {
    expect_lte(norm(a[[j]], type = "2"), 1 + tol)
    expect_lte(sum(abs(a[[j]])), const[j])
  }
}

verify <- function(A, sparsity, init, C, dg, tol = 1e-14) {
  J <- length(A)
  n <- nrow(A[[1]])
  pjs <- vapply(A, ncol, FUN.VALUE = integer(1L))
  const <- sparsity * pjs
  tmp <- sgcca_init(A, init, TRUE, TRUE, const, pjs, J, n)
  crit_old <- sum(C * g(cov2(tmp$Y, bias = TRUE)))
  tmp <- sgcca_update(A, tmp$a, tmp$Y, TRUE, TRUE, const, J, n, dg, C)
  verify_norm_constraint(tmp$a, const, tol)
  crit <- sum(C * g(cov2(tmp$Y, bias = TRUE)))
  expect_true(crit - crit_old > -tol)
}

data(Russett)
X_agric <- as.matrix(Russett[, c("gini", "farm", "rent")])
X_ind <- as.matrix(Russett[, c("gnpr", "labo")])
X_polit <- as.matrix(Russett[, c("demostab", "dictator")])
A <- list(X_agric, X_ind, X_polit)
A <- scaling(A, scale = T, bias = T, scale_block = T)
C <- 1 - diag(3)

test_that("sgcca_update generates vectors a that satisfy the norm
          constraints", {
  verify(A, sparsity = c(1, 1, 1), init = "svd", C, dg)
  verify(A, sparsity = c(1, 1, 1), init = "random", C, dg)

  verify(A, sparsity = c(0.71, 0.9, 0.87), init = "svd", C, dg)
  verify(A, sparsity = c(0.71, 0.9, 0.87), init = "random", C, dg)
})
