set.seed(0)

verify_norm_constraint <- function(a, M) {
  for (j in seq_along(a)) {
    expect_equal(drop(t(a[[j]]) %*% M[[j]] %*% a[[j]]), 1)
  }
}

verify <- function(A, tau, init = "svd", tol = 1e-12) {
  n <- nrow(A[[1]])
  M <- lapply(seq_along(A), function(j) {
    tau[j] * diag(NCOL(A[[j]])) + ((1 - tau[j])) * 1 / n *
      (pm(t(A[[j]]), A[[j]], na.rm = TRUE))
  })
  init_object <- rgcca_init(A, init, TRUE, TRUE, tau)
  verify_norm_constraint(init_object$a, M)
}

### Test primal case
data(Russett)
X_agric <- as.matrix(Russett[, c("gini", "farm", "rent")])
X_ind <- as.matrix(Russett[, c("gnpr", "labo")])
X_polit <- as.matrix(Russett[, c("demostab", "dictator")])
A <- list(X_agric, X_ind, X_polit)
A <- scaling(A, scale = T, bias = T, scale_block = T)

test_that("rgcca_init generates vectors a that satisfy the norm constraints
          in the primal case", {
  verify(A, tau = c(1, 1, 1), init = "svd")
  verify(A, tau = c(1, 1, 1), init = "random")

  verify(A, tau = c(0.5, 0.5, 0.5), init = "svd")
  verify(A, tau = c(0.5, 0.5, 0.5), init = "random")

  verify(A, tau = c(0, 0, 0), init = "svd")
  verify(A, tau = c(0, 0, 0), init = "random")
})

### Test dual case
A <- list(
  matrix(rnorm(48), nrow = 4, ncol = 12),
  matrix(rnorm(40), nrow = 4, ncol = 10)
)
A <- scaling(A, scale = T, bias = T, scale_block = T)

test_that("rgcca_init generates vectors a that satisfy the norm constraints
          in the dual case", {
  verify(A, tau = c(1, 1, 1), init = "svd")
  verify(A, tau = c(1, 1, 1), init = "random")

  verify(A, tau = c(0.5, 0.5, 0.5), init = "svd")
  verify(A, tau = c(0.5, 0.5, 0.5), init = "random")

  verify(A, tau = c(0, 0, 0), init = "svd")
  verify(A, tau = c(0, 0, 0), init = "random")
})
