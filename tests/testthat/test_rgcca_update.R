set.seed(0)

g <- function(x) x^2
dg <- Deriv::Deriv(g, env = parent.frame())

verify_norm_constraint <- function(a, M) {
  for (j in seq_along(a)) {
    expect_equal(drop(t(a[[j]]) %*% M[[j]] %*% a[[j]]), 1)
  }
}

verify <- function(A, tau, init, C, dg, tol = 1e-12) {
  # Initialize
  n <- nrow(A[[1]])
  M <- lapply(seq_along(A), function(j) {
    tau[j] * diag(NCOL(A[[j]])) + ((1 - tau[j])) * 1 / n *
      (pm(t(A[[j]]), A[[j]], na.rm = TRUE))
  })
  init_object <- rgcca_init(A, init, TRUE, TRUE, tau)
  a <- init_object$a
  Y <- init_object$Y
  crit_old <- sum(C * g(cov2(init_object$Y, bias = TRUE)))
  # Compute update
  update_object <- rgcca_update(A, TRUE, TRUE, tau, dg, C, a, Y, init_object)
  verify_norm_constraint(update_object$a, M)
  crit <- sum(C * g(cov2(update_object$Y, bias = TRUE)))
  expect_true(crit - crit_old > -tol)
}

### Test primal case
data(Russett)
X_agric <- as.matrix(Russett[, c("gini", "farm", "rent")])
X_ind <- as.matrix(Russett[, c("gnpr", "labo")])
X_polit <- as.matrix(Russett[, c("demostab", "dictator")])
A <- list(X_agric, X_ind, X_polit)
A <- scaling(A, scale = TRUE, bias = TRUE, scale_block = TRUE)
C <- 1 - diag(3)


test_that("rgcca_update generates vectors a that satisfy the norm constraints
          in the primal case and increase the criterion", {
  verify(A, tau = c(1, 1, 1), init = "svd", C, dg)
  verify(A, tau = c(1, 1, 1), init = "random", C, dg)

  verify(A, tau = c(0.5, 0.5, 0.5), init = "svd", C, dg)
  verify(A, tau = c(0.5, 0.5, 0.5), init = "random", C, dg)

  verify(A, tau = c(0, 0, 0), init = "svd", C, dg)
  verify(A, tau = c(0, 0, 0), init = "random", C, dg)
})

### Test dual case
A <- list(
  matrix(rnorm(48), nrow = 4, ncol = 12),
  matrix(rnorm(40), nrow = 4, ncol = 10)
)
A <- scaling(A, scale = TRUE, bias = TRUE, scale_block = TRUE)
C <- 1 - diag(2)

test_that("rgcca_update generates vectors a that satisfy the norm constraints
          in the dual case and increase the crition", {
  verify(A, tau = c(1, 1, 1), init = "svd", C, dg)
  verify(A, tau = c(1, 1, 1), init = "random", C, dg)

  verify(A, tau = c(0.5, 0.5, 0.5), init = "svd", C, dg)
  verify(A, tau = c(0.5, 0.5, 0.5), init = "random", C, dg)

  verify(A, tau = c(0, 0, 0), init = "svd", C, dg)
  verify(A, tau = c(0, 0, 0), init = "random", C, dg)
})
