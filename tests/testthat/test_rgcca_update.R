set.seed(0)

g <- function(x) x^2
dg <- Deriv::Deriv(g, env = parent.frame())

verify_norm_constraint <- function(res, M) {
  for (j in seq_along(M)) {
    expect_equal(drop(t(res$a[[j]]) %*% M[[j]] %*% res$a[[j]]), 1)
  }
}

verify <- function(A, tau, init, C, dg, tol = 1e-12) {
  J <- length(A)
  n <- nrow(A[[1]])
  pjs <- vapply(A, ncol, FUN.VALUE = integer(1L))
  which.primal <- which((n >= pjs) == 1)
  which.dual <- which((n < pjs) == 1)
  M <- lapply(1:J, function(j) {
    tau[j] * diag(pjs[j]) + ((1 - tau[j])) * 1 / n *
      (pm(t(A[[j]]), A[[j]], na.rm = TRUE))
  })
  tmp <- rgcca_init(
    A, init, TRUE, TRUE, tau, pjs, which.primal,
    which.dual, J, n
  )
  crit_old <- sum(C * g(cov2(tmp$Y, bias = TRUE)))
  tmp <- rgcca_update(
    A, tmp$a, tmp$alpha, tmp$Y, tmp$M, tmp$K, tmp$Minv,
    TRUE, TRUE, tau, which.primal, which.dual, J, n, dg, C
  )
  verify_norm_constraint(tmp, M)
  crit <- sum(C * g(cov2(tmp$Y, bias = TRUE)))
  expect_true(crit - crit_old > -tol)
}

### Test primal case
data(Russett)
X_agric <- as.matrix(Russett[, c("gini", "farm", "rent")])
X_ind <- as.matrix(Russett[, c("gnpr", "labo")])
X_polit <- as.matrix(Russett[, c("demostab", "dictator")])
A <- list(X_agric, X_ind, X_polit)
A <- scaling(A, scale = T, bias = T, scale_block = T)
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
A <- scaling(A, scale = T, bias = T, scale_block = T)
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
