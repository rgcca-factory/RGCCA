set.seed(0)

verify_norm_constraint <- function(a, const, tol) {
  for (j in seq_along(a)) {
    expect_lte(norm(a[[j]], type = "2"), 1 + tol)
    expect_lte(sum(abs(a[[j]])), const[j])
  }
}

verify <- function(A, sparsity, init = "svd", tol = 1e-14) {
  pjs <- vapply(A, ncol, FUN.VALUE = integer(1L))
  const <- sparsity * pjs
  init_object <- sgcca_init(
    A, init, TRUE, TRUE, sparsity, response = NULL, disjunction = NULL
  )
  verify_norm_constraint(init_object$a, const, tol)
}

data(Russett)
X_agric <- as.matrix(Russett[, c("gini", "farm", "rent")])
X_ind <- as.matrix(Russett[, c("gnpr", "labo")])
X_polit <- as.matrix(Russett[, c("demostab", "dictator")])
A <- list(X_agric, X_ind, X_polit)
A <- scaling(A, scale = TRUE, bias = TRUE, scale_block = TRUE)

test_that("sgcca_init generates vectors a that satisfy the norm constraints", {
  verify(A, sparsity = c(1, 1, 1), init = "svd")
  verify(A, sparsity = c(1, 1, 1), init = "random")

  verify(A, sparsity = c(0.71, 0.9, 0.87), init = "svd")
  verify(A, sparsity = c(0.71, 0.9, 0.87), init = "random")
})
