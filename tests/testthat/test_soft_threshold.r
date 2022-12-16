set.seed(0)

### Test soft_threshold
verify_constraints <- function(x, s, tol = 1e-12) {
  expect_lte(norm(x, type = "2"), 1 + tol)
  expect_lte(sum(abs(x)), s + tol)
}

test_that("soft_threshold returns a vector that satisfies the constraints", {
  verify_constraints(soft_threshold(rep(0, 12), 1), 1)
  verify_constraints(soft_threshold(rnorm(42), 10), 10)
  verify_constraints(soft_threshold(rnorm(1500), 127), 127)
})
