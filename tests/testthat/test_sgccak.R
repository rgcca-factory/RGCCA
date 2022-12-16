set.seed(0)
blocks <- list(
  matrix(rnorm(100 * 41), nrow = 100),
  matrix(rnorm(100 * 8), nrow = 100),
  matrix(rnorm(100 * 24), nrow = 100)
)
connection <- 1 - diag(3)
blocks <- scaling(blocks, scale = FALSE, scale_block = FALSE)
fit_sgccak <- sgccak(
  blocks, connection,
  sparsity = c(0.2, 0.5, 0.3),
  scheme = "factorial", tol = 1e-8
)

test_that("criterion is increasing", {
  n_iter <- length(fit_sgccak$crit)
  crit_old <- fit_sgccak$crit[1:(n_iter - 1)]
  crit_next <- fit_sgccak$crit[2:n_iter]
  expect_true(all(crit_next - crit_old >= 0))
})
