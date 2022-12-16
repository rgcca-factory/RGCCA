### Test format_bootstrap_list
data(Russett)
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)

n_boot <- 20
ncomp <- c(1, 2, 2)
fit_rgcca <- rgcca(blocks, ncomp = ncomp)
W <- lapply(seq(n_boot), function(i) {
  rgcca_bootstrap_k(fit_rgcca, sample(seq_len(nrow(blocks[[1]]))))
})

test_that("format_bootstrap_list reorders raw bootstrap results", {
  res <- format_bootstrap_list(W, fit_rgcca, n_boot, n = 1)
  expect_equal(length(res), max(ncomp))
  expect_true(
    all(vapply(res, length, FUN.VALUE = integer(1)) == length(blocks))
  )
  expect_true(
    all(vapply(res, function(x) ncol(x[[1]]), FUN.VALUE = integer(1)) == n_boot)
  )
  n <- 1
  b <- 17
  expect_equal(res[[n]][[1]][, b], W[[b]][[1]][[1]][, n])
})
