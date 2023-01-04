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

test_that("format_bootstrap_list creates a data.frame with bootstrap results", {
  res <- format_bootstrap_list(W, fit_rgcca)
  expect_true(nrow(res) == 2 * n_boot * sum(
    vapply(fit_rgcca$a, length, FUN.VALUE = 1L)
  ))
  n <- 1
  b <- 17
  w <- res %>%
    dplyr::filter(
      comp == n, boot == b, type == "weights", block == "agriculture"
    )
  expect_equal(w$value, unname(W[[b]][[1]][[1]][, n]))
})
