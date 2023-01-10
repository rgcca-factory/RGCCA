### Test format_bootstrap_list
set.seed(0)
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

extract_from_list <- function(W, block, type, comp, boot) {
  type <- ifelse(type == "weights", "W", "L")
  return(unname(W[[boot]][[type]][[block]][, comp]))
}
extract_from_df <- function(df, block_, type_, comp_, boot_) {
  return(subset(
    df, block == block_ & type == type_ & boot == boot_ & comp == comp_
  )$value)
}

test_that("format_bootstrap_list creates a data.frame with bootstrap results", {
  res <- format_bootstrap_list(W, fit_rgcca)
  expect_true(nrow(res) == 2 * n_boot * sum(
    vapply(fit_rgcca$a, length, FUN.VALUE = 1L)
  ))
  n <- 1
  for (i in seq_len(10)) {
    boot <- sample(seq(n_boot), 1)
    type <- sample(c("weights", "loadings"), 1)
    block <- sample(names(blocks), 1)
    expect_equal(
      extract_from_list(W, block, type, n, boot),
      extract_from_df(res, block, type, n, boot)
    )
  }

  n <- 2
  for (i in seq_len(10)) {
    boot <- sample(seq(n_boot), 1)
    type <- sample(c("weights", "loadings"), 1)
    block <- sample(names(blocks)[-1], 1)
    expect_equal(
      extract_from_list(W, block, type, n, boot),
      extract_from_df(res, block, type, n, boot)
    )
  }
})
