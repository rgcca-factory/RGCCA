#' # print.cval
#'''
set.seed(0)
data(Russett)
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)

test_that("print_cval", {
  local_edition(3)
  expect_snapshot({
    res <- rgcca_cv(blocks,
      response = 3, method = "rgcca", par_type = "tau",
      par_value = c(0, 0.2, 0.3), n_run = 1, n_cores = 1,
      par_length = 2, verbose = FALSE
    )
    print(res)
  })

  expect_snapshot({
    res <- rgcca_cv(blocks,
      validation = "loo",
      response = 3, method = "sgcca", par_type = "sparsity",
      n_run = 1, n_cores = 1, par_length = 2, verbose = FALSE
    )
    print(res, bars = "sd")
  })

  expect_snapshot({
    blocks2 <- c(blocks, blocks)
    names(blocks2) <- NULL
    res <- rgcca_cv(blocks2,
      validation = "kfold", k = 2,
      response = 3, method = "rgcca", par_type = "tau",
      n_run = 1, n_cores = 1, par_length = 2, verbose = FALSE
    )
    print(res, bars = "stderr")
  })
})
