#' # print.permutation
#'''
set.seed(0)
data(Russett)
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)

test_that("print.permutation", {
  local_edition(3)
  expect_snapshot({
    res <- rgcca_permutation(blocks,
      par_type = "tau", par_length = 2,
      n_perms = 5, n_cores = 1, verbose = FALSE
    )
    print(res)
  })

  expect_snapshot({
    blocks2 <- c(blocks, blocks)
    names(blocks2) <- NULL
    res <- rgcca_permutation(blocks2,
      par_type = "ncomp", par_length = 2,
      n_perms = 2, n_cores = 1, verbose = FALSE
    )
    print(res)
  })
})
