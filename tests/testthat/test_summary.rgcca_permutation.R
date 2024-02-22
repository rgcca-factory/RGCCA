#' # summary.rgcca_permutation
#'''
set.seed(0)
data(Russett)
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)

test_that("summary.rgcca_permutation prints the expected text", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  local_edition(3)
  expect_snapshot({
    res <- rgcca_permutation(blocks,
      par_type = "tau", par_length = 2,
      n_perms = 5, n_cores = 1, verbose = FALSE
    )
    summary(res)
  })
})

test_that("summary.rgcca_permutation prints the expected text 2", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  local_edition(3)
  expect_snapshot({
    blocks2 <- rep(blocks, 3)
    names(blocks2) <- NULL
    res <- rgcca_permutation(blocks2,
      par_type = "ncomp", par_length = 2,
      n_perms = 2, n_cores = 1, verbose = FALSE
    )
    summary(res)
  })
})
