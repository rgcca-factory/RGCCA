#' # summary.rgcca
#'''
data(Russett)
X_agric <- as.matrix(Russett[, c("gini", "farm", "rent")])
X_ind <- as.matrix(Russett[, c("gnpr", "labo")])
X_polit <- as.matrix(Russett[, c("demostab", "dictator")])
blocks <- list(X_agric, X_ind, X_polit)
blocks2 <- list(X_agric, X_ind, as.factor(Russett[, "demostab"]))

test_that("summary.rgcca produces the expected text", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  local_edition(3)
  expect_snapshot({
    res <- rgcca(blocks, ncomp = 1, tau = c(1, 0.7, 0.9))
    summary(res)
  })
})

test_that("summary.rgcca produces the expected text 2", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  local_edition(3)
  expect_snapshot({
    tau <- matrix(c(1, 1, 1, 0, 0, 0), nrow = 2, byrow = TRUE)
    res <- rgcca(blocks, ncomp = 2, tau = tau)
    summary(res)
  })
})

test_that("summary.rgcca produces the expected text 3", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  local_edition(3)
  expect_snapshot({
    res <- rgcca(blocks, ncomp = 1, sparsity = c(0.7, 1, 0.9), method = "sgcca")
    summary(res)
  })
})

test_that("summary.rgcca produces the expected text 4", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  local_edition(3)
  expect_snapshot({
    sparsity <- matrix(c(0.8, 1, 0.9, 1, 0.9, 0.8), nrow = 2, byrow = TRUE)
    res <- rgcca(blocks, ncomp = 2, sparsity = sparsity, method = "sgcca")
    summary(res)
  })
})

test_that("summary.rgcca produces the expected text 5", {
  skip_if_not(as.logical(Sys.getenv("TEST_SNAPSHOTS")))
  local_edition(3)
  expect_snapshot({
    res <- rgcca(blocks2, ncomp = 1, sparsity = c(0.7, 1, 0.9), response = 3)
    summary(res)
  })
})
