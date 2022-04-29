#' # print.rgcca
#'''
data(Russett)
X_agric <- as.matrix(Russett[, c("gini", "farm", "rent")])
X_ind <- as.matrix(Russett[, c("gnpr", "labo")])
X_polit <- as.matrix(Russett[, c("demostab", "dictator")])
blocks <- list(X_agric, X_ind, X_polit)

test_that("print.rgcca", {
  local_edition(3)
  expect_snapshot({
    res <- rgcca(blocks, ncomp = 1, tau = c(1, 0.7, 0.9))
    print(res)
  })

  expect_snapshot({
    tau <- matrix(c(1, 1, 1, 0, 0, 0), nrow = 2, byrow = TRUE)
    res <- rgcca(blocks, ncomp = 2, tau = tau)
    print(res)
  })

  expect_snapshot({
    res <- rgcca(blocks, ncomp = 1, sparsity = c(0.7, 1, 0.9), method = "sgcca")
    print(res)
  })

  expect_snapshot({
    sparsity <- matrix(c(0.8, 1, 0.9, 1, 0.9, 0.8), nrow = 2, byrow = TRUE)
    res <- rgcca(blocks, ncomp = 2, sparsity = sparsity, method = "sgcca")
    print(res)
  })
})
