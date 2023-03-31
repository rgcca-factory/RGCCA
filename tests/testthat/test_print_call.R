#' # print_call
#'''
data(Russett)
X_agric <- as.matrix(Russett[, c("gini", "farm", "rent")])
X_ind <- as.matrix(Russett[, c("gnpr", "labo")])
X_polit <- as.matrix(Russett[, c("demostab", "dictator")])
blocks <- list(X_agric, X_ind, X_polit)

test_that("print_call prints the expected text", {
  local_edition(3)
  expect_snapshot({
    res <- rgcca(blocks,
      ncomp = 1, tau = 1, scheme = "horst",
      connection = 1 - diag(3), scale = FALSE,
      scale_block = "lambda1", superblock = FALSE,
      response = NULL, NA_method = "na.omit"
    )
    print_call(res$call)
  })
})

test_that("print_call prints the expected text 2", {
  local_edition(3)
  expect_snapshot({
    res <- rgcca(blocks,
      ncomp = 2, tau = 1, scheme = "centroid", scale = TRUE,
      scale_block = TRUE, superblock = FALSE,
      response = 3, NA_method = "na.ignore"
    )
    print_call(res$call)
  })
})
