data(Russett)
X_agric <- as.matrix(Russett[, c("gini", "farm", "rent")])
X_ind <- as.matrix(Russett[, c("gnpr", "labo")])
X_polit <- as.matrix(Russett[, c("demostab", "dictator")])
blocks <- list(X_agric, X_ind, X_polit)
blocks <- lapply(blocks, scale)
blocks[[4]] <- do.call(cbind, blocks)

C1 <- 1 - diag(3)
C2 <- matrix(c(
  0, 1, 0,
  1, 0, 1,
  0, 1, 0
), nrow = 3, byrow = TRUE)
C3 <- matrix(c(
  0, 0, 0, 1,
  0, 0, 0, 1,
  0, 0, 0, 1,
  1, 1, 1, 0
), nrow = 4, byrow = TRUE)

test_that("rgccad produces cumulated AVE that are below 1", {
  res <- rgccad(blocks[-4], connection = C1, ncomp = rep(2, 3), verbose = FALSE)
  expect_true(all(unlist(lapply(res$AVE$AVE_X_cor, sum)) <= 1))

  res <- rgccad(blocks[-4],
    connection = C2, ncomp = rep(3, 3),
    response = 2, verbose = FALSE
  )
  expect_true(all(unlist(lapply(res$AVE$AVE_X_cor, sum)) <= 1))

  res <- rgccad(blocks,
    connection = C3, ncomp = rep(6, 4),
    superblock = TRUE, verbose = FALSE
  )
  expect_true(all(unlist(lapply(res$AVE$AVE_X_cor, sum)) <= 1))
})

test_that("rgccad returns equal AVE and corrected AVE if components are
          not correlated", {
  res <- rgccad(blocks[-4], connection = C1, ncomp = rep(2, 3), verbose = FALSE)
  expect_true(all.equal(res$AVE$AVE_X_cor, res$AVE$AVE_X))
})
