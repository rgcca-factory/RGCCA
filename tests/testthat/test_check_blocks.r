data(Russett)
X_agric <- as.matrix(Russett[, c("gini", "farm", "rent")])
X_ind <- as.matrix(Russett[, c("gnpr", "labo")])
X_polit <- as.matrix(Russett[, "demostab"])
X_quali <- colnames(Russett)[9:11][apply(Russett[, 9:11], 1, which.max)]

test_that("check_blocks returns a list of blocks", {
  expect_true(is.list(check_blocks(X_agric)))
  expect_true(is.list(check_blocks(list(X_agric, X_ind))))
})

test_that("check_blocks raises an error if a list of blocks cannot be
          created", {
  expect_error(check_blocks(NA), "blocks must be a list.", fixed = TRUE)
})

test_that("check_blocks returns a list of matrices", {
  blocks <- list(as.matrix(X_agric), as.data.frame(X_ind), as.vector(X_polit))
  expect_true(
    all(vapply(check_blocks(blocks), is.matrix, FUN.VALUE = logical(1)))
  )
})

test_that("check_blocks returns an error if a block is qualitative and
          is not the response block", {
  blocks <- list(X_agric, X_ind, X_quali)
  expect_error(
    check_blocks(blocks), "unsupported qualitative block.",
    fixed = TRUE
  )
  expect_error(check_blocks(blocks, response = 3), NA)
})

test_that("check_blocks returns an error if a block has multiple variates
          with at least a qualitative one", {
  blocks <- list(X_agric, X_ind, cbind(X_quali, X_ind))
  expect_error(
    check_blocks(blocks, response = 3),
    "unsupported multivariate qualitative block.",
    fixed = TRUE
  )
})

test_that("check_blocks renames blocks if names are missing", {
  expect_equal(names(check_blocks(list(X_agric, X_ind))), c("block1", "block2"))
  expect_equal(
    names(check_blocks(list(agriculture = X_agric, X_ind))),
    c("agriculture", "block2")
  )
  expect_equal(
    names(check_blocks(list(agriculture = X_agric, industry = X_ind))),
    c("agriculture", "industry")
  )
})

test_that("check_blocks add colnames with blocks with no colnames", {
  blocks <- list(agri = X_agric, polit = X_polit)
  expect_equal(colnames(check_blocks(blocks)[[2]]), "polit")
  expect_equal(colnames(check_blocks(blocks)[[1]]), colnames(X_agric))
  colnames(blocks[[1]]) <- NULL
  expect_equal(
    colnames(check_blocks(blocks)[[1]]), paste0("agri_", seq_len(NCOL(X_agric)))
  )
})

test_that("check_blocks add prefixes to avoid duplicated colnames", {
  blocks <- list(agri = X_agric, ind = X_ind)
  expect_message(check_blocks(blocks, quiet = FALSE), NA)
  colnames(blocks[[2]]) <- c("gini", "labo")
  expect_message(
    check_blocks(blocks, quiet = FALSE),
    "Duplicated dimnames across blocks are modified to avoid confusion.",
    fixed = TRUE
  )
  blocks2 <- check_blocks(blocks, quiet = FALSE)
  expect_equal(
    colnames(blocks2[[1]]), paste("agri", colnames(blocks[[1]]), sep = "_")
  )
  expect_equal(
    colnames(blocks2[[2]]), paste("ind", colnames(blocks[[2]]), sep = "_")
  )
})

test_that("check_blocks raises an error if there are duplicated rownames", {
  blocks <- list(rbind(X_agric, X_agric), rbind(X_ind, X_ind))
  expect_error(
    check_blocks(blocks), "blocks have duplicated names on dimension 1.",
    fixed = TRUE
  )
})

test_that("check_blocks creates rownames if no block has rownames", {
  expect_equal(
    rownames(check_blocks(X_polit)[[1]]), paste0("S", seq_along(X_polit))
  )
  expect_message(check_blocks(list(polit = X_polit), quiet = FALSE), NA)
})

test_that("check_blocks add rownames if a block lacks rownames and other
          rownames are compatible", {
  blocks <- list(agric = X_agric, ind = X_ind, polit = X_polit)
  expect_equal(rownames(check_blocks(blocks)[[3]]), rownames(blocks[[1]]))
  expect_message(check_blocks(blocks, quiet = FALSE), NA)
})

test_that("check_blocks raises an error if a block lacks rownames and other
          rownames do not match", {
  blocks <- list(agric = X_agric, ind = X_ind, polit = X_polit)
  rownames(blocks[[2]]) <- rownames(blocks[[1]])[c(2, 1, seq(3, 47))]
  expect_error(check_blocks(blocks), paste0(
    "some blocks are missing names on dimension 1, and the other blocks' ",
    "names on dimension 1 are not consistent."
  ), fixed = TRUE)
})

test_that("check_blocks returns blocks with the same rownames in the same
          order", {
  blocks <- list(agric = X_agric, ind = X_ind)
  rownames(blocks[[2]]) <- c("xxx", rownames(blocks[[1]])[seq(2, 47)])
  blocks2 <- check_blocks(blocks)
  expect_equal(rownames(blocks2[[1]]), rownames(blocks2[[2]]))
  expect_equal(rownames(blocks2[[1]]), union(rownames(blocks[[1]]), "xxx"))
})

### Make sure check_blocks works for arrays
n <- nrow(X_agric)
array_block <- array(rnorm(n * 12 * 15), dim = c(n, 12, 15))

test_that("check_blocks still works when a block is an array", {
  blocks <- list(agric = X_agric, ind = X_ind, array = array_block)
  blocks <- check_blocks(blocks)
  expect_equal(rownames(blocks[[1]]), rownames(blocks[[3]]))
  expect_equal(
    dimnames(check_blocks(blocks)[[3]])[[2]],
    paste0("array_2_", seq_len(dim(array_block)[2]))
  )
  expect_equal(
    dimnames(check_blocks(blocks)[[3]])[[3]],
    paste0("array_3_", seq_len(dim(array_block)[3]))
  )

  blocks <- list(agric = X_agric, ind = X_ind, array = array_block)
  dimnames(blocks[[3]])[c(2, 3)] <- list(
    seq_len(dim(array_block)[2]), seq_len(dim(array_block)[3])
  )
  expect_message(
    check_blocks(blocks, quiet = FALSE),
    "Duplicated dimnames are modified to avoid confusion.", fixed = TRUE
  )

  dimnames(blocks[[3]])[c(2, 3)] <- list(
    paste0("2", seq_len(dim(array_block)[2])),
    paste0("3", seq_len(dim(array_block)[3]))
  )
  blocks2 <- check_blocks(blocks)
  expect_equal(dimnames(blocks2[[3]])[-1], dimnames(blocks[[3]])[-1])

  rownames(blocks[[3]]) <- sample(
    rownames(blocks[[1]]), size = n, replace = FALSE
  )
  blocks2 <- check_blocks(blocks)
  expect_equal(rownames(blocks2[[1]]), rownames(blocks2[[3]]))
  expect_equal(
    blocks2[[3]][rownames(blocks2[[3]]), , ],
    blocks[[3]][rownames(blocks2[[3]]), , ]
  )

  rownames(blocks[[3]])[1] <- "xxx"
  blocks2 <- check_blocks(blocks)
  expect_equal(rownames(blocks2[[1]]), union(rownames(blocks[[1]]), "xxx"))

  rownames(blocks[[2]]) <- NULL
  expect_error(check_blocks(blocks), paste0(
    "some blocks are missing names on dimension 1, and the other blocks' ",
    "names on dimension 1 are not consistent."
  ), fixed = TRUE)
})
