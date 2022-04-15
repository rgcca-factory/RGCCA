data(Russett)
X_agric <- as.matrix(Russett[, c("gini", "farm", "rent")])
X_ind <- as.matrix(Russett[, c("gnpr", "labo")])
X_polit <- as.matrix(Russett[, c("demostab")])

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

  # Check for messages as well
  expect_message(
    check_blocks(list(agriculture = X_agric, X_ind), quiet = FALSE),
    "Missing block names are automatically labeled.",
    fixed = TRUE
  )
  expect_message(
    check_blocks(list(agriculture = X_agric, industry = X_ind), quiet = FALSE),
    NA
  )
})

test_that("check_blocks add colnames with blocks with no colnames", {
  blocks <- list(agri = X_agric, polit = X_polit)
  expect_equal(colnames(check_blocks(blocks)[[2]]), c("polit"))
  expect_equal(colnames(check_blocks(blocks)[[1]]), colnames(X_agric))
  colnames(blocks[[1]]) <- NULL
  expect_equal(
    colnames(check_blocks(blocks)[[1]]), paste0("V1_", seq(NCOL(X_agric)))
  )
  expect_message(check_blocks(blocks, quiet = FALSE),
    "Missing colnames are automatically labeled.",
    fixed = TRUE
  )
  expect_message(check_blocks(list(agric = X_agric), quiet = FALSE), NA)
})

test_that("check_blocks add prefixes to avoid duplicated colnames", {
  blocks <- list(agri = X_agric, ind = X_ind)
  expect_message(check_blocks(blocks, quiet = FALSE), NA)
  colnames(blocks[[2]]) <- c("gini", "labo")
  expect_message(
    check_blocks(blocks, quiet = FALSE),
    "Duplicated colnames are modified to avoid confusion.",
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
    check_blocks(blocks), "blocks have duplicated rownames.",
    fixed = TRUE
  )
})

test_that("check_blocks creates rownames if no block has rownames", {
  expect_equal(
    rownames(check_blocks(X_polit)[[1]]), paste0("S", seq_along(X_polit))
  )
  expect_message(
    check_blocks(list(polit = X_polit), quiet = FALSE),
    "Missing rownames are automatically labeled.",
    fixed = TRUE
  )
  expect_message(check_blocks(list(agric = X_agric), quiet = FALSE), NA)
  expect_error(
    check_blocks(X_polit, allow_unnames = FALSE), "blocks must have rownames.",
    fixed = TRUE
  )
})

test_that("check_blocks add rownames if a block lacks rownames and other
          rownames are compatible", {
  blocks <- list(agric = X_agric, ind = X_ind, polit = X_polit)
  expect_equal(rownames(check_blocks(blocks)[[3]]), rownames(blocks[[1]]))
  expect_message(
    check_blocks(blocks, quiet = FALSE),
    "Missing rownames are automatically labeled.",
    fixed = TRUE
  )
  expect_message(check_blocks(blocks[-3], quiet = FALSE), NA)
})

test_that("check_blocks raises an error if a block lacks rownames and other
          rownames do not match", {
  blocks <- list(agric = X_agric, ind = X_ind, polit = X_polit)
  rownames(blocks[[2]]) <- rownames(blocks[[1]])[c(2, 1, seq(3, 47))]
  expect_error(check_blocks(blocks), paste0(
    "some blocks are missing rownames, and the other blocks' ",
    "rownames are not consistent."
  ), fixed = TRUE)
})

test_that("check_blocks allows only blocks with the same rownames if
          add_NAlines is FALSE", {
  blocks <- list(agric = X_agric, ind = X_ind)
  rownames(blocks[[2]]) <- c("xxx", rownames(blocks[[1]])[seq(2, 47)])
  expect_error(check_blocks(blocks, add_NAlines = FALSE),
    "blocks must have the same rownames",
    fixed = TRUE
  )
  rownames(blocks[[2]]) <- rownames(blocks[[1]])[c(2, 1, seq(3, 47))]
  expect_error(check_blocks(blocks, add_NAlines = FALSE), NA)
})

test_that("check_blocks returns blocks with the same rownames in the same
          order", {
  blocks <- list(agric = X_agric, ind = X_ind)
  rownames(blocks[[2]]) <- c("xxx", rownames(blocks[[1]])[seq(2, 47)])
  blocks2 <- check_blocks(blocks, add_NAlines = TRUE)
  expect_equal(rownames(blocks2[[1]]), rownames(blocks2[[2]]))
  expect_equal(rownames(blocks2[[1]]), union(rownames(blocks[[1]]), "xxx"))
})

test_that("check_blocks removes null variance columns if init is TRUE", {
  blocks <- list(agric = cbind(X_agric, 0))
  colnames(blocks[[1]])[4] <- "null_variance_col"
  expect_equal(check_blocks(blocks, init = TRUE)[[1]], X_agric)
  expect_equal(check_blocks(blocks, init = FALSE), blocks)
})
