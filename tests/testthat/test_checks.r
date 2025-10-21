# Load data
data(Russett)
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)

# Test check_blockx
test_that("check_blockx raises an error if x is greater than the number of
          blocks", {
  expect_error(check_blockx("x", 7, blocks),
    "x must be lower than the number of blocks, i.e. 3.",
    fixed = TRUE
  )
})
test_that("check_blockx raises an error for invalid x", {
  expect_error(check_blockx("x", -1, blocks))
  expect_error(check_blockx("x", c(1, 2, 3), blocks))
  expect_error(check_blockx("x", 1.7, blocks))
  expect_error(check_blockx("x", NA, blocks))
})
test_that("check_blockx passes and returns x when x is valid", {
  expect_equal(check_blockx("x", 2, blocks), 2)
})

# Test check_boolean
test_that("check_boolean raises an error if x contains NA", {
  expect_error(check_boolean("x", c(FALSE, NA)),
    "x must not be NA.",
    fixed = TRUE
  )
})
test_that("check_boolean raises an error if x is not logical", {
  expect_error(check_boolean("x", 0),
    "x must be TRUE or FALSE.",
    fixed = TRUE
  )
})
test_that("check_boolean raises an error if type is scalar and x not of
          length 1", {
  expect_error(check_boolean("x", c(TRUE, FALSE), type = "scalar"),
    "x must be of length 1.",
    fixed = TRUE
  )
})
test_that("check_boolean passes when x is valid", {
  expect_error(check_boolean(TRUE), NA)
  expect_error(check_boolean(c(TRUE, FALSE), type = "vector"), NA)
})

# Test check_colors
test_that("check_colors raises an error if colors contains unknown colors", {
  expect_error(check_colors(c("white", "blueen")),
    paste0(
      "Unrecognized colors. Colors must be in colors() ",
      "or a rgb character."
    ),
    fixed = TRUE
  )
  expect_error(check_colors(c("white", "#79eff16342")),
    paste0(
      "Unrecognized colors. Colors must be in colors() ",
      "or a rgb character."
    ),
    fixed = TRUE
  )
})
test_that("check_colors passes when colors is valid", {
  expect_error(check_colors(c("#79eff1", "white", "yellow", NA)), NA)
})

# Test check_compx
test_that("check_compx raises an error if x is greater than the number
          of components", {
  expect_error(check_compx("x", 7, c(3, 3, 3), 1), paste0(
    "not existing component. Trying to extract component 7",
    " for block 1 , but only 3 component(s) are available for ",
    "this block."
  ),
  fixed = TRUE
  )
})
test_that("check_compx raises an error for invalid x", {
  expect_error(check_compx("x", -1, c(3, 3, 3), 1))
  expect_error(check_compx("x", c(1, 2, 3), c(3, 3, 3), 1))
  expect_error(check_compx("x", 1.7, c(3, 3, 3), 1))
  expect_error(check_compx("x", NA, c(3, 3, 3), 1))
})
test_that("check_compx passes and returns x when x is valid", {
  expect_equal(check_compx("x", 2, c(3, 3, 3), 1), 2)
})

# Test check_connection
test_that("check_connection raises an error if C is not a matrix", {
  expect_error(check_connection(print, blocks),
    "connection matrix C must be a matrix.",
    fixed = TRUE
  )
})
test_that("check_connection raises an error if C is not symmetric", {
  expect_error(check_connection(matrix(1:4, 2, 2), blocks),
    "connection matrix C must be symmetric.",
    fixed = TRUE
  )
})
test_that("check_connection raises an error if C contains NA values", {
  C <- diag(2)
  C[1, 1] <- NA
  expect_error(check_connection(C, blocks),
    "connection matrix C must not contain NA values.",
    fixed = TRUE
  )
})
test_that("check_connection raises an error if C contains values outside
          [0, 1]", {
  expect_error(check_connection(diag(2) * 2, blocks),
    "connection matrix C must contain numbers between 0 and 1.",
    fixed = TRUE
  )
})
test_that("check_connection raises an error if C is the null matrix", {
  expect_error(check_connection(matrix(0, 2, 2), blocks),
    "connection matrix C must not contain only 0.",
    fixed = TRUE
  )
})
test_that("check_connection raises an error if the dimensions of C do not
          match the number of blocks", {
  expect_error(check_connection(diag(2), blocks),
    paste0(
      "connection matrix must have the same number of ",
      "columns (actually 2) as the number of blocks (3)."
    ),
    fixed = TRUE
  )
})
test_that("check_connection raises an error if the dimnames of C do not match
          block names", {
  C <- diag(3)
  rownames(C) <- colnames(C) <- paste0("V", 1:3)
  expect_error(check_connection(C, blocks),
    paste0(
      "connection matrix C must have the rownames and the ",
      "colnames that match with the names of the blocks."
    ),
    fixed = TRUE
  )
})
test_that("check_connection passes and returns C when C is valid", {
  C <- diag(3)
  rownames(C) <- colnames(C) <- names(blocks)
  expect_equal(check_connection(diag(3), blocks), C)
})

# Test check_integer
test_that("check_integer raises an error if x is not numeric", {
  expect_error(check_integer("x", "toto"),
    "x must be numeric.",
    fixed = TRUE
  )
})

test_that("check_integer raises an error if x contains NA", {
  expect_error(check_integer("x", c(42, NA)),
    "x must not be NA.",
    fixed = TRUE
  )
})

test_that("check_integer raises an error if type is scalar and x not of
          length 1", {
  expect_error(check_integer("x", c(42, 7), type = "scalar"),
    "x must be of length 1.",
    fixed = TRUE
  )
})

test_that("check_integer raises an error any element of x is a float but float
          is false", {
  expect_error(check_integer("x", c(1, 1.7, 2), type = "vector", float = FALSE),
    "x must be an integer.",
    fixed = TRUE
  )
})

test_that("check_integer raises an error if any element of x is below min", {
  expect_error(check_integer("x", c(0, 1, 2), type = "vector", min = 1),
    "x must be higher than or equal to 1.",
    fixed = TRUE
  )
})
test_that("check_integer raises an error if any element of x is above max", {
  expect_error(check_integer("x", c(1, 3), type = "vector", max = 2),
    "x must be lower than or equal to 2.",
    fixed = TRUE
  )
  expect_error(
    check_integer("x", c(1, 3), type = "vector", max = 2, message = "error"),
    "error",
    fixed = TRUE
  )
})

test_that("check_integer passes and returns x when x is valid", {
  expect_equal(check_integer(1), 1)
  expect_equal(check_integer(c(1, 2, 3), type = "vector"), c(1, 2, 3))
  expect_equal(check_integer(1.7, float = TRUE), 1.7)
  x <- matrix(1:4, 2, 2)
  rownames(x) <- paste0("R", 1:2)
  colnames(x) <- paste0("C", 1:2)
  expect_equal(check_integer(x, type = "matrix"), x)
  x <- as.data.frame(x)
  expect_equal(check_integer(x, type = "data.frame"), x)
})

# Test check_method
test_that("check_method raises an error when method", {
  expect_error(check_method("toto"))
})
test_that("check_method passes when method is valid", {
  for (method in available_methods()) {
    expect_error(check_method(method), NA)
  }
})

# Test check_nblocks
test_that("check_nblocks raises an error if method is pca and number of blocks
          different from 1", {
  expect_error(check_nblocks(blocks, method = "pca"),
    paste0(
      "3 blocks were provided but the number of",
      " blocks for pca must be 1."
    ),
    fixed = TRUE
  )
})
test_that("check_nblocks raises an error if method is cca and number of blocks
          is different from 2", {
  expect_error(check_nblocks(blocks, method = "cca"),
    paste0(
      "3 blocks were provided but the number of",
      " blocks for cca must be 2."
    ),
    fixed = TRUE
  )
})
test_that("check_nblocks passes and returns blocks when blocks is valid", {
  A <- list(blocks[[1]])
  expect_equal(check_nblocks(A, method = "pca"), A)
  A <- list(blocks[[1]], blocks[[2]])
  expect_equal(check_nblocks(A, method = "cca"), A)
})

# Test check_ncomp
test_that("check_ncomp raises an error if there is a superblock and ncomp
          contains at least two distinct values", {
  expect_error(check_ncomp(c(1, 2, 1), blocks, superblock = TRUE),
    paste0(
      "only one number of components must be ",
      "specified (superblock)."
    ),
    fixed = TRUE
  )
})
test_that("check_ncomp raises an error if blocks and ncomp have different
          lengths and length of ncomp is greater than 1", {
  expect_error(check_ncomp(c(2, 2), blocks),
    paste0(
      "ncomp must have the same size (actually 2) ",
      "as the number of blocks (3)."
    ),
    fixed = TRUE
  )
})
test_that("check_ncomp raises an error if any element of ncomp is greater than
          the number of variables in the corresponding block", {
  expect_error(check_ncomp(c(2, 7, 2), blocks),
    paste0(
      "ncomp[2] must be lower than the number of variables for block 2,",
      " i.e. 2."
    ),
    fixed = TRUE
  )
})
test_that("check_ncomp raises an error if ncomp is greater than the number of
          columns in the superblock", {
  expect_error(check_ncomp(12, blocks, superblock = TRUE),
    paste0(
      "the number of components must be lower ",
      "than the number of variables in the superblock,",
      " i.e. 11."
    ),
    fixed = TRUE
  )
})
test_that("check_ncomp raises an error for invalid ncomp", {
  expect_error(check_ncomp(c(-1, 1, 1), blocks))
  expect_error(check_ncomp(c(1, 1, 1.7), blocks))
  expect_error(check_ncomp(c(NA, 1, 2), blocks))
})
test_that("check_ncomp passes and returns ncomp when ncomp is valid", {
  expect_equal(check_ncomp(2, blocks), c(2, 2, 2))
  expect_equal(check_ncomp(c(2, 2, 2), blocks), c(2, 2, 2))
})

# Test check_sign_comp
test_that("check_sign_comp changes the sign of weight vector if correlation
          with reference is negative", {
  fit_rgcca <- rgcca(blocks, ncomp = 1)
  a <- fit_rgcca$a
  a[[1]] <- -a[[1]]
  expect_identical(fit_rgcca$a, check_sign_comp(fit_rgcca, a))

  fit_rgcca <- rgcca(blocks, ncomp = 2)
  a <- fit_rgcca$a
  a[[1]][, 1] <- -a[[1]][, 1]
  expect_identical(fit_rgcca$a, check_sign_comp(fit_rgcca, a))
})

# Test check_size_blocks
test_that("check_size_blocks raises an error when number of columns of x does
          not match length of blocks", {
  expect_error(check_size_blocks(blocks, "x", diag(2)),
    paste0(
      "x must have the same number of columns",
      " (actually 2) as the number of blocks (3)."
    ),
    fixed = TRUE
  )
})
test_that("check_size_blocks raises an error when number of rows of x does
          not match specified n_row", {
  expect_error(check_size_blocks(blocks, "x", diag(3), n_row = 2),
    paste0("x must have 2 rows."),
    fixed = TRUE
  )
})
test_that("check_size_blocks raises an error when size of x does
          not match length of blocks", {
  expect_error(check_size_blocks(blocks, "x", c(2, 2)),
    paste0(
      "x must have the same size (actually 2) ",
      "as the number of blocks (3)."
    ),
    fixed = TRUE
  )
})
test_that("check_size_blocks passes when x is valid", {
  expect_error(check_size_blocks(blocks, "x", diag(3)), NA)
  expect_error(check_size_blocks(blocks, "x", c(2, 2, 2)), NA)
})

# Test check_penalty
test_that("check_penalty raises an error if blocks and penalty have different
          lengths and length of penalty is greater than 1", {
  expect_error(check_penalty(c(1, 1), blocks),
    paste0(
      "tau must have the same size ",
      "(actually 2) as the number of blocks (3)."
    ),
    fixed = TRUE
  )
})
test_that("check_penalty raises an error if penalty has two rows but
          ncomp is 1", {
  expect_error(check_penalty(matrix(1, 2, 3), blocks, ncomp = 1),
    paste0("tau must have 1 rows."),
    fixed = TRUE
  )
})
test_that("check_penalty raises an error if any element of sparsity is lower
          than the inverse of the square root of the number of variables in the
          corresponding block", {
  min_sparsity <- 1 / sqrt(NCOL(blocks[[3]]))
  expect_error(check_penalty(c(1, 1, 0.2), blocks, method = "sgcca"),
    paste0(
      "too low sparsity. Sparsity parameter equals 0.2. For SGCCA, ",
      "it must be greater than 1/sqrt(number_column)",
      " (i.e., ", round(min_sparsity, 4), " for block 3)."
    ),
    fixed = TRUE
  )
})
test_that("check_penalty raises an error for invalid penalty", {
  expect_error(check_penalty(c(-1, 1, 1), blocks, method = "rgcca"))
  expect_error(check_penalty(c(1, 1, 2), blocks, method = "rgcca"))
  expect_error(check_penalty(c(NA, 1, 1), blocks, method = "rgcca"))
  expect_error(check_penalty("toto", blocks, method = "rgcca"))
  expect_error(check_penalty(c(-1, 1, 1), blocks, method = "sgcca"))
  expect_error(check_penalty(c(1, 1, 2), blocks, method = "sgcca"))
  expect_error(check_penalty(c(NA, 1, 1), blocks, method = "sgcca"))
  expect_error(check_penalty("optimal", blocks, method = "sgcca"))
})
test_that("check_penalty passes and returns penalty when penalty is valid", {
  expect_equal(check_penalty(1, blocks, method = "rgcca"), c(1, 1, 1))
  expect_equal(
    check_penalty(c(0.8, 1, 0.5), blocks, method = "rgcca"), c(0.8, 1, 0.5)
  )
  expect_equal(
    check_penalty("optimal", blocks, method = "rgcca"),
    c("optimal", "optimal", "optimal")
  )
  expect_equal(
    check_penalty(matrix(1, 5, 3), blocks, method = "rgcca"), matrix(1, 5, 3)
  )
  expect_equal(
    check_penalty(rep(1, 4), blocks, method = "rgcca", superblock = TRUE),
    rep(1, 4)
  )
  expect_equal(check_penalty(1, blocks, method = "sgcca"), c(1, 1, 1))
  expect_equal(
    check_penalty(c(0.8, 1, 0.5), blocks, method = "sgcca"), c(0.8, 1, 0.5)
  )
  expect_equal(
    check_penalty(matrix(1, 5, 3), blocks, method = "sgcca"), matrix(1, 5, 3)
  )
})

# Test check_spars
test_that("check_spars raises an error for invalid sparsity", {
  expect_error(check_spars(0.2, blocks[[3]]))
  expect_error(check_spars(2, blocks[[1]]))
  expect_error(check_spars(NA, blocks[[2]]))
})
test_that("check_spars passes and returns sparsity when sparsity is valid", {
  expect_equal(check_spars(1, blocks[[1]], 1), 1)
  expect_equal(check_spars(0.5, blocks[[3]], 0.5), 0.5)
})

# Test check_tau
test_that("check_tau raises an error for invalid tau", {
  expect_error(check_tau(-1))
  expect_error(check_tau(2))
  expect_error(check_tau(NA))
  expect_error(check_tau("toto"))
})
test_that("check_tau passes and returns tau when tau is valid", {
  expect_equal(check_tau(1), 1)
  expect_equal(check_tau(0.3), 0.3)
  expect_equal(check_tau("optimal"), "optimal")
})

# Test check_scheme
test_that("check_scheme raises an error for invalid scheme", {
  expect_error(check_scheme("toto"),
    paste0(
      "scheme must be one of the following schemes: 'horst', ",
      "'centroid', 'factorial' or a function."
    ),
    fixed = TRUE
  )
})
test_that("check_scheme passes when scheme is valid", {
  expect_error(check_scheme("horst"), NA)
  expect_error(check_scheme("centroid"), NA)
  expect_error(check_scheme("factorial"), NA)
  g <- function(x) x^2
  expect_error(check_scheme(g), NA)
})

# Test check_rank
n <- nrow(blocks[[1]])
array_blocks <- blocks
array_blocks[[4]] <- array(seq(n * 10 * 5), dim = c(n, 10, 5))
test_that("check_rank raises an error for invalid rank", {
  expect_error(
    check_rank("toto", array_blocks, rep(1, 4), 1),
    "rank must be numeric.", fixed = TRUE
  )
  expect_error(
    check_rank(rep(1, 5), array_blocks, rep(1, 4), 1),
    paste0("rank must have the same number of columns (actually 5) ",
           "as the number of blocks (4)."),
    fixed = TRUE
  )
  expect_error(
    check_rank(c(1, 1, 1, 7), array_blocks, c(1, 1, 1, 2), 1),
    paste0("rank[4] should be comprise between 1 and 5 (that is the number of ",
           "variables of the mode bearing the orthogonality for block 4)."),
    fixed = TRUE
  )
  expect_error(
    check_rank(matrix(rep(1, 8), nrow = 2), array_blocks, rep(1, 4), 1),
    "rank must have 1 rows.", fixed = TRUE
  )
})
test_that("check_rank passes when rank is valid", {
  expect_equal(check_rank(1, array_blocks, rep(1, 4), 1), matrix(1, 1, 4))
  r <- matrix(1, 2, 4)
  expect_equal(check_rank(r, array_blocks, rep(1, 4), 2), r)
  r[2, 4] <- 7
  expect_equal(check_rank(r, array_blocks, rep(1, 4), 2), r)
})
