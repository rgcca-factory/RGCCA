# Load data
data(Russett)
blocks = list(agriculture = Russett[, seq(3)],
              industry = Russett[, 4:5],
              politic = Russett[, 6:11])

# Test check_blockx
test_that("check_blockx raises an error if x is greater than the number of
          blocks", {
  expect_error(check_blockx("x", 7, blocks),
               "x should be lower than 3 (that is the number of blocks).",
               fixed = TRUE)
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
               "x should not be NA.",
               fixed = TRUE)
})
test_that("check_boolean raises an error if x is not logical", {
  expect_error(check_boolean("x", 0),
               "x should be TRUE or FALSE.",
               fixed = TRUE)
})
test_that("check_boolean raises an error if type is scalar and x not of
          length 1", {
  expect_error(check_boolean("x", c(TRUE, FALSE), type = "scalar"),
               "x should be of length 1.",
               fixed = TRUE)
})
test_that("check_boolean passes when x is valid", {
  expect_error(check_boolean(TRUE), NA)
  expect_error(check_boolean(c(TRUE, FALSE), type = "vector"), NA)
})

# Test check_colors
test_that("check_colors raises an error if colors contains unknown colors", {
  expect_error(check_colors(c("white", "blueen")),
               "colors must be in colors() or a rgb character.",
               fixed = TRUE)
  expect_error(check_colors(c("white", "#79eff16342")),
               "colors must be in colors() or a rgb character.",
               fixed = TRUE)
})
test_that("check_colors passes when colors is valid", {
  expect_error(check_colors(c("#79eff1", "white", "yellow", NA)), NA)
})

# Test check_compx
test_that("check_compx raises an error if x is greater than the number
          of components", {
  expect_error(check_compx("x", 7, c(3, 3, 3), 1),
               paste0("x equals to 7 and should be less than or equal to 3 ",
               "(the number of block component)."),
               fixed = TRUE)
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
test_that("check_connection raises an error if C is not symmetric", {
  expect_error(check_connection(matrix(1:4, 2, 2), blocks),
               "The design matrix C should be symmetric.",
               fixed = TRUE)
})
test_that("check_connection raises an error if C contains NA values", {
  C = diag(2)
  C[1, 1] = NA
  expect_error(check_connection(C, blocks),
               "The design matrix C should not contain NA values.",
               fixed = TRUE)
})
test_that("check_connection raises an error if C contains values outside
          [0, 1]", {
  expect_error(check_connection(diag(2) * 2, blocks),
               "The design matrix C should contain numbers between 0 and 1.",
               fixed = TRUE)
})
test_that("check_connection raises an error if C is the null matrix", {
  expect_error(check_connection(matrix(0, 2, 2), blocks),
               "The design matrix C should not contain only 0.",
               fixed = TRUE)
})
test_that("check_connection raises an error if the dimensions of C do not
          match the number of blocks", {
  expect_error(check_connection(diag(2), blocks),
               paste0("connection matrix should have the same number of ",
                      "columns (actually 2) than the number of blocks (3)."),
               fixed = TRUE)
})
test_that("check_connection raises an error if the dimnames of C do not match
          block names", {
  C = diag(3)
  rownames(C) <- colnames(C) <- paste0("V", 1:3)
  expect_error(check_connection(C, blocks),
               paste0("The design matrix C should have the rownames and the ",
                      "colnames that match with the names of the blocks."),
               fixed = TRUE)
})
test_that("check_connection passes and returns C when C is valid", {
  C = diag(3)
  rownames(C) <- colnames(C) <- names(blocks)
  expect_equal(check_connection(diag(3), blocks), C)
})

# Test check_file
test_that("check_file raises an error when file does not exist", {
  expect_error(check_file("a_file_that_does_not_exist"),
               "a_file_that_does_not_exist does not exist.",
               fixed = TRUE)
})
test_that("check_file passes when file is valid", {
  expect_error(check_file("./test_checks.r"), NA)
})

# Test check_integer
test_that("check_integer raises an error if x is not numeric", {
  expect_error(check_integer("x", "toto"),
               "x should be numeric.",
               fixed = TRUE)
})
test_that("check_integer raises an error if x contains NA", {
  expect_error(check_integer("x", c(42, NA)),
               "x should not be NA.",
               fixed = TRUE)
})
test_that("check_integer raises an error if type is scalar and x not of
          length 1", {
  expect_error(check_integer("x", c(42, 7), type = "scalar"),
               "x should be of length 1.",
               fixed = TRUE)
})
test_that("check_integer raises an error any element of x is a float but float
          is false", {
  expect_error(check_integer("x", c(1, 1.7, 2), type = "vector", float = FALSE),
               "x should be an integer.",
               fixed = TRUE)
})
test_that("check_integer raises an error if any element of x is below min", {
  expect_error(check_integer("x", c(0, 1, 2), type = "vector", min = 1),
               "x should be higher than or equal to 1.",
               fixed = TRUE)
})
test_that("check_integer raises an error if any element of x is above max", {
  expect_error(check_integer("x", c(1, 3), type = "vector", max = 2),
               "x should be lower than or equal to 2.",
               fixed = TRUE)
  expect_error(check_integer("x", c(1, 3), type = "vector", max = 2, message = "error"),
               "error",
               fixed = TRUE)
})
test_that("check_integer passes and returns x when x is valid", {
  expect_equal(check_integer(1), 1)
  expect_equal(check_integer(c(1, 2, 3), type = "vector"), c(1, 2, 3))
  expect_equal(check_integer(1.7, float = TRUE), 1.7)
  x = matrix(1:4, 2, 2)
  rownames(x) = paste0("R", 1:2)
  colnames(x) = paste0("C", 1:2)
  expect_equal(check_integer(x, type = "matrix"), x)
  x = as.data.frame(x)
  expect_equal(check_integer(x, type = "data.frame"), x)
})

# Test check_method
test_that("check_method raises an error when method", {
  expect_error(check_method("toto"))
})
test_that("check_method passes when method is valid", {
  analysis <- c("rgcca", "sgcca", "pca", "spca", "pls", "spls",
                "cca", "ifa", "ra", "gcca", "maxvar", "maxvar-b",
                "maxvar-a", "mcoa","cpca-1", "cpca-2", "cpca-4",
                "hpca", "maxbet-b", "maxbet", "maxdiff-b", "maxdiff",
                "maxvar-a", "sabscor", "ssqcor", "ssqcor", "ssqcov-1",
                "ssqcov-2", "ssqcov", "sumcor", "sumcov-1", "sumcov-2",
                "sumcov", "sabscov", "sabscov-1", "sabscov-2")
  for (method in analysis) {
    expect_error(check_method(method), NA)
  }
})

# Test check_nblocks
test_that("check_nblocks raises an error if type is pca and number of blocks
          different from 1", {
            expect_error(check_nblocks(blocks, type = "pca"),
                         paste0("3 blocks were provided but the number of",
                                " blocks for pca must be 1."),
                         fixed = TRUE)
          })
test_that("check_nblocks raises an error if type is cca and number of blocks
          is different from 2", {
            expect_error(check_nblocks(blocks, type = "cca"),
                         paste0("3 blocks were provided but the number of",
                                " blocks for cca must be 2."),
                         fixed = TRUE)
          })
test_that("check_nblocks passes and returns blocks when blocks is valid", {
  A = list(blocks[[1]])
  expect_equal(check_nblocks(A, type = "pca"), A)
  A = list(blocks[[1]], blocks[[2]])
  expect_equal(check_nblocks(A, type = "cca"), A)
})

# Test check_ncol (Is this the intended behaviour of check_ncol?)
test_that("check_ncol raises an error x[[i_block]] has less than 2 rows", {
  x = list(matrix(1:10, 1, 10))
  expect_error(check_ncol(x, i_block = 1),
               "This output is available only for more than one variable",
               fixed = TRUE)
})
test_that("check_ncol passes when x[[i_block]] has at least 2 rows", {
  x = list(matrix(1:10, 2, 10))
  expect_error(check_ncol(x, i_block = 1), NA)
})

# Test check_ncomp
test_that("check_ncomp raises an error if blocks and ncomp have different
          lengths and length of ncomp is greater than 1", {
            expect_error(check_ncomp(c(2, 2), blocks),
                         paste0("ncomp should have the same size (actually 2) ",
                                "than the number of blocks (3)."),
                         fixed = TRUE)
          })
test_that("check_ncomp raises an error if any element of ncomp is greater than
          the number of variables in the corresponding block", {
            expect_error(check_ncomp(c(2, 7, 2), blocks),
                         paste0("ncomp[2] should be lower than 2 (that is ",
                                "the number of variables of block 2)."),
                         fixed = TRUE)
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

# Test check_quantitative

# Test check_response

# Test check_sign_comp

# Test check_size_blocks
test_that("check_size_blocks raises an error when number of columns of x does
          not match length of blocks", {
            expect_error(check_size_blocks(blocks, "x", diag(2)),
                         paste0("x should have the same number of columns",
                                " (actually 2) than the number of blocks (3)."),
                         fixed = TRUE)
          })
test_that("check_size_blocks raises an error when size of x does
          not match length of blocks", {
            expect_error(check_size_blocks(blocks, "x", c(2, 2)),
                         paste0("x should have the same size (actually 2) ",
                                "than the number of blocks (3)."),
                         fixed = TRUE)
          })
test_that("check_size_blocks passes when x is valid", {
  expect_error(check_size_blocks(blocks, "x", diag(3)), NA)
  expect_error(check_size_blocks(blocks, "x", c(2, 2, 2)), NA)
})

# Test check_size_file

# Test check_spars
test_that("check_spars raises an error if blocks and tau have different
          lengths and length of tau is greater than 1", {
            expect_error(check_spars(blocks, c(1, 1), type = "sgcca"),
                         paste0("sparsity should have the same size ",
                                "(actually 2) than the number of blocks (3)."),
                         fixed = TRUE)
          })
test_that("check_spars raises an error if any element of tau is lower than
          the inverse of the square root of the number of variables in the
          corresponding block", {
            min_sparsity <- lapply(blocks, function(x) 1 / sqrt(NCOL(x)))
            msg = toString(unlist(
              lapply(min_sparsity, function(x)
                ceiling(x * 100) / 100)
            ))
            expect_error(check_spars(blocks, c(1, 1, 0.2), type = "sgcca"),
                         paste0("Sparsity parameter equals to 0.2. For SGCCA, ",
                                "it must be greater than 1/sqrt(number_column)",
                                " (i.e., ", msg, ")."),
                         fixed = TRUE)
          })
test_that("check_spars raises an error for invalid tau", {
  expect_error(check_spars(blocks, c(-1, 1, 1), type = "sgcca"))
  expect_error(check_spars(blocks, c(1, 1, 2), type = "sgcca"))
  expect_error(check_spars(blocks, c(NA, 1, 1), type = "sgcca"))
})
test_that("check_spars passes and returns tau when tau is valid", {
  expect_equal(check_spars(blocks, 1, type = "sgcca"), c(1, 1, 1))
  expect_equal(check_spars(blocks, c(0.8, 1, 0.5), type = "sgcca"), c(0.8, 1, 0.5))
})

# Test check_superblock

# Test check_tau
test_that("check_tau raises an error if blocks and tau have different
          lengths and length of tau is greater than 1", {
            expect_error(check_tau(c(1, 1), blocks),
                         paste0("tau should have the same size ",
                                "(actually 2) than the number of blocks (3)."),
                         fixed = TRUE)
          })
test_that("check_tau raises an error for invalid tau", {
  expect_error(check_tau(c(-1, 1, 1), blocks))
  expect_error(check_tau(c(1, 1, 2), blocks))
  expect_error(check_tau(c(NA, 1, 1), blocks))
})
test_that("check_tau passes and returns tau when tau is valid", {
  expect_equal(check_tau(1, blocks), c(1, 1, 1))
  expect_equal(check_tau(c(0.8, 1, 0.5), blocks), c(0.8, 1, 0.5))
  expect_equal(check_tau(diag(3), blocks), diag(3))
  expect_equal(check_tau(c(1, 1, 1, 1), blocks, superblock = T), c(1, 1, 1, 1))
})
