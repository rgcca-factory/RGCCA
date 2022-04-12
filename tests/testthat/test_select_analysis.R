data(Russett)
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)
J <- length(blocks)

available_methods <- c(
  "rgcca", "sgcca", "pca", "spca", "pls", "spls", "cca",
  "ifa", "ra", "gcca", "maxvar", "maxvar-b", "maxvar-a",
  "mcoa", "cpca-1", "cpca-2", "cpca-4", "hpca", "maxbet-b",
  "maxbet", "maxdiff-b", "maxdiff", "sabscor",
  "ssqcor", "ssqcov-1", "ssqcov-2", "ssqcov",
  "sumcor", "sumcov-1", "sumcov-2", "sumcov", "sabscov-1",
  "sabscov-2"
)
one_block_methods <- c("pca", "spca")
two_block_methods <- c("cca", "ra", "ifa", "pls", "spls")
superblock_methods <- c(
  "pca", "spca", "gcca", "maxvar", "maxvar-b", "maxvar-a",
  "mcoa", "cpca-1", "cpca-2", "cpca-4", "hpca"
)
cov_methods <- c(
  "pca", "spca", "pls", "ifa", "maxvar-a", "mcoa", "cpca-1",
  "cpca-2", "cpca-4", "hpca", "maxbet-b", "maxbet", "maxdiff-b",
  "maxdiff", "ssqcov-1", "ssqcov-2", "ssqcov", "sumcov-1",
  "sumcov-2", "sumcov", "sabscov-1", "sabscov-2"
)
cor_methods <- c(
  "cca", "gcca", "maxvar", "maxvar-b", "sabscor", "ssqcor",
  "sumcor"
)
horst_methods <- c(
  "pca", "spca", "pls", "spls", "cca", "ifa", "ra", "cpca-1",
  "maxbet", "maxdiff", "sumcor", "sumcov-1", "sumcov-2",
  "sumcov"
)
factorial_methods <- c(
  "gcca", "maxvar", "maxvar-b", "maxvar-a", "mcoa",
  "cpca-2", "maxbet-b", "maxdiff-b", "ssqcor", "ssqcor",
  "ssqcov-1", "ssqcov-2", "ssqcov"
)
centroid_methods <- c("sabscor", "sabscov-1", "sabscov-2")
x4_methods <- c("cpca-4", "hpca")

run_selection <- function(method, quiet = TRUE, ...) {
  if (method %in% one_block_methods) {
    res <- select_analysis(list(agriculture = blocks[[1]]),
      method = method,
      quiet = quiet, ...
    )
    J <- 1
  } else if (method %in% two_block_methods) {
    res <- select_analysis(
      list(agriculture = blocks[[1]], industry = blocks[[2]]),
      method = method,
      quiet = quiet, ...
    )
    J <- 2
  } else {
    res <- select_analysis(blocks, method = method, quiet = quiet, ...)
    J <- length(blocks)
  }
  return(list(J = J, res = res))
}

### Utility functions to create the connection matrices
name_c <- function(C, names_blocks) {
  rownames(C) <- colnames(C) <- names_blocks
  return(C)
}
c_all <- function(J, blocks) {
  return(name_c(matrix(1, J, J), names(blocks)))
}
c_response <- function(J, blocks, resp = J) {
  names_blocks <- names(blocks)
  if (J > length(blocks)) names_blocks <- c(names_blocks, "superblock")
  x <- matrix(0, J, J)
  x[, resp] <- x[resp, ] <- 1
  x[resp, resp] <- 0
  return(name_c(x, names_blocks))
}
c_pair <- function(J, blocks) {
  return(name_c(1 - diag(J), names(blocks)))
}

test_that("superblock methods sets all attributes of a superblock", {
  for (method in available_methods) {
    tmp <- run_selection(method, superblock = TRUE)
    res <- tmp$res
    J <- tmp$J

    if (method %in% c(superblock_methods, "rgcca", "sgcca")) {
      expect_true(res$superblock)
      expect_equal(res$connection, c_response(J + 1, blocks[seq(J)]))
      expect_equal(length(res$ncomp), J + 1)
      expect_equal(length(unique(res$ncomp)), 1)
      expect_equal(length(res$penalty), J + 1)
    } else {
      expect_false(res$superblock)
      expect_equal(dim(res$connection), c(J, J))
      expect_equal(length(res$penalty), J)
    }
  }

  method <- "rgcca"
  tau <- matrix(runif(6), 2, 3)
  tmp <- run_selection(method, superblock = TRUE, tau = tau, ncomp = 2)
  res <- tmp$res
  expect_equal(res$penalty, cbind(tau, 1))
  tau <- matrix(runif(8), 2, 4)
  tmp <- run_selection(method, superblock = TRUE, tau = tau, ncomp = 2)
  res <- tmp$res
  expect_equal(res$penalty, tau)

  method <- "spca"
  res <- run_selection(method, superblock = TRUE, sparsity = 0.7)$res
  expect_equal(res$method, "spca")
  expect_equal(res$penalty, c(0.7, 0.7))
  expect_equal(res$par, "sparsity")
})

test_that("cov methods set penalty to 1", {
  for (method in cov_methods) {
    tmp <- run_selection(method)
    res <- tmp$res
    J <- tmp$J

    expect_equal(res$penalty[1:J], rep(1, J))
  }
})

test_that("cor methods set penalty to 0", {
  for (method in cor_methods) {
    tmp <- run_selection(method)
    res <- tmp$res
    J <- tmp$J

    expect_equal(res$penalty[1:J], rep(0, J))
  }
})

test_that("horst methods set scheme to 'horst'", {
  for (method in horst_methods) {
    res <- run_selection(method)$res
    expect_equal(res$scheme, "horst")
  }
})

test_that("factorial methods set scheme to 'factorial'", {
  for (method in factorial_methods) {
    res <- run_selection(method)$res
    expect_equal(res$scheme, "factorial")
  }
})

test_that("centroid methods set scheme to 'centroid'", {
  for (method in centroid_methods) {
    res <- run_selection(method)$res
    expect_equal(res$scheme, "centroid")
  }
})

test_that("x4 methods set scheme to x^4", {
  for (method in x4_methods) {
    res <- run_selection(method)$res
    expect_equal(mode(res$scheme), "function")
    vec <- rnorm(10)
    expect_equal(res$scheme(vec), vec^4)
  }
})

test_that("warnings are produced if quiet is FALSE and params have been
          modified", {
  method <- "gcca"
  expect_warning(run_selection(method, quiet = TRUE), regexp = NA)
  expect_warning(run_selection(method, quiet = FALSE),
    regexp = paste0(
      "Choice of method 'gcca' overwrote ",
      "parameters 'ncomp', 'scheme', 'tau', ",
      "'connection'."
    )
  )
})
