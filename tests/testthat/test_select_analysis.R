data(Russett)
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)
J <- length(blocks)

run_selection <- function(method, quiet = TRUE, ...) {
  J <- length(blocks)
  if (method %in% one_block_methods()) J <- 1
  if (method %in% two_block_methods()) J <- 2

  rgcca_args <- list(
    tau = rep(1, J),
    rank = rep(1, J),
    ncomp = rep(1, J),
    quiet = quiet,
    scheme = "centroid",
    method = method,
    response = NULL,
    sparsity = rep(1, J),
    mode_orth = rep(1, J),
    connection = 1 - diag(J),
    superblock = FALSE
  )
  rgcca_args <- modifyList(rgcca_args, list(...), keep.null = TRUE)

  if (method %in% one_block_methods()) {
    res <- select_analysis(rgcca_args, list(agriculture = blocks[[1]]))
  } else if (method %in% two_block_methods()) {
    res <- select_analysis(
      rgcca_args,
      list(agriculture = blocks[[1]], industry = blocks[[2]])
    )
  } else {
    res <- select_analysis(rgcca_args, blocks)
  }
  return(list(J = J, res = res))
}

test_that("superblock methods sets all attributes of a superblock", {
  for (method in available_methods()) {
    tmp <- run_selection(method, superblock = TRUE)
    res <- tmp$res
    J <- tmp$J

    if (method %in% c(superblock_methods(), "rgcca", "sgcca")) {
      expect_true(res$rgcca_args$superblock)
      expect_equal(
        res$rgcca_args$connection,
        connection_matrix(blocks[seq(J)], type = "response", J = J + 1)
      )
      expect_equal(length(res$rgcca_args$ncomp), J + 1)
      expect_equal(length(unique(res$rgcca_args$ncomp)), 1)
      expect_equal(length(res$rgcca_args[[res$opt$param]]), J + 1)
    } else {
      expect_false(res$rgcca_args$superblock)
      expect_equal(dim(res$rgcca_args$connection), c(J, J))
      expect_equal(length(res$rgcca_args[[res$opt$param]]), J)
    }
  }

  method <- "rgcca"
  tau <- matrix(stats::runif(6), 2, 3)
  tmp <- run_selection(method, superblock = TRUE, tau = tau, ncomp = 2)
  res <- tmp$res
  expect_equal(unname(res$rgcca_args[[res$opt$param]]), cbind(tau, 1))
  tau <- matrix(stats::runif(8), 2, 4)
  tmp <- run_selection(method, superblock = TRUE, tau = tau, ncomp = 2)
  res <- tmp$res
  expect_equal(unname(res$rgcca_args[[res$opt$param]]), tau)

  method <- "spca"
  res <- run_selection(method, superblock = TRUE, sparsity = 0.7)$res
  expect_equal(res$rgcca_args$method, "spca")
  expect_equal(res$rgcca_args[[res$opt$param]], c(0.7, 0.7))
  expect_equal(res$opt$param, "sparsity")
})

test_that("cov methods set penalty to 1", {
  for (method in cov_methods()) {
    tmp <- run_selection(method)
    res <- tmp$res
    J <- tmp$J

    expect_equal(res$rgcca_args[[res$opt$param]][1:J], rep(1, J))
  }
})

test_that("cor methods set penalty to 0", {
  for (method in cor_methods()) {
    tmp <- run_selection(method)
    res <- tmp$res
    J <- tmp$J

    expect_equal(res$rgcca_args[[res$opt$param]][1:J], rep(0, J))
  }
})

test_that("horst methods set scheme to 'horst'", {
  for (method in horst_methods()) {
    res <- run_selection(method)$res
    expect_equal(res$rgcca_args$scheme, "horst")
  }
})

test_that("factorial methods set scheme to 'factorial'", {
  for (method in factorial_methods()) {
    res <- run_selection(method)$res
    expect_equal(res$rgcca_args$scheme, "factorial")
  }
})

test_that("centroid methods set scheme to 'centroid'", {
  for (method in centroid_methods()) {
    res <- run_selection(method)$res
    expect_equal(res$rgcca_args$scheme, "centroid")
  }
})

test_that("x4 methods set scheme to x^4", {
  for (method in x4_methods()) {
    res <- run_selection(method)$res
    expect_equal(mode(res$rgcca_args$scheme), "function")
    vec <- rnorm(10)
    expect_equal(res$rgcca_args$scheme(vec), vec^4)
  }
})

test_that("messages are produced if quiet is FALSE and params have been
          modified", {
  method <- "gcca"
  expect_warning(run_selection(method, quiet = TRUE), regexp = NA)
  expect_message(run_selection(method, quiet = FALSE),
    regexp = paste0(
      "Choice of method 'gcca' overwrote ",
      "parameters 'ncomp', 'scheme', 'tau', ",
      "'superblock', 'connection'."
    )
  )
})
