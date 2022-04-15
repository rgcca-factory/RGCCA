set.seed(1)
# Building blocks
data("Russett")
blocks <- list(
  agriculture = Russett[, 1:3],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)

C <- connection <- 1 - diag(3)
A <- lapply(blocks, function(x) x[1:32, ])
A_scaled <- scaling(A)
A_test <- lapply(blocks, function(x) x[39:47, ])
ncomp <- c(3, 2, 4)

fit.rgcca <- rgcca(A,
  connection = C, tau = c(0.7, 0.8, 0.7),
  ncomp = ncomp, superblock = FALSE
)

#-------------------------------------------------------------------------
# Checking error cases
#-------------------------------------------------------------------------
test_that("rgcca_transform raises an error if rgcca_res is not of type rgcca", {
  expect_error(rgcca_transform(42, A))
})

test_that("rgcca_transform raises an error if X has no names", {
  expect_error(rgcca_transform(fit.rgcca, list(42)),
    "Please provide names for blocks_test.",
    fixed = TRUE
  )
})

test_that("rgcca_transform raises an error if block names do not match", {
  expect_error(rgcca_transform(fit.rgcca, list("wrong_name" = 42)),
    paste0(
      "At least one block from blocks_test was not found in the training",
      " blocks. Please check block names."
    ),
    fixed = TRUE
  )
})

test_that("rgcca_transform raises an error if block dimensions do not match", {
  expect_error(rgcca_transform(fit.rgcca, list("agriculture" = 42)),
    paste0(
      "Dimensions of blocks do not match for block ",
      "agriculture"
    ),
    fixed = TRUE
  )
})

#-------------------------------------------------------------------------
# Checking Y matches the projection on training samples
#-------------------------------------------------------------------------
# Without permutation
projection <- rgcca_transform(fit.rgcca, A)

test_that("rgcca_transform retrieves Y when projecting the training samples", {
  for (j in 1:length(projection)) {
    expect_true(max(abs(projection[[j]] - fit.rgcca$Y[[j]])) < 1e-14)
  }
})

test_that("rgcca_transform retrieves Y with different scaling scenarios", {
  fit.rgcca <- rgcca(A,
    tau = c(0.7, 0.8, 0.7), ncomp = ncomp,
    scale = FALSE, scale_block = FALSE
  )
  projection <- rgcca_transform(fit.rgcca, A)
  for (j in 1:length(projection)) {
    expect_true(max(abs(projection[[j]] - fit.rgcca$Y[[j]])) < 1e-14)
  }

  fit.rgcca <- rgcca(A,
    tau = c(0.7, 0.8, 0.7), ncomp = ncomp,
    scale = FALSE, scale_block = "inertia"
  )
  projection <- rgcca_transform(fit.rgcca, A)
  for (j in 1:length(projection)) {
    expect_true(max(abs(projection[[j]] - fit.rgcca$Y[[j]])) < 1e-14)
  }

  fit.rgcca <- rgcca(A,
    tau = c(0.7, 0.8, 0.7), ncomp = ncomp,
    scale = FALSE, scale_block = "lambda1"
  )
  projection <- rgcca_transform(fit.rgcca, A)
  for (j in 1:length(projection)) {
    expect_true(max(abs(projection[[j]] - fit.rgcca$Y[[j]])) < 1e-14)
  }

  fit.rgcca <- rgcca(A,
    tau = c(0.7, 0.8, 0.7), ncomp = ncomp,
    scale = TRUE, scale_block = "lambda1"
  )
  projection <- rgcca_transform(fit.rgcca, A)
  for (j in 1:length(projection)) {
    expect_true(max(abs(projection[[j]] - fit.rgcca$Y[[j]])) < 1e-14)
  }
})

# With permutations
A_perm <- list(
  agriculture = A[[1]][, c(3, 2, 1)],
  industry = A[[2]],
  politic = A[[3]][, c(2:6, 1)]
)
projection <- rgcca_transform(fit.rgcca, A_perm)

test_that("rgcca_transform retrieves Y when projecting the training samples
          with permuted columns", {
  for (j in 1:length(projection)) {
    expect_true(max(abs(projection[[j]] - fit.rgcca$Y[[j]])) < 1e-14)
  }
})

# With less blocks
projection <- rgcca_transform(fit.rgcca, A[-3])

test_that("rgcca_transform retrieves Y when projecting a subset of the
          training blocks", {
  for (j in 1:length(projection)) {
    expect_true(max(abs(projection[[j]] - fit.rgcca$Y[[j]])) < 1e-14)
  }
})

# With missing values
RussettWithNA <- Russett
RussettWithNA[1:2, 1:3] <- NA
RussettWithNA[3, 4:5] <- NA
RussettWithNA[3, 1] <- NA
blocksNA <- list(
  agriculture = RussettWithNA[, seq(3)],
  industry = RussettWithNA[, 4:5],
  politic = RussettWithNA[, 6:11]
)
fit.rgccaNA <- rgcca(blocksNA,
  connection = C, tau = c(0.7, 0.8, 0.7),
  ncomp = c(1, 1, 1), superblock = FALSE
)
projection <- rgcca_transform(fit.rgccaNA, blocksNA)

test_that("rgcca_transform retrieves Y when projecting the training
          samples with missing values", {
  for (j in 1:length(projection)) {
    expect_true(max(abs(projection[[j]] - fit.rgccaNA$Y[[j]])) < 1e-14)
  }
})

#-------------------------------------------------------------------------
# Checking rgcca_transform works on unseen data
#-------------------------------------------------------------------------
projection <- rgcca_transform(fit.rgcca, A_test)

test_that("rgcca_transform creates projection with the right number of
          dimensions", {
  expect_equal(length(projection), length(A_test))
  for (j in 1:length(projection)) {
    expect_true(all(dim(projection[[j]]) == c(nrow(projection[[j]]), ncomp[j])))
  }
})
