tol <- 1e-6
tol2 <- 1e-5

data(Russett)
X_agric <- as.matrix(Russett[, c("gini", "farm", "rent")])
X_ind <- as.matrix(Russett[, c("gnpr", "labo")])
X_polit <- as.matrix(Russett[, c("demostab", "dictator")])
##### Retrieve scaled PCA with RGCCA #####
X <- Russett[, seq(5)]
fit.pca <- prcomp(X, scale = TRUE)

test_that("RGCCA is equivalent to scaled PCA with duplicated X and default
          parameters", {
  fit.rgcca <- rgcca(list(X, X))
  if (sign(fit.rgcca$a[[1]][1]) != sign(fit.pca$rotation[1])) {
    fit.rgcca$a[[1]] <- -fit.rgcca$a[[1]]
  }
  expect_lte(max(abs(fit.pca$rotation[, 1] - fit.rgcca$a[[1]])), tol)
})

test_that("RGCCA is equivalent to scaled PCA with columns split into blocks
          and connected to a superblock", {
  A <- list(X[, 1], X[, 2], X[, 3], X[, 4], X[, 5])
  fit.rgcca <- rgcca(
    blocks = A, superblock = TRUE, tau = 0, scheme = "factorial"
  )
  if (sign(fit.rgcca$Y[[6]][1]) != sign(fit.pca$x[1])) {
    fit.rgcca$Y[[6]] <- -fit.rgcca$Y[[6]]
  }
  expect_equal(cor(fit.rgcca$Y[[6]], fit.pca$x[, 1])[1], 1,
    tolerance = tol
  )
})

test_that("RGCCA is equivalent to scaled PCA when method = 'pca'", {
  fit.rgcca <- rgcca(list(X), method = "pca", ncomp = 5)
  if (sign(fit.rgcca$a[[1]][1, 1]) != sign(fit.pca$rotation[1, 1])) {
    fit.rgcca$a[[1]][, 1] <- -fit.rgcca$a[[1]][, 1]
  }
  expect_lte(max(abs(fit.pca$rotation[, 1] - fit.rgcca$a[[1]][, 1])), tol)
  if (sign(fit.rgcca$a[[1]][1, 2]) != sign(fit.pca$rotation[1, 2])) {
    fit.rgcca$a[[1]][, 2] <- -fit.rgcca$a[[1]][, 2]
  }
  expect_lte(max(abs(fit.pca$rotation[, 2] - fit.rgcca$a[[1]][, 2])), tol)
  # Check AVE
  AVE_pca <- summary(fit.pca)$importance[2, ]
  AVE_rgcca <- fit.rgcca$AVE$AVE_X[[1]]
  expect_lte(max(abs(AVE_pca - AVE_rgcca)), tol2)
})

##### Retrieve unscaled PCA with RGCCA #####
X <- Russett[, seq(5)]
fit.pca <- prcomp(X, scale = FALSE)

test_that("RGCCA is equivalent to unscaled PCA with duplicated X and default
          parameters", {
  fit.rgcca <- rgcca(list(X, X), scale = FALSE, scale_block = FALSE)
  if (sign(fit.rgcca$a[[1]][1]) != sign(fit.pca$rotation[1])) {
    fit.rgcca$a[[1]] <- -fit.rgcca$a[[1]]
  }
  expect_lte(max(abs(fit.pca$rotation[, 1] - fit.rgcca$a[[1]])), tol)
})

test_that("RGCCA is equivalent to unscaled PCA with columns split into blocks
          and connected to a superblock", {
  A <- list(X[, 1], X[, 2], X[, 3], X[, 4], X[, 5])
  fit.rgcca <- rgcca(
    blocks = A, superblock = TRUE, tau = 0,
    scheme = "factorial", scale = FALSE,
    scale_block = FALSE
  )
  if (sign(fit.rgcca$Y[[6]][1]) != sign(fit.pca$x[1])) {
    fit.rgcca$Y[[6]] <- -fit.rgcca$Y[[6]]
  }
  expect_equal(cor(fit.rgcca$Y[[6]], fit.pca$x[, 1])[1], 1,
    tolerance = 1e-5
  )
})

test_that("RGCCA is equivalent to unscaled PCA when method = 'pca'", {
  fit.rgcca <- rgcca(list(X),
    method = "pca", ncomp = 2, scale = FALSE,
    scale_block = FALSE
  )
  if (sign(fit.rgcca$a[[1]][1, 1]) != sign(fit.pca$rotation[1, 1])) {
    fit.rgcca$a[[1]][, 1] <- -fit.rgcca$a[[1]][, 1]
  }
  expect_lte(max(abs(fit.pca$rotation[, 1] - fit.rgcca$a[[1]][, 1])), tol)
  if (sign(fit.rgcca$a[[1]][1, 2]) != sign(fit.pca$rotation[1, 2])) {
    fit.rgcca$a[[1]][, 2] <- -fit.rgcca$a[[1]][, 2]
  }
  expect_lte(max(abs(fit.pca$rotation[, 2] - fit.rgcca$a[[1]][, 2])), tol)
})

##### Retrieve PLS with RGCCA #####
A <- list(agriculture = X_agric, industry = X_ind)
fit.svd <- svd(t(scale(X_agric)) %*% scale(X_ind))

test_that("RGCCA is equivalent to PLS when there are two blocks and tau = 1", {
  fit.rgcca <- rgcca(
    blocks = A, scale = TRUE, scale_block = FALSE,
    bias = FALSE, scheme = "horst", tol = 1e-16
  )
  expect_lte(max(abs(fit.svd$u[, 1] - fit.rgcca$a[[1]])), tol)
  expect_lte(max(abs(fit.svd$v[, 1] - fit.rgcca$a[[2]])), tol)
})

test_that("RGCCA is equivalent to PLS when method = 'pls'", {
  fit.rgcca <- rgcca(
    blocks = A, method = "pls", scale = TRUE,
    scale_block = FALSE, tol = 1e-16, bias = FALSE
  )
  expect_lte(max(abs(fit.svd$u[, 1] - fit.rgcca$a[[1]][, 1])), tol)
  expect_lte(max(abs(fit.svd$v[, 1] - fit.rgcca$a[[2]][, 1])), tol)
})

##### Retrieve CCA with RGCCA #####
Sigma_11 <- cov(X_agric)
Sigma_12 <- cov(X_agric, X_ind)
Sigma_22 <- cov(X_ind)

sqrt_matrix <- function(M) {
  eig <- eigen(M)
  M_sqrt <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  Minv_sqrt <- eig$vectors %*% diag(eig$values^(-1 / 2)) %*% t(eig$vectors)
  return(list(M_sqrt = M_sqrt, Minv_sqrt = Minv_sqrt))
}

fit.svd <- svd(
  sqrt_matrix(Sigma_11)[[2]] %*% Sigma_12 %*% sqrt_matrix(Sigma_22)[[2]]
)

test_that("RGCCA is equivalent to CCA when there are two blocks and tau = 0", {
  fit.rgcca <- rgcca(
    blocks = A, scale = FALSE, tau = 0, scheme = "horst",
    scale_block = FALSE, tol = 1e-16, bias = FALSE
  )
  expect_lte(max(abs(
    fit.svd$u[, 1] - sqrt_matrix(Sigma_11)[[1]] %*% fit.rgcca$a[[1]]
  )), tol)
  expect_lte(max(abs(
    fit.svd$v[, 1] - sqrt_matrix(Sigma_22)[[1]] %*% fit.rgcca$a[[2]]
  )), tol)
})

test_that("RGCCA is equivalent to CCA when method = 'cca'", {
  fit.rgcca <- rgcca(
    blocks = A, method = "cca", scale = FALSE,
    scale_block = FALSE, tol = 1e-16, bias = FALSE
  )
  expect_lte(max(abs(
    fit.svd$u[, 1] - sqrt_matrix(Sigma_11)[[1]] %*% fit.rgcca$a[[1]]
  )), tol)
  expect_lte(max(abs(
    fit.svd$v[, 1] - sqrt_matrix(Sigma_22)[[1]] %*% fit.rgcca$a[[2]]
  )), tol)
})

##### Perform OLS with RGCCA #####
y <- Russett[, 4]
fit.lm <- lm(y ~ as.matrix(X_agric))
r2_ols <- summary(fit.lm)$r.squared

test_that("The criterion of RGCCA is the R-squared for the OLS problem", {
  fit.rgcca <- rgcca(
    blocks = list(y, X_agric), scheme = "factorial",
    tau = 0, tol = 1e-16, bias = FALSE
  )
  crit_rgcca <- fit.rgcca$crit[length(fit.rgcca$crit)] / 2
  expect_equal(r2_ols, crit_rgcca, tolerance = tol)
})

##### Verify theoretical relations between superblock and block components #####
A <- list(Agric = X_agric, Ind = X_ind, Polit = X_polit)
J <- length(A)
fit <- rgcca(
  blocks = A, tau = 1, scheme = "factorial", scale = FALSE,
  scale_block = FALSE, bias = FALSE, ncomp = 1, superblock = TRUE
)

test_that("Block weights can be retrieved using the superblock component", {
  for (j in seq_len(J)) {
    a <- t(A[[j]]) %*% fit$Y[[J + 1]]
    if (sign(a[1]) != sign(fit$a[[j]][1])) a <- -a
    expect_lte(max(abs(fit$a[[j]] - a / norm(a, type = "2"))), tol)
  }
})

test_that("Block weights can be retrieved using the superblock weights", {
  idx <- c(0, cumsum(vapply(A, ncol, FUN.VALUE = 1L)))
  for (j in seq_len(J)) {
    a <- fit$a[[J + 1]][seq(1 + idx[j], idx[j + 1])]
    if (sign(a[1]) != sign(fit$a[[j]][1])) a <- -a
    expect_lte(max(abs(fit$a[[j]] - a / norm(a, type = "2"))), tol)
  }
})

##### Retrieve MFA with RGCCA #####
if ("FactoMineR" %in% rownames(installed.packages())) {
  df <- Russett[, c(
    "gini", "farm", "rent", "gnpr", "labo",
    "inst", "ecks", "death", "demostab", "dictator"
  )]
  fit.mfa <- FactoMineR::MFA(df,
                             group = c(3, 2, 5), ncp = 2, type = rep("s", 3),
                             graph = FALSE
  )

  X_agric <- Russett[, c("gini", "farm", "rent")]
  X_ind <- Russett[, c("gnpr", "labo")]
  X_polit <- Russett[, c(
    "inst", "ecks", "death",
    "demostab", "dictator"
  )]
  A <- list(Agric = X_agric, Ind = X_ind, Polit = X_polit)

  test_that("RGCCA is equivalent to MFA with right parameters", {
    fit.mcoa <- rgcca(
      blocks = A, tau = 1, scheme = "factorial", scale = TRUE,
      scale_block = "lambda1", bias = TRUE, ncomp = 2,
      superblock = TRUE, tol = 1e-16
    )

    expect_lte(max(abs(fit.mcoa$Y[[4]][, 1] - fit.mfa$ind$coord[, 1])), tol)
    expect_lte(max(abs(fit.mcoa$Y[[4]][, 2] - fit.mfa$ind$coord[, 2])), tol)
  })
}

##### Test AVE #####
X_agric <- Russett[, c("gini", "farm", "rent")]
X_ind <- Russett[, c("gnpr", "labo")]
X_polit <- Russett[, c(
  "inst", "ecks", "death",
  "demostab", "dictator"
)]
A <- list(Agric = X_agric, Ind = X_ind, Polit = X_polit)

test_that("rgcca produces cumulated AVE that are below 1", {
  res <- rgcca(A, ncomp = rep(2, 3))
  expect_true(all(unlist(lapply(res$AVE$AVE_X_cor, sum)) <= 1 + tol))

  res <- rgcca(A, ncomp = rep(3, 3), response = 2)
  expect_true(all(unlist(lapply(res$AVE$AVE_X_cor, sum)) <= 1 + tol))

  res <- rgcca(A, ncomp = rep(6, 4), superblock = TRUE)
  expect_true(all(unlist(lapply(res$AVE$AVE_X_cor, sum)) <= 1 + tol))
})

test_that("rgcca returns equal AVE and corrected AVE if components are
          not correlated", {
  res <- rgcca(A, ncomp = rep(2, 3))
  expect_true(all.equal(res$AVE$AVE_X_cor, res$AVE$AVE_X))
})

test_that("rgcca does not report AVE for qualitative response block", {
  A[[3]] <- as.factor(A[[3]][, 5])
  res <- rgcca(A, ncomp = rep(2, 3), response = 3)
  expect_equal(names(res$AVE$AVE_X), names(A)[-3])
})
