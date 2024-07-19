set.seed(1)
# Building the blocks
data("Russett")
blocks <- list(
  agriculture = Russett[, 1:3],
  industry = Russett[, 4:5],
  politic = Russett[, 6:8]
)

blocks_classif <- list(
  agriculture = Russett[, 1:3],
  industry = Russett[, 4:5],
  politic = matrix(Russett[, 11], ncol = 1)
)
blocks_classif[["politic"]][blocks_classif[["politic"]][, 1] == 1, ] <- "demo"
blocks_classif[["politic"]][blocks_classif[["politic"]][, 1] == 0, ] <- "ndemo"

fit_rgcca <- rgcca(blocks, response = 3)

test_that("rgcca_predict raises an error if rgcca_res has no response", {
  expect_error(rgcca_predict(rgcca(blocks), blocks),
    "RGCCA must use a response block.",
    fixed = TRUE
  )
})

test_that("rgcca_predict raises an error if an unknown prediction model is
          given", {
  expect_error(rgcca_predict(fit_rgcca, blocks, prediction_model = "toto"),
    "unknown model.",
    fixed = TRUE
  )
})

test_that("rgcca_predict raises an error if blocks_test has no names", {
  expect_error(
    rgcca_predict(fit_rgcca, unname(blocks)),
    "Please provide names for blocks_test."
  )
})

test_that("rgcca_predict raises an error if response block dimensions do not
          match", {
  blocks_test <- blocks
  blocks_test[[3]] <- blocks_test[[3]][, 1]
  expect_error(rgcca_predict(fit_rgcca, blocks_test),
    "Dimensions of response do not match",
    fixed = TRUE
  )
})

test_that("rgcca_predict raises a warning if the projected blocks are constant
          within classes", {
  skip_if_not_installed("gliomaData")
  skip_if_not_installed("klaR")
  data("ge_cgh_locIGR", package = "gliomaData")
  blocks <- ge_cgh_locIGR$multiblocks
  Loc <- factor(ge_cgh_locIGR$y)
  levels(Loc) <- colnames(ge_cgh_locIGR$multiblocks$y)
  blocks[[3]] <- Loc
  fit_rgcca <- rgcca(blocks, tau = 0, response = 3)
  expect_warning(
    rgcca_predict(fit_rgcca, blocks, prediction_model = "rda"),
    "overfitting risk.",
    fixed = TRUE
  )
})

# Regression
#-----------
test_that("rgcca_predict with lm predictor gives the same prediction as
          applying lm directly on rgcca score Y", {
  A <- lapply(blocks, function(x) x[1:32, ])
  response <- 3
  fit_rgcca <- rgcca(A,
    tau = c(0.7, 0.8, 0.7), ncomp = c(3, 2, 3),
    response = response
  )
  res_predict <- rgcca_predict(
    rgcca_res = fit_rgcca, blocks_test = A
  )
  res_lm <- apply(fit_rgcca$call$blocks[[response]], 2, function(x) {
    lm(x ~ fit_rgcca$Y[[1]][, 1:3] + fit_rgcca$Y[[2]][, 1:2])$residuals
  })
  score_lm <- mean(apply(res_lm, 2, function(x) {
    return(sqrt(mean(x^2, na.rm = TRUE)))
  }))
  names(score_lm) <- "RMSE"
  expect_equal(res_predict$score, score_lm)
  expect_equal(as.matrix(A[[response]] - res_predict$prediction$test), res_lm)
})

test_that("rgcca_predict returns an empty probs in regression", {
  res_predict <- rgcca_predict(rgcca_res = fit_rgcca)
  expect_equal(nrow(res_predict$probs$train), 0)
  expect_equal(nrow(res_predict$probs$test), 0)
})

# Classification
#---------------
test_that("rgcca_predict with lda predictor gives the same prediction as
          applying lda directly on rgcca score Y", {
  A <- lapply(blocks_classif, function(x) x[1:32, ])
  response <- 3
  fit_rgcca <- rgcca(A, tau = 1, ncomp = c(3, 2, 1), response = response)
  res_predict <- rgcca_predict(fit_rgcca,
    blocks_test = A[-3],
    prediction_model = "lda"
  )
  Y <- data.frame(cbind(fit_rgcca$Y[[1]][, 1:3], fit_rgcca$Y[[2]][, 1:2]))
  res_lda <- MASS::lda(fit_rgcca$call$blocks[[response]] ~ as.matrix(Y))
  prediction_lda <- predict(res_lda, Y)$class
  expect_equal(
    res_predict$prediction$test,
    data.frame(politic = prediction_lda)
  )
})

test_that("rgcca_predict returns probs in classification with adequate model", {
  A <- lapply(blocks_classif, function(x) x[1:32, ])
  B <- lapply(blocks_classif, function(x) x[33:47, ])
  response <- 3
  fit_rgcca <- rgcca(A, tau = 1, ncomp = c(3, 2, 1), response = response)
  res_predict <- rgcca_predict(fit_rgcca,
                               blocks_test = B[-3],
                               prediction_model = "lda"
  )
  expect_equal(nrow(res_predict$probs$train), 32)
  expect_equal(nrow(res_predict$probs$test), 15)
})
