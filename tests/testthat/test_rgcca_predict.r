library(gliomaData)

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

test_that("rgcca_predict raises an error if an unknown prediction model is
          given", {
  expect_error(rgcca_predict(fit_rgcca, blocks, 3, prediction_model = "toto"),
    "unknown model.",
    fixed = TRUE
  )
})

test_that("rgcca_predict raises an error if blocks_test has no names", {
  expect_error(
    rgcca_predict(fit_rgcca, unname(blocks), 3),
    "Please provide names for blocks_test."
  )
})

test_that("rgcca_predict raises an error if response is an integer but
          associated block names do not match between train and test", {
  blocks_test <- blocks
  names(blocks_test)[[3]] <- "wrong_name"
  expect_error(rgcca_predict(fit_rgcca, blocks_test, 3),
    paste0(
      "Block to predict was provided as an integer but ",
      "associated block names do not match between train and ",
      "test blocks. Reorder your blocks or use a name."
    ),
    fixed = TRUE
  )
})

test_that("rgcca_predict raises an error if response block is not present in
          training blocks", {
  blocks_test <- blocks
  names(blocks_test)[[3]] <- "response"
  expect_error(rgcca_predict(fit_rgcca, blocks_test, response = "response"),
    paste0(
      "The block to predict is not among both train and ",
      "test blocks. Please provide an appropriate one."
    ),
    fixed = TRUE
  )
})

test_that("rgcca_predict raises an error if response block dimensions do not
          match", {
  blocks_test <- blocks
  blocks_test[[3]] <- blocks_test[[3]][, 1]
  expect_error(rgcca_predict(fit_rgcca, blocks_test, 3),
    "Dimensions of response do not match",
    fixed = TRUE
  )
})

test_that("rgcca_predict raises an error if the projected blocks are constant
          within classes", {
  data(ge_cgh_locIGR)
  blocks <- ge_cgh_locIGR$multiblocks
  Loc <- factor(ge_cgh_locIGR$y)
  levels(Loc) <- colnames(ge_cgh_locIGR$multiblocks$y)
  blocks[[3]] <- Loc
  fit_rgcca <- rgcca(blocks, tau = 0, response = 3)
  expect_error(
    rgcca_predict(fit_rgcca, blocks, response = 3, prediction_model = "lda"),
    "overfitting model.",
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
    rgcca_res = fit_rgcca, blocks_test = A,
    response = response
  )
  res_lm <- apply(fit_rgcca$call$raw[[response]], 2, function(x) {
    lm(x ~ fit_rgcca$Y[[1]][, 1:3] + fit_rgcca$Y[[2]][, 1:2])$residuals
  })
  score_lm <- mean(apply(res_lm, 2, function(x) {
    return(sqrt(mean(x^2, na.rm = T)))
  }))
  expect_equal(res_predict$score, score_lm)
  expect_equal(as.matrix(A[[response]] - res_predict$prediction), res_lm)
})

# Classification
#---------------
test_that("rgcca_predict with lda predictor gives the same prediction as
          applying lda directly on rgcca score Y", {
  A <- lapply(blocks_classif, function(x) x[1:32, ])
  response <- 3
  fit_rgcca <- rgcca(A, tau = 1, ncomp = c(3, 2, 1), response = response)
  res_predict <- rgcca_predict(fit_rgcca,
    blocks_test = A,
    prediction_model = "lda", response = "politic"
  )
  Y <- data.frame(cbind(fit_rgcca$Y[[1]][, 1:3], fit_rgcca$Y[[2]][, 1:2]))
  res_lda <- MASS::lda(fit_rgcca$call$raw[[response]] ~ as.matrix(Y))
  prediction_lda <- predict(res_lda, Y)$class
  expect_equal(drop(res_predict$prediction), as.numeric(prediction_lda))
})
