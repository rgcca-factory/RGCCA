set.seed(1)

data("Russett")
blocks <- list(
  agriculture = Russett[, seq(3)],
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

RussettWithNA <- Russett
RussettWithNA[1:2, 1:3] <- NA
RussettWithNA[3, 4:5] <- NA
RussettWithNA[3, 1] <- NA
blocksNA <- list(
  agriculture = RussettWithNA[, seq(3)],
  industry = RussettWithNA[, 4:5],
  politic = RussettWithNA[, 6:8]
)

blocks_null_sd <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)
blocks_null_sd[[3]][, 6] <- c(1, rep(0, 46))

custom_rgcca_cv_k <- function(blocks, response, prediction_model = "lm", ...) {
  # Step 1: extract line 1 of each block
  A_minus_1 <- lapply(blocks, function(x) {
    return(x[-1, ])
  })
  A_1 <- lapply(blocks, function(x) {
    return(x[1, ])
  })

  # Step 2: Fit RGCCA on block A_minus_1
  rgcca_out_1 <- rgcca(A_minus_1, response = response, ...)

  # Step 3: Evaluate model on A_1
  pred_1 <- rgcca_predict(rgcca_out_1, A_1,
    response = response,
    prediction_model = prediction_model
  )
  return(pred_1$score)
}

test_that("rgcca_cv_k gives the same scores as custom_rgcca_cv_k", {
  # Regression
  response <- 1
  rgcca_out <- rgcca(blocks, response = response, ncomp = 1)
  res_cv_k <- rgcca_cv_k(
    rgcca_out$call, inds = 1, prediction_model = "lm",
    par_type = "tau", par_value = rep(1, length(blocks))
  )
  res_custom_cv_k <- custom_rgcca_cv_k(blocks, response)
  expect_equal(res_cv_k$score, res_custom_cv_k)

  rgcca_out <- rgcca(blocks, response = response, ncomp = 2)
  res_cv_k <- rgcca_cv_k(
    rgcca_out$call, inds = 1, prediction_model = "lm",
    par_type = "tau", par_value = rep(1, length(blocks))
  )
  res_custom_cv_k <- custom_rgcca_cv_k(blocks, response, ncomp = 2)
  expect_equal(res_cv_k$score, res_custom_cv_k)

  # Classification
  response <- 3
  rgcca_out <- rgcca(blocks_classif, response = response, ncomp = 1)
  res_cv_k <- rgcca_cv_k(
    rgcca_out$call, inds = 1, prediction_model = "lda",
    par_type = "tau", par_value = rep(1, length(blocks_classif))
  )
  res_custom_cv_k <- custom_rgcca_cv_k(blocks_classif, response, "lda")
  expect_equal(res_cv_k$score, res_custom_cv_k)

  rgcca_out <- rgcca(blocks_classif, response = response, ncomp = 2)
  res_cv_k <- rgcca_cv_k(
    rgcca_out$call, inds = 1, prediction_model = "lda",
    par_type = "tau", par_value = rep(1, length(blocks_classif))
  )
  res_custom_cv_k <- custom_rgcca_cv_k(blocks_classif, response, "lda",
    ncomp = 2
  )
  expect_equal(res_cv_k$score, res_custom_cv_k)

  # With missing values
  response <- 1
  rgcca_out <- rgcca(blocksNA,
    response = response, ncomp = 1,
    NA_method = "nipals"
  )
  res_cv_k <- rgcca_cv_k(
    rgcca_out$call, inds = 1, prediction_model = "lm",
    par_type = "tau", par_value = rep(1, length(blocksNA))
  )
  res_custom_cv_k <- custom_rgcca_cv_k(blocksNA, response, NA_method = "nipals")
  expect_equal(res_cv_k$score, res_custom_cv_k)

  rgcca_out <- rgcca(blocksNA,
    response = response, ncomp = 1,
    NA_method = "complete"
  )
  res_cv_k <- rgcca_cv_k(
    rgcca_out$call, inds = 1, prediction_model = "lm",
    par_type = "tau", par_value = rep(1, length(blocksNA))
  )
  res_custom_cv_k <- custom_rgcca_cv_k(blocksNA, response,
    NA_method = "complete"
  )
  expect_equal(res_cv_k$score, res_custom_cv_k)

  # With a column of null variance in the training set
  response <- 1
  rgcca_out <- rgcca(blocks_null_sd, response = response, ncomp = 1)
  res_cv_k <- rgcca_cv_k(
    rgcca_out$call, inds = 1, prediction_model = "lm",
    par_type = "tau", par_value = rep(1, length(blocks_null_sd))
  )
  res_custom_cv_k <- custom_rgcca_cv_k(blocks_null_sd, response)
  expect_equal(res_cv_k$score, res_custom_cv_k)
})
