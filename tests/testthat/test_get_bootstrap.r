data("Russett")
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)

rgcca_out <- rgcca(blocks)
boot <- rgcca_bootstrap(rgcca_out, n_boot = 5, n_cores = 1, verbose = FALSE)

test_that("get_bootstrap returns the expected statistics", {
  df <- get_bootstrap(boot,
    type = "weights", block = 2, comp = 1,
    empirical = TRUE
  )
  df_test <- data.frame(
    estimate = rgcca_out$a[[2]][, 1],
    mean = apply(boot$bootstrap$W[[1]][[2]], 1, mean, na.rm = TRUE),
    sd = apply(boot$bootstrap$W[[1]][[2]], 1, sd, na.rm = TRUE),
    lower_bound = apply(boot$bootstrap$W[[1]][[2]], 1, quantile, 0.025),
    upper_bound = apply(boot$bootstrap$W[[1]][[2]], 1, quantile, 0.975)
  )
  expect_equal(
    df_test, df[, c("estimate", "mean", "sd", "lower_bound", "upper_bound")]
  )

  df <- get_bootstrap(boot,
    type = "loadings", block = 3, comp = 1,
    empirical = TRUE
  )
  df_test <- data.frame(
    estimate = cor(rgcca_out$blocks[[3]], rgcca_out$Y[[3]][, 1]),
    mean = apply(boot$bootstrap$L[[1]][[3]], 1, mean, na.rm = TRUE),
    lower_bound = apply(boot$bootstrap$L[[1]][[3]], 1, quantile, 0.025),
    upper_bound = apply(boot$bootstrap$L[[1]][[3]], 1, quantile, 0.975)
  )
  expect_equal(
    df_test, df[, c("estimate", "mean", "lower_bound", "upper_bound")]
  )

  tail <- qnorm(1 - .05 / 2)
  df <- get_bootstrap(boot,
    type = "weights", block = 1, comp = 1,
    empirical = FALSE
  )
  std <- apply(boot$bootstrap$W[[1]][[1]], 1, sd, na.rm = TRUE)
  df_test <- data.frame(
    lower_bound = rgcca_out$a[[1]][, 1] - std * tail,
    upper_bound = rgcca_out$a[[1]][, 1] + std * tail
  )
  expect_equal(df_test, df[, c("lower_bound", "upper_bound")])

  df <- get_bootstrap(boot,
    type = "loadings", block = 1, comp = 1,
    empirical = FALSE
  )
  r <- 0.5 * log(
    (1 + boot$bootstrap$L[[1]][[1]]) / (1 - boot$bootstrap$L[[1]][[1]])
  )
  std <- apply(r, 1, function(x) sd(x, na.rm = TRUE))
  estimate <- cor(rgcca_out$blocks[[1]], rgcca_out$Y[[1]][, 1])
  df_test <- data.frame(
    sd = std,
    lower_bound = estimate - std * tail,
    upper_bound = estimate + std * tail
  )
  expect_equal(df_test, df[, c("sd", "lower_bound", "upper_bound")])
})
