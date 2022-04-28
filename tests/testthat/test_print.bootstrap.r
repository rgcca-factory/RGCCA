#' # print.bootstrap
set.seed(0)
data(Russett)
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5]
)

fit.rgcca <- rgcca(blocks, ncomp = c(1, 2), method = "rgcca", tau = 1)

test_that("print.bootstrap outputs the expected message", {
  local_reproducible_output(width = 80)
  res <- bootstrap(fit.rgcca, n_boot = 5, n_cores = 1, verbose = FALSE)
  out <- capture.output(print(res))
  expect_equal(out[seq(8)], capture.output(print_call(fit.rgcca$call)))
  expect_equal(out[9], "")
  expect_equal(
    out[10],
    "Extracted statistics on the block-weight vectors from 5 bootstrap samples "
  )
  expect_equal(out[11], "Component: 1 ")
  expect_equal(out[12], paste0(
    "       estimate       mean         sd lower_bound ",
    "upper_bound bootstrap_ratio"
  ))
  expect_equal(out[13], paste0(
    "gini  0.6283955  0.6094169 0.07770200   0.5041928  ",
    "0.69705679       8.0872498"
  ))
  expect_equal(out[14], paste0(
    "farm  0.7635057  0.7631563 0.06131671   0.7116233  ",
    "0.85564322      12.4518377"
  ))
  expect_equal(out[15], paste0(
    "rent -0.1489231 -0.1229535 0.17046994  -0.3168363  ",
    "0.06241138      -0.8736034"
  ))
  expect_equal(out[16], paste0(
    "gnpr  0.7397265  0.7489023 0.06103677   0.6888806  ",
    "0.81202239      12.1193585"
  ))
  expect_equal(out[17], paste0(
    "labo -0.6729076 -0.6574286 0.07029069  -0.7248742 ",
    "-0.58362603      -9.5732109"
  ))
  expect_equal(out[18], "          pval adjust.pval")
  expect_equal(out[19], "gini 0.0000000   0.0000000")
  expect_equal(out[20], "farm 0.0000000   0.0000000")
  expect_equal(out[21], "rent 0.6666667   0.6666667")
  expect_equal(out[22], "gnpr 0.0000000   0.0000000")
  expect_equal(out[23], "labo 0.0000000   0.0000000")
  expect_equal(out[24], "Component: 2 ")
  expect_equal(out[25], paste0(
    "      estimate      mean         sd lower_bound ",
    "upper_bound bootstrap_ratio"
  ))
  expect_equal(out[26], paste0(
    "gnpr 0.6729076 0.6574286 0.07029069   0.5836260   ",
    "0.7248742        9.573211"
  ))
  expect_equal(out[27], paste0(
    "labo 0.7397265 0.7489023 0.06103677   0.6888806   ",
    "0.8120224       12.119358"
  ))
  expect_equal(out[28], "     pval adjust.pval")
  expect_equal(out[29], "gnpr    0           0")
  expect_equal(out[30], "labo    0           0")
  expect_equal(length(out), 30)

  res <- bootstrap(fit.rgcca, n_boot = 2, n_cores = 1, verbose = FALSE)
  out <- capture.output(print(res, type = "loadings"))
  expect_equal(out[10], paste0(
    "Extracted statistics on the block-loadings vectors from 2 ",
    "bootstrap samples "
  ))
  expect_equal(length(out), 30)
})
