#' # print.cval
#'''
set.seed(0)
data(Russett)
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)

test_that("print.cval outputs the expected message", {
  res <- rgcca_cv(blocks,
    response = 3, method = "rgcca", par_type = "tau",
    par_value = c(0, 0.2, 0.3), n_run = 1, n_cores = 1,
    par_length = 2, verbose = FALSE
  )
  out <- capture.output(print(res))
  expect_equal(out[seq(9)], capture.output(print_call(res$call)))
  expect_equal(out[10], "")
  expect_equal(out[11], "Tuning parameters (tau) used: ")
  expect_equal(out[12], "  agriculture industry politic")
  expect_equal(out[13], "1           0      0.2     0.3")
  expect_equal(out[14], "2           0      0.0     0.0")
  expect_equal(out[15], "")
  expect_equal(out[16], "Validation: kfold with 5 folds and 1 run(s)) ")
  expect_equal(out[17], "Prediction model: lm ")
  expect_equal(out[18], "")
  expect_equal(out[19], "  Tuning parameters Median error  2.5% 97.5%")
  expect_equal(out[20], "1       0.0/0.2/0.3        0.772 0.709 0.874")
  expect_equal(out[21], "2       0.0/0.0/0.0        0.722 0.661 0.872")
  expect_equal(out[22], "")
  expect_equal(
    out[23], "The best combination is: 0 0 0 for a mean CV error of 0.722 "
  )
  expect_equal(length(out), 23)

  res <- rgcca_cv(blocks,
    validation = "loo",
    response = 3, method = "sgcca", par_type = "sparsity",
    n_run = 1, n_cores = 1, par_length = 2, verbose = FALSE
  )
  out <- capture.output(print(res, bars = "sd"))
  expect_equal(out[11], "Tuning parameters (sparsity) used: ")
  expect_equal(out[12], "  agriculture industry politic")
  expect_equal(out[13], "1       1.000    1.000   1.000")
  expect_equal(out[14], "2       0.577    0.707   0.408")
  expect_equal(out[16], "Validation: loo ")
  expect_equal(out[19], "  Tuning parameters Mean error Mean - Sd Mean + Sd")
  expect_equal(length(out), 23)

  blocks2 <- c(blocks, blocks)
  names(blocks2) <- NULL
  res <- rgcca_cv(blocks2,
    validation = "kfold", k = 2,
    response = 3, method = "rgcca", par_type = "tau",
    n_run = 1, n_cores = 1, par_length = 2, verbose = FALSE
  )
  out <- capture.output(print(res, bars = "stderr"))
  expect_equal(
    out[22],
    "      Tuning parameters Mean error Mean - Std Error Mean + Std Error"
  )
  expect_equal(
    out[23],
    "1 Tuning parameter set       0.692            0.682            0.702"
  )
  expect_equal(
    out[24],
    "2 Tuning parameter set       0.808            0.687            0.929"
  )
  expect_equal(length(out), 26)
})
