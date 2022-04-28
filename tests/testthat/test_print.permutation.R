#' # print.permutation
#'''
set.seed(0)
data(Russett)
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11]
)

test_that("print.permutation outputs the expected message", {
  local_reproducible_output(width = 80)
  res <- rgcca_permutation(blocks,
    par_type = "tau", par_length = 2,
    n_perms = 5, n_cores = 1, verbose = FALSE
  )
  out <- capture.output(print(res))
  expect_equal(out[seq(9)], capture.output(print_call(res$call)))
  expect_equal(out[10], "")
  expect_equal(out[11], "Tuning parameters (tau) used: ")
  expect_equal(out[12], "  agriculture industry politic")
  expect_equal(out[13], "1           1        1       1")
  expect_equal(out[14], "2           0        0       0")
  expect_equal(out[15], "")
  expect_equal(
    out[16],
    "  Tuning parameters Criterion Permuted criterion sd     zstat  p-value"
  )
  expect_equal(
    out[17],
    "1 1/1/1              0.717     0.150              0.058  9.770  0.000 "
  )
  expect_equal(
    out[18],
    "2 0/0/0              2.422     0.783              0.162 10.121  0.000 "
  )
  expect_equal(out[19], "")
  expect_equal(out[20], paste0(
    "The best combination is: 0, 0, 0 for a z ",
    "score of 10.1 and a p-value of 0."
  ))
  expect_equal(length(out), 20)

  blocks2 <- c(blocks, blocks)
  names(blocks2) <- NULL
  res <- rgcca_permutation(blocks2,
    par_type = "ncomp", par_length = 2,
    n_perms = 2, n_cores = 1, verbose = FALSE
  )
  out <- capture.output(print(res))
  expect_equal(out[14], "Tuning parameters (ncomp) used: ")
  expect_equal(
    out[19],
    "  Tuning parameters      Criterion Permuted criterion sd       zstat   "
  )
  expect_equal(
    out[20],
    "1 Tuning parameter set 1   6.2353    0.3610             0.0365 161.0061"
  )
  expect_equal(
    out[21],
    "2 Tuning parameter set 2   5.8807    0.5198             0.0297 180.6396"
  )
  expect_equal(out[22], "  p-value ")
  expect_equal(out[23], "1   0.0000")
  expect_equal(out[24], "2   0.0000")
  expect_equal(length(out), 26)
})
