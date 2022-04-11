set.seed(0)
eps <- 0.001

### Test proj_l1_l2
# Case ||a||_1 / ||a||_2 <= s
first_test <- function(a) {
  test_that(paste0(
    "lambda is 0 and norm is saturated for a = (",
    paste0(a, collapse = ", "), "), ||a||_1 / ||a||_2 <= s "
  ), {
    s <- sum(abs(a)) / norm(a, type = "2") + eps
    res <- proj_l1_l2(a, s)
    expect_equal(res$lambda, 0)
    expect_true(res$l2_sat)
  })
}

first_test(c(1, 2, 3))
first_test(c(1, 2, 2))
first_test(c(1, 1, 1, 1))
first_test(rnorm(15))

# Case ||a||_1 / ||a||_2 > s and 1 <= s < sqrt(n_max)
second_test <- function(a) {
  q <- sum(abs(a)) / norm(a, type = "2")
  s <- min(q, sqrt(sum(abs(a) == max(abs(a))))) - eps
  if (s < 1) {
    return()
  }
  test_that(paste0(
    "norm is not saturated for a = (",
    paste0(a, collapse = ", "), "), ||a||_1 / ||a||_2 > s ",
    "and 1 <= s < sqrt(n_max)"
  ), {
    res <- proj_l1_l2(a, s)
    expect_false(res$l2_sat)
    expect_equal(sum(abs(res$sol)), s)
  })
}

second_test(c(1, 2, 3))
second_test(c(1, 2, 2))
second_test(c(1, 1, 1, 1))

# Case ||a||_1 / ||a||_2 > s and sqrt(n_max) <= s <= sqrt(p)
third_test <- function(a) {
  n_max <- sum(abs(a) == max(abs(a)))
  s <- min(
    sum(abs(a)) / norm(a, type = "2") - eps,
    sqrt(n_max)
  )
  if (s < sqrt(n_max)) {
    return()
  }

  test_that(paste0(
    "norm is saturated for a = (",
    paste0(a, collapse = ", "), "), ||a||_1 / ||a||_2 > s ",
    "and sqrt(n_max) <= s <= sqrt(p)"
  ), {
    res <- proj_l1_l2(a, s)
    expect_true(res$l2_sat)
    x <- soft(a, res$lambda)
    x <- (x / norm(x, type = "2"))
    expect_equal(sum(abs(x)), s)
  })
}

third_test(c(1, 2, 3))
third_test(c(1, 2, 2))
third_test(c(1, 1, 1, 1))
third_test(rnorm(15))
