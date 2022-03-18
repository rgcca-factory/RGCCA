vec <- c("demo", "ndemo", "demo")
test_that("dim_without_level", {
  expect_equal(dim(as_disjonctive(vec)), c(3, 2))
})
test_that("dim_with_level", {
  expect_equal(dim(as_disjonctive(vec, levs = c("demo", "ndemo", "nsp"))), c(3, 3))
})
vec2 <- c("demo")
test_that("dim_with_level_for_one", {
  expect_equal(dim(as_disjonctive(vec2, levs = c("demo", "ndemo"))), c(1, 2))
})
