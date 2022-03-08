
res <- color_group(seq(10))
test_that("test_color_group", {
  expect_true(length(res) == 10)
})
