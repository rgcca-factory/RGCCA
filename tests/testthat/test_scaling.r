data("Russett")
blocks <- list(
  agriculture = Russett[, seq(3)],
  industry = Russett[, 4:5],
  politic = Russett[, 6:11])
tol = 1e-14

# Test matrix blocks
blocks3=lapply(blocks,scale)
blocks3=lapply(blocks3,function(x) return(x/sqrt(ncol(x))))
blocks2=scaling(blocks,scale=TRUE,scale_block=TRUE,bias=FALSE)
test_that("scaling_default_1", {
  expect_true(sum(abs(blocks3[[1]] - blocks2[[1]])) < tol)
  expect_true(sum(abs(blocks3[[2]] - blocks2[[2]])) < tol)
  expect_true(sum(abs(blocks3[[3]] - blocks2[[3]])) < tol)
})
test_that("reverse_scaling", {
  center          = lapply(blocks2, function(x) -attr(x, "scaled:center") / attr(x, "scaled:scale"))
  scale           = lapply(blocks2, function(x) 1/attr(x, "scaled:scale"))
  unscaled_blocks = apply_scaling(blocks2, center, scale)
  expect_true(max(abs(blocks[[1]] - unscaled_blocks[[1]])) < tol)
  expect_true(max(abs(blocks[[2]] - unscaled_blocks[[2]])) < tol)
  expect_true(max(abs(blocks[[3]] - unscaled_blocks[[3]])) < tol)
})

# Test tensor blocks
set.seed(0)
blocks = helper.generate_blocks(list(
  c(40, 20, 30), c(40, 18, 27, 12)
))
scaled_blocks  = scaling(blocks,scale=TRUE,scale_block=TRUE,bias=FALSE)
scaled_blocks2 = lapply(blocks, function(x) {
  block  = matrix(x, nrow = NROW(x))
  scale  = apply(block, 1, function(x) sqrt(sum(x ^ 2))) * sqrt(NCOL(block))
  block  = apply(block, -1, function(x) x / scale)
  center = apply(block, 2, mean)
  block  = t(apply(block, 1, function(y) y - center))
  array(block, dim = dim(x))
})
test_that("scaling_default_tensor", {
  expect_true(sum(abs(scaled_blocks[[1]] - scaled_blocks2[[1]])) < tol)
  expect_true(sum(abs(scaled_blocks[[2]] - scaled_blocks2[[2]])) < tol)
})
test_that("reverse_scaling_tensor", {
  center          = lapply(scaled_blocks, function(x) -attr(x, "scaled:center"))
  scale           = lapply(scaled_blocks, function(x) 1/attr(x, "scaled:scale"))
  unscaled_blocks = apply_scaling(scaled_blocks, center, scale)
  expect_true(max(abs(blocks[[1]] - unscaled_blocks[[1]])) < tol)
  expect_true(max(abs(blocks[[2]] - unscaled_blocks[[2]])) < tol)
})
