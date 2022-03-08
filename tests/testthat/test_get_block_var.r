#' # get_block_var
#'''
rgcca_out <- list(a = rep(NA, 4))
names(rgcca_out$a) <- LETTERS[seq(4)]
g <- get_bloc_var(rgcca_out$a)
bool <- all.equal(g, LETTERS[seq(3)])
test_that("get_block_var", {
  expect_true(bool)
})
