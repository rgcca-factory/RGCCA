AVE <- list(c(0.6, 0.5), c(0.7, 0.45))
rgcca_out <- list(AVE = list(AVE_X = AVE))
# For the superblock (or the last block)
rgcca_out$call$method <- "rgcca"
print_comp(rgcca_out, 1)
# "Axis 1 (70%)"
# For the first block
print_comp(rgcca_out, 2, 1)
# "Axis 2 (50%)"

setAVE <- function() {
  AVE <- list(c(0.6, 0.5), c(0.7, 0.45))
  rgcca_out <- list(AVE = list(AVE_X = AVE))
  rgcca_out$call$method <- "rgcca"
  return(rgcca_out)
}

test_that(
  "print_comp by default",
  expect_equal(print_comp(setAVE()), "Comp. 1 (70%)")
)


test_that(
  "print_comp for the 2nd component and the 1rst block",
  expect_equal(print_comp(setAVE(), 2, 1), "Comp. 2 (50%)")
)

test_that("print_comp for the outer AVE", {
  rgcca_out <- setAVE()
  rgcca_out$AVE$AVE_outer <- 0.87
  expect_equal(print_comp(rgcca_out, outer = TRUE), "First outer comp. : 87%")
  rgcca_out$AVE$AVE_outer[2] <- 0.9
  expect_equal(print_comp(rgcca_out, outer = TRUE), "First outer comp. : 87% & 90%")
})
