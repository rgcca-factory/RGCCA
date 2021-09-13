vec=c("demo","ndemo","demo")
test_that("dim_without_level",{expect_true(sum(dim(asDisjonctive(vec))==c(3,2))==2)})
test_that("dim_with_level",{expect_true(sum(dim(asDisjonctive(vec,levs=c("demo","ndemo","nsp")))==c(3,3))==2)})
vec2=c("demo")
test_that("dim_with_level_for_one",{expect_true(sum(dim(asDisjonctive(vec2, levs=c("demo","ndemo")))==c(1,2))==2
)})
