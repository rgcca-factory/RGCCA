#'# defl_select

#'''
#'

set.seed(42); X1 = matrix(rnorm(15),3,5);
set.seed(22); X2 = matrix(rnorm(12),3,4);
set.seed(2);  X3 = matrix(rnorm(21),3,7);
A = list(X1,X2,X3)
yy = list(
  matrix(c(1, 0, 0), 3, 1),
  matrix(c(0, 0, 1), 3, 1),
  matrix(c(1/sqrt(2), 1/sqrt(2), 0), 3, 1)
) # projection vectors are identity
res = defl.select(yy = yy, rr = A, nncomp = c(1,1,1), nn = 1)

test_that("defl.select returns rr when nn > nncomp and projection matrix is 0",
          {
            res = defl.select(yy = yy, rr = A, nncomp = c(1,1,1), nn = 2)
            P = lapply(A, function(x) rep(0, NCOL(x)))
            expect_identical(res$resdefl, A)
            expect_equal(res$pdefl, P)
          })

test_that("defl.select returns a residual that is orthogonal to yy and that rr
          can be reconstructed from projection matrix and residual", {
            res = defl.select(yy = yy, rr = A, nncomp = c(1,1,1), nn = 1)
            for (q in 1:length(A)) {
              expect_true(all(abs(t(yy[[q]]) %*% res$resdefl[[q]]) < 1e-14))
              expect_equal(A[[q]], res$resdefl[[q]] + yy[[q]] %*% t(res$pdefl[[q]]))
            }
          })
