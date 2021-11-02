data(Russett)
X_agric = as.matrix(Russett[, c("gini","farm","rent")]);
X_ind = as.matrix(Russett[, c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
blocks = list(X_agric, X_ind, X_polit);
connection = matrix(c(0, 0, 1,
             0, 0, 1,
             1, 1, 0), 3, 3);
blocks = scaling(blocks, scale = F, scale_block = F)
fit_rgccak = rgccak(blocks, connection, tau = c(1, 1, 1), scheme = "factorial")

test_that("criterion is increasing", {
  n_iter    = length(fit_rgccak$crit)
  crit_old  = fit_rgccak$crit[1:(n_iter - 1)]
  crit_next = fit_rgccak$crit[2:n_iter]
  expect_true(all(crit_next - crit_old >= 0))
})
