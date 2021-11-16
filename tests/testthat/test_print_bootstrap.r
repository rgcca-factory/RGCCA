#'# print.rgcca
data(Russett)
X_agric = as.matrix(Russett[, c("gini", "farm", "rent")])
X_ind = as.matrix(Russett[ , c("gnpr", "labo")])
X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
A = list(agri=X_agric, ind=X_ind, pol=X_polit)
C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
res = rgcca(A, ncomp = rep(1,3), tau = c(1, 1, 1),
            scheme = "factorial", scale = TRUE,
            verbose = FALSE, superblock = FALSE)
boot = bootstrap(res, n_cores=1)
print(boot, block=1)
plot(boot, block=1)

res = rgcca(A, ncomp = rep(1, 3), tau = c(1, 1, 1),
            scheme = "factorial", scale = TRUE,
            verbose = FALSE, superblock = TRUE)
boot = bootstrap(res, n_cores = 1)
print(boot)
plot(boot)

res = rgcca(A, ncomp = c(2, 1, 2), tau = c(1, 1, 1),
            scheme = "factorial", scale = TRUE,
            verbose = FALSE, superblock = TRUE)
boot = bootstrap(res, n_cores = 1)
print(boot)
plot(boot)


