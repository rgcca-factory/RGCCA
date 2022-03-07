#' # print.rgcca
#'''
data(Russett)
X_agric <- as.matrix(Russett[, c("gini", "farm", "rent")])
X_ind <- as.matrix(Russett[, c("gnpr", "labo")])
X_polit <- as.matrix(Russett[, c("demostab", "dictator")])
A <- list(X_agric, X_ind, X_polit)
C <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
res <- rgcca(A, ncomp = rep(1, 3), tau = c(1, 1, 1), scheme = "factorial", scale = TRUE, verbose = FALSE, superblock = FALSE)
print(res)
