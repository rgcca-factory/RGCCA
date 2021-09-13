#'# test bootstrap

#'''
library("plotly")
 data("Russett")
 A = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
     politic = Russett[, 6:11] )
 perm <- rgcca_permutation(A,par_type="tau", n_perms = 2, n_cores = 1)
 plot_permut_3D(perm)
 perm <- rgcca_permutation(A, par_type = "sparsity", n_perms = 2, n_cores = 1)
 plot_permut_3D(perm)
c1s <- expand.grid(
     lapply(
         seq(length(A)),
         function(x) seq(1 / sqrt(ncol(A[[x]])), 1, by = 0.1)
     )
 )
 perm <- rgcca_permutation(A, par_type = "sparsity", par_value = c1s, n_perms = 2, n_cores = 1)
 plot_permut_3D(perm)
