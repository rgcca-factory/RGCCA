#'# test bootstrap
#-----------------------------
data(Russett)
blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
              politic = Russett[, 6:11] )
res_permut=rgcca_permutation(blocks=blocks,par_type = "sparsity",n_perms=5,n_cores=1, init = "random")
plot.permutation(res_permut)
plot(res_permut,bars="sd")
