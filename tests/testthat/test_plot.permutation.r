#'# test bootstrap
#-----------------------------
data(Russett)
blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
              politic = Russett[, 6:11] )
res_permut=rgcca_permutation(blocks=blocks,p_spars=TRUE,nperm=5,n_cores=1)
plot(res_permut)


