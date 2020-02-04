#'# test bootstrap

#'''

 data("Russett")
 blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
     politic = Russett[, 6:11] )
 rgcca_out = rgcca(blocks)
 bootstrap(rgcca_out, n_boot = 2, n_cores = 1)
 bootstrap(rgcca_out, n_boot = 2, n_cores = 1, blocks = lapply(blocks, scale),
  superblock = FALSE)

 
 blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
               politic = Russett[, 6:11] )
 resRGCCA=rgcca(blocks,ncomp=c(2,2,2))

 resBootstrap=bootstrap(rgcca=resRGCCA)
 select_var=get_bootstrap(resRGCCA,resBootstrap)
#  plot_bootstrap_1D(select_var)
#  
# testthat("bootstrap_1",{expect_true(abs(select_var["labo", 1]-mean(c(resBootstrap[[1]][[4]]["labo", 1], resBootstrap[[2]][[4]]["labo", 1])))<2e-16)})
# 
# testthat("bootstrap_2",{expect_true(select_var["labo", 2]==resRGCCA$a[[4]]["labo", 1])})
#  
# 
# blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#               politic = Russett[, 6:11] )
# blocks[[1]][1:3,1]=NA
# blocks[[1]][4,]=NA
# resRGCCA=rgcca(blocks,ncomp=c(2,2,2))
# set.seed(seed=18)
# resBootstrap=bootstrap( rgcca=resRGCCA,n_boot = 2)
# select_var=get_bootstrap(resRGCCA,resBootstrap)
# plot_bootstrap_1D(select_var)
# 
# select_var["labo", 1]==mean(c(resBootstrap[[1]][[4]]["labo", 1], resBootstrap[[2]][[4]]["labo", 1]))
# 
# select_var["labo", 2]==resRGCCA$a[[4]]["labo", 1]
