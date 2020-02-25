data(Russett)
blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
              politic = Russett[, 6:11] )
resRGCCA=rgcca(blocks,ncomp=c(2,2,2),scheme=function(x) x^4, type="sgcca",sparsity = c(.6, .75, .5, 1))
resBootstrap=bootstrap(resRGCCA,n_boot = 5,n_cores=1)
#selected.var=get_bootstrap(resBootstrap)
plot(resBootstrap,i_block=4,n_cores=1)
plot(resBootstrap,type="2D",n_cores=1)
