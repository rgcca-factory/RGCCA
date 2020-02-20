
 data("Russett")
 blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
 politic = Russett[, 6:11] )
 rgcca_out = rgcca(blocks, sparsity = 0.75, type = "sgcca")
 boot = bootstrap(rgcca_out, 2, n_cores = 1)
 selected.var = get_bootstrap(boot, n_cores = 1)
 plot_bootstrap_2D(selected.var)
 rgcca_out = rgcca(blocks)
 boot = bootstrap(rgcca_out, 2, n_cores = 1)
 selected.var = get_bootstrap(boot, n_cores = 1)
 plot_bootstrap_2D(selected.var)

