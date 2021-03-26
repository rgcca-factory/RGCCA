data("Russett")
blocks <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11] )

# Rgcca
#---------------
rgcca_out <- rgcca(blocks, method = "rgcca")
boot <- bootstrap(rgcca_out, 100, n_cores = 1)
plot(boot)
p<-plot(boot,block=1)
p<-plot(boot,block=2)
p+scale_fill_manual(values=c("blue","red"))

# Sparsity
#---------------
rgcca_out <- rgcca(blocks, sparsity = 0.71, method = "sgcca")
rgcca_stab <- rgcca_stability(rgcca_out, n_cores = 1)
boot <- bootstrap(rgcca_stab, 100, n_cores = 1)
i_block=1
selected.var <- get_bootstrap(boot,display_order=TRUE,block=i_block)
n_boot=attributes(selected.var)$n_boot
p1<-plot_bootstrap_1D(boot, type = "weight")
