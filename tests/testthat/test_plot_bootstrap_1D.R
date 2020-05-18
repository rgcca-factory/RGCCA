data("Russett")
blocks <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11] )
rgcca_out <- rgcca(blocks, sparsity = 0.71, type = "sgcca")
boot <- bootstrap(rgcca_out, 10, n_cores = 1)
i_block=1
selected.var <- get_bootstrap(boot, n_cores = 1,display_order=TRUE,i_block=i_block)
n_boot=attributes(selected.var)$n_boot
nvar=length(boot$bootstrap[[1]][[i_block]][,1])
avg_p_occ=mean(selected.var$occurrences)
probComp=avg_p_occ/nvar
p1 <- plot(boot)
q1=qbinom(size=n_boot,prob=probComp,p=0.05,lower.tail = FALSE)
q2=qbinom(size=n_boot,prob=probComp,p=0.01,lower.tail = FALSE)
q3=qbinom(size=n_boot,prob=probComp,p=0.05/nvar,lower.tail = FALSE)

p1<-plot_bootstrap_1D(boot,x="occurrences",y="estimate")

