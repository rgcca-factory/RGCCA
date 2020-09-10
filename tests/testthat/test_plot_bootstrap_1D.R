data("Russett")
blocks <- list(
    agriculture = Russett[, seq(3)],
    industry = Russett[, 4:5],
    politic = Russett[, 6:11] )
# Rgcca
#---------------
rgcca_out <- rgcca(blocks, type = "rgcca")
boot <- bootstrap(rgcca_out, 100, n_cores = 1)
plot(boot)
p<-plot(boot,block=1)
p<-plot(boot,block=2)
p+scale_fill_manual(values=c("blue","red"))
# Sparsity
#---------------
rgcca_out <- rgcca(blocks, sparsity = 0.71, type = "sgcca")
boot <- bootstrap(rgcca_out, 100, n_cores = 1)
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

a=matrix(rnorm(60),10,6);colnames(a)=paste0("V",1:6);rownames(a)=paste("S",1:10)
b=matrix(rnorm(70),10,7);colnames(b)=paste0("W",1:7);rownames(b)=paste("S",1:10)
A=list(a=a,b=b)
res=rgcca(A,type="sgcca",sparsity=0.7)
b=bootstrap(res,n_boot=100,n_cores=1)
plot_bootstrap_1D(b,x="occurrences")

