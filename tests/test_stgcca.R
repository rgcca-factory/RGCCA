
library(RGCCA)
#install.packages('sNPLS')
library(sNPLS)

data(bread)
print(bread)
Xbread <-bread$Xbread
Ybread <- bread$Ybread





bread=list(X=bread$Xbread, response=bread$Ybread)
print('data size')
print(dim(bread$X))
print('response size')
print(length(bread$response))
fit <- rgcca(blocks = bread, tau=0,response = 2,method = "stgcca", ncomp = 2,separable=TRUE,rank=1,
             sparsity=list(c(1,0.5),1),sparse_lambda=0.5,scale=1,
             verbose = TRUE)
summary(fit)






#different color per block if there is more than one
plot(fit, type = "weights",
     display_order = TRUE, cex = 0.7)

#different color per mode if we are focused in one
plot(fit, type = "weights",block=1,
     display_order = TRUE, cex = 0.7)

# CV for lambda sparsity
cv_out <- rgcca_cv(blocks = bread, tau=1,response = 2,method = "stgcca", ncomp = 2,
                   par_value=matrix(c(1,0.8,1,1),nrow=2,ncol=2),scale=TRUE,
                   verbose=TRUE,par_type='sparse_lambda',k=3,rank=2, 
                   
                   upsample=FALSE,metric='RMSE',
                   n_iter_max=2000,n_run=10,n_cores = 1)#prediction_model = "knn",tuning=data.frame(k = c(3,5)),
summary(cv_out)

#CV for sparsity
cv_out <- rgcca_cv(blocks = bread, tau=1,response = 2,method = "stgcca", ncomp = 2,
                   par_value=matrix(list(c(0.8,0.8),c(0.3,0.3),1,1),nrow=2,ncol=2),scale=TRUE,
                   verbose=TRUE,par_type='sparsity',k=3,rank=1, 
                   
                   upsample=FALSE,metric='RMSE',
                   n_iter_max=2000,n_run=10,n_cores = 1)#prediction_model = "knn",tuning=data.frame(k = c(3,5)),
summary(cv_out)




#bootstrap without fit stability

boot_out <- rgcca_bootstrap(fit, n_boot = 500)
summary(boot_out)



fit_stab <- rgcca_stability(fit, method='tgcca',
 ,n_boot = 100, verbose = TRUE, n_cores = 2)

boot_out <- rgcca_bootstrap(fit_stab, n_boot = 500)
summary(boot_out)

#if a block is selected, colors per mode
plot(boot_out,type='weights',block=1,
display_order = FALSE,
n_mark = 50, cex = 0.7, cex_sub = 6,
show_star = TRUE)


#otherwise just by blocks
plot(boot_out,type='weights',
     display_order = FALSE,
     n_mark = 50, cex = 0.7, cex_sub = 6,
     show_star = TRUE)


