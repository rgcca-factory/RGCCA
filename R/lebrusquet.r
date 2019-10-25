leb=function(x_k,missing,z,sameBlockWeight=TRUE,weight=NULL,argmax=TRUE)
{#leb(x_k=A[[j]][,k],missing,z=Z[,j])
    n=length(x_k)
    indices_miss=which(missing)
    n_miss=length(indices_miss); 
    n_obs=n-n_miss
    n_inc=n_miss+1
    missing=1:n%in%indices_miss
    indices_obs=(1:n)[!missing]
    x_obs=x_k[indices_obs]
    M1=matrix(0,n,n_miss+1) # n_mis+1 colonne, n lignes
    M1[indices_obs,1]=x_obs
    M1[indices_miss,2:(n_miss+1)]=diag(n_miss)
    v1=rep(0,n);v1[indices_obs]=1
    M2=M1-matrix(v1,ncol=1)%*%matrix(rep(1,n),nrow=1)%*%M1/n_obs
    resSvd=svd(M2)
    #(t(M2)%*%M2)-t(M3)%*%M3
    M3=diag(resSvd$d)%*%t(resSvd$v)
    #solve(M3)-resSvd$v%*%diag(1/(resSvd$d))
    normyres=sqrt(t((t(solve(M3))%*%t(M2)%*%matrix(z,ncol=1)))%*%(t(solve(M3))%*%t(M2)%*%matrix(z,ncol=1)))
    if(argmax)
    {
        yres=sqrt(n)*(t(solve(M3))%*%t(M2)%*%matrix(z,ncol=1))/as.numeric(normyres)
    }
    if(!argmax)
    {
        yres=-sqrt(n)*(t(solve(M3))%*%t(M2)%*%matrix(z,ncol=1))/as.numeric(normyres)
    }
      if(sameBlockWeight){yres=yres/weight}
    # t(yres)%*%yres
    ures=solve(M3)%*%yres
    xres=M2%*%ures
    return(xres)
}
# 
# leb=function(x_k,missing,z,sameBlockWeight=TRUE,p=NULL)
# {#leb(x_k=A[[j]][,k],missing,z=Z[,j])
#     print("criteria before")
#      print(cov(x_k,z))
#     #x_k=scale(x_k)
#     n=length(x_k)
#     indices_miss=which(missing)
#     n_miss=length(indices_miss); 
#     n_obs=n-n_miss
#     n_inc=n_miss+1
#     missing=1:n%in%indices_miss
#     indices_obs=(1:n)[!missing]
#     x_obs=x_k[indices_obs]
#     M1=matrix(0,n,n_miss+1) # n_mis+1 colonne, n lignes
#     M1[indices_obs,1]=x_obs
#     M1[indices_miss,2:(n_miss+1)]=diag(n_miss)
#     v1=rep(0,n);v1[indices_obs]=1
#     M2=M1-matrix(v1,ncol=1)%*%matrix(rep(1,n),nrow=1)%*%M1/n_obs
#     resSvd=svd(M2)
#     #(t(M2)%*%M2)-t(M3)%*%M3
#     M3=diag(resSvd$d)%*%t(resSvd$v)
#     #solve(M3)-resSvd$v%*%diag(1/(resSvd$d))
#     normyres=sqrt(t((t(solve(M3))%*%t(M2)%*%matrix(z,ncol=1)))%*%(t(solve(M3))%*%t(M2)%*%matrix(z,ncol=1)))
#     yres=sqrt(n)*(t(solve(M3))%*%t(M2)%*%matrix(z,ncol=1))/as.numeric(normyres)
#     if(sameBlockWeight){yres=yres/sqrt(p)}
#     # t(yres)%*%yres
#     ures=solve(M3)%*%yres
#     xres=M2%*%ures
#     print("criteria after")
#     print(cov(xres,z))
#     return(xres)
# }

# En refaisant tourner la fonction "à la main". 
# test function
# leb2=function(x_k,missing,z)
# {#leb(x_k=A[[j]][,k],missing,z=Z[,j])
#     
#     # test function
#     # n=100;x_obs_ref=rnorm(n);x_obs_ref=scale(x_obs_ref)
#     # prop=0.2;indices_miss=sample(1:n,round(prop*n))
#     # n_miss=length(indices_miss); n_obs=n-n_miss
#     # n_inc=n_miss+1
#     # indices_obs=(1:n)[!(1:n)%in%indices_miss]
#     # x_obs=x_obs_ref[indices_obs]
#     # z=x_obs_ref+0.05*rnorm(n)
#     print("criteria before")
#     # print(cov(x_k,z))
#     #x_k=scale(x_k)
#     print(cov(x_k,z))
#     x_obs=x_k[!missing]
#     x_mis=x_k[missing]
#     n_obs=length(x_obs)
#     n_mis=length(x_mis)
#     n=n_obs+n_mis
#     M1=matrix(0,n,n_mis+1) # n_mis+1 colonne, n lignes
#     M1[1:n_obs,1]=x_obs
#     M1[(n_obs+1):n,2:(n_mis+1)]=diag(n_mis)
#     M2=M1-matrix(c(rep(1,n_obs),rep(0,n_mis)),ncol=1)%*%matrix(rep(1,n),nrow=1)%*%M1/n_obs
#     resSvd=svd(M2)
#     #(t(M2)%*%M2)-t(M3)%*%M3
#     M3=diag(resSvd$d)%*%t(resSvd$v)
#     #solve(M3)-resSvd$v%*%diag(1/(resSvd$d))
#     normyres=sqrt(t((t(solve(M3))%*%t(M2)%*%matrix(Z[, j],ncol=1)))%*%(t(solve(M3))%*%t(M2)%*%matrix(Z[, j],ncol=1)))
#     yres=sqrt(n)*(t(solve(M3))%*%t(M2)%*%matrix(Z[, j],ncol=1))/as.numeric(normyres)
#     # t(yres)%*%yres
#     ures=solve(M3)%*%yres
#     xres=M2%*%ures
#     # added for a condition = 1/sqrt(p)
#     #xres=xres/sqrt(j)
#     x_k_res=x_k
#     x_k_res[missing]=xres[(n_obs+1):n]
#     x_k_res[!missing]=xres[1:n_obs]
#     print("criteria after")
#     print(cov(x_k_res,z))
#     return(x_k_res)
# }