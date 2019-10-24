source("")
n=100;x_obs_ref=rnorm(n);
prop=0.2;
x_obs_ref=scale(x_obs_ref)
z=x_obs_ref+0.05*rnorm(n)
indices_miss=sample(1:n,round(prop*n))
missing=1:n%in%indices_miss

n_miss=length(indices_miss); 
n_obs=n-n_miss
n_inc=n_miss+1
indices_obs=(1:n)[!missing]
x_obs=x_obs_ref[indices_obs]
n= n_obs+n_miss
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
yres=sqrt(n)*(t(solve(M3))%*%t(M2)%*%matrix(z,ncol=1))/as.numeric(normyres)
# t(yres)%*%yres
ures=solve(M3)%*%yres
xres=M2%*%ures

# sans la fonction :  
x_knew=xres # fonctionne
plot(indices_obs,x_obs_ref[indices_obs],pch=16,col="blue")
points(indices_miss,x_obs_ref[indices_miss],pch=16,col="black")
points(indices_obs,x_knew[indices_obs],pch=1,col="blue")
points(indices_miss,x_knew[indices_miss],pch=1,col="black")
plot(x_obs_ref[indices_obs],x_knew[indices_obs])
points(x_obs_ref[indices_miss],x_knew[indices_miss],pch=15,col="red")

# avec la fonction leb
x_knew=leb(x_k=x_obs_ref,missing=missing,z=z,sameBlockWeight=TRUE,weight=2)
x_knew=leb(x_k=x_obs_ref,missing=missing,z=z,sameBlockWeight=FALSE)
plot(indices_obs,x_obs_ref[indices_obs],pch=16,col="blue")
points(indices_miss,x_obs_ref[indices_miss],pch=16,col="black")
points(indices_obs,x_knew[indices_obs],pch=1,col="blue")
points(indices_miss,x_knew[indices_miss],pch=1,col="black")
lm(x_obs_ref[indices_obs]~x_knew[indices_obs])
plot(x_obs_ref[indices_obs],x_knew[indices_obs])
points(x_obs_ref[indices_miss],x_knew[indices_miss],pch=15,col="red")

