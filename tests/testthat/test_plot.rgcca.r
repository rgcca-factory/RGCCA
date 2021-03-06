#' plot.rgcca
#'''
data(Russett)
X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
A = list(agri=X_agric,ind=X_ind,polit=X_polit);

res=rgcca(A,method="rgcca",ncomp=2)
library(ggplot2)
library(gridExtra)
library(grid)
plot(res,colors="blue")
plot.rgcca(res,type="var")
plot(res,type="ind")
plot(res,type="ave",colors=c("blue","red"))
plot(res,type="network")


res=rgcca(A,method="rgcca",superblock=TRUE,ncomp=2)
plot(res,type="ind")
plot(res,type="var")
plot(res,type="ave",colors=c("blue","red"))
plot(res,type="network")


# Response
#------------
res=rgcca(A,method="rgcca",ncomp=2)
library(ggplot2)
library(gridExtra)
library(grid)
plot(res)
plot(res)
plot(res,type="both")
plot(res,type="var")
plot(res,type="ind",colors=c("blue","green"),resp=A[[3]][,1])
plot(res,type="ave",colors=c("blue","red"))
plot(res,type="network")




vec_colors=c(rep(letters[1:9],5),"a","b")
names(vec_colors)=rownames(A[[1]])
plot(res,type="ind",resp=vec_colors)


vec_colors=c(rep(letters[1:11],4),"a","b","c")
names(vec_colors)=rownames(A[[1]])
plot(res,type="ind",resp=vec_colors,colors=rainbow(11))
