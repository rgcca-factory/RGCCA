#' plot.rgcca
#'''
data(Russett)
X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
A = list(agri=X_agric,ind=X_ind,polit=X_polit);

res=rgcca(A,type="rgcca",ncomp=2)
library(ggplot2)
library(gridExtra)
library(grid)
plot(res,colors="blue")
plot.rgcca(res,type="var")
plot(res,type="ind")
plot(res,type="ave",colors=c("blue","red"))
plot(res,type="network")

