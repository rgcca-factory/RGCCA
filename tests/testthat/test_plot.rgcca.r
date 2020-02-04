#'# plot.rgcca
#'''
data(Russett)
X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
X_ind = as.matrix(Russett[,c("gnpr","labo")]);
X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
A = list(agri=X_agric,ind=X_ind,polit=X_polit);

res=rgcca(A,type="rgcca")
library(ggplot2)
library(gridExtra)
library(grid)
plot.rgcca(res)

