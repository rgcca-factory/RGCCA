 df = as.data.frame(matrix(runif(20*2, min = -1), 20, 2))
 AVE = lapply(seq(4), function(x) runif(2))
 rgcca_out = list(AVE = list(AVE_X = AVE), call = list(type = "rgcca"))
 plot2D(rgcca_out, df, "Samples", rep(c("a","b"), each=10), "Response")
 data(Russett)
 blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
 politic = Russett[, 6:11] )
 rgcca_out = rgcca(blocks,ncomp=2)
 plot2D(rgcca_out, df)




 X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
 X_ind = as.matrix(Russett[,c("gnpr","labo")]);
 X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
 A = list(agri=X_agric,ind=X_ind,polit=X_polit)
