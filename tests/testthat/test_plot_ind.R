#' Plot the two components of a RGCCA
 coord = lapply(seq(3),
  function(x) matrix(runif(15 * 2, min = -1), 15, 2))
 AVE_X = lapply(seq(3), function(x) runif(2))
 for (i in 1:length(coord))
 row.names(coord[[i]]) = seq(15)
 rgcca_out = list(Y = coord, AVE = list(AVE_X = AVE_X), call = list(blocks = coord, ncomp = rep(2,3)))
 rgcca_out$call$type="rgcca"
 class(rgcca_out) <- "rgcca"
 # Using a superblock
 resp = as.matrix(rep(LETTERS[seq(3)], each = 5))
row.names(resp) = seq(15)

p1<-plot2D ( rgcca_res=rgcca_out,   df=as.data.frame(rgcca_out$Y[[1]]))
p1 <- p1+geom_point(aes(x=1,y=1,colour="blue"))
plot_ind(rgcca_out, resp)
# Using the first block
 resp = as.matrix(runif(15, min=-15, max = 15))
rownames(resp) = seq(15)
#rgcca_out$call$type=
 plot_ind(rgcca_out, resp, 1, 2, 1)
 

# Using plot_ind
 data(Russett)
 blocks = list(agri=Russett[,1:3],ind=Russett[,4:5],polit=Russett[,8:11])
 rgcca_out = rgcca(blocks = blocks, response = 3)
 plot_ind( rgcca_out,resp=1:47,i_block=1)
 
# Using the predict 
data(Russett)
 blocks = list(agri=Russett[,1:3],ind=Russett[,4:5],polit=Russett[,8:11])
 rgcca_out = rgcca(blocks = blocks, response = 3)
 loo <- rgcca_crossvalidation(rgcca_out, n_cores = 2)

 plot_ind(rgcca_out, predicted = loo)
 
 # using a color response vector
 plot_ind(rgcca_out,resp=resp)
 