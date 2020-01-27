#' plot 
#' Plots 
#' @title Regularized Generalized Canonical Correlation Analysis (RGCCA) 
#' @param ... parameters from plot_ind or plot_var_2D
#' @param resRgcca Result of rgcca.analyze function
#' @examples
#' data(Russett)
#' X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
#' X_ind = as.matrix(Russett[,c("gnpr","labo")]);
#' X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
#' A = list(X_agric, X_ind, X_polit);
#' C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3);
#' resRgcca=rgcca.analyze(blocks=A,connection=C,
#' tau=rep(1,3),ncomp=rep(2,3),superblock=FALSE)
# plot.rgcca.analyze(resRgcca,i_block=1,nrow=2,ncol=1)
#' @importFrom gridExtra grid.arrange
#' @export
plot.rgcca.analyze=function(resRgcca,...)
{
    p1<-plot_ind(resRgcca,...)
    p2<-plot_var_2D(resRgcca,...)
    
   # p3<-plot_ave(resRgcca)
   # p4<-plot_network(resRgcca)
    #p5<-grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)
       p5<-grid.arrange(p1,p2,nrow=2,ncol=1,cex=1.2)
    return(p5)
}