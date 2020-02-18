#' plot 
#' Plots 
#' @title Regularized Generalized Canonical Correlation Analysis (RGCCA) 
#' @param ... parameters from plot_ind or plot_var_2D
#' @param x Result of rgcca function
#' @inheritParams plot_ind
#' @inheritParams plot_var_2D
#' @examples
#' data(Russett)
#' X_agric =as.matrix(Russett[,c("gini","farm","rent")]);
#' X_ind = as.matrix(Russett[,c("gnpr","labo")]);
#' X_polit = as.matrix(Russett[ , c("demostab", "dictator")]);
#' A = list(X_agric, X_ind, X_polit);
#' C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3);
#' resRgcca=rgcca(blocks=A,connection=C,
#' tau=rep(1,3),ncomp=rep(2,3),superblock=FALSE)
#' plot(resRgcca)
#' @importFrom gridExtra grid.arrange
#' @export
plot.rgcca=function(x,resp=rep(1, NROW(x$Y[[1]])),i_block=1,compx=1,compy=2,remove_var=FALSE,...)
{
    if(is.null(i_block)){i_block=length(x$call$blocks)}
  
    p1<-plot_ind(x,i_block=i_block,compx=compx,compy=compy,cex_sub=10,cex_main=14,cex_lab=12,resp=resp,...)
    p2<-plot_var_2D(x,i_block=i_block,compx=compx,compy=compy,cex_sub=10,cex_main=14,cex_lab=12,remove_var=remove_var,...)
    
   # p3<-plot_ave(x)
   # p4<-plot_network(x)
    #p5<-grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)
   
   titlePlot=toupper(names(x$call$blocks)[i_block])
  #  titlePlot=textGrob("Daily QC: Blue",gp=gpar(fontsize=20,font=3))
     #  p5<-grid.arrange(p1,p2,nrow=1,ncol=2,top = titlePlot)
   p5<-grid.arrange(p1,p2,nrow=1,ncol=2,top = titlePlot)
   
    return(p5)
}