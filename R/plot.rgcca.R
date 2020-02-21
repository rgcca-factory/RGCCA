#' plot 
#' Plots 
#' @title Regularized Generalized Canonical Correlation Analysis (RGCCA) 
#' @param x Result of rgcca function
#' @param text_ind A bolean to represent the individuals with their row names (TRUE)
#' or with circles (FALSE)
#' @param text_var A bolean to represent the variables with their row names (TRUE)
#' or with circles (FALSE)
#' @param title_ind Character for the title of the individual space 
#' @param title_var Character for the title of the variable space 
#' @param ... Further graphical parameters applied to both (individual and variable) spaces
#' @inheritParams plot_ind
#' @inheritParams plot2D
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
plot.rgcca=function(x,resp=rep(1, NROW(x$Y[[1]])),i_block=1,compx=1,compy=2,remove_var=FALSE,text_var=TRUE,text_ind=TRUE,response_name= "Response",no_overlap=FALSE,title_var="Variable correlations with",title_ind= "Sample space",n_mark=100,collapse=FALSE,cex_sub=10,cex_main=14,cex_lab=12,...)
{
    if(is.null(i_block)){i_block=length(x$call$blocks)}
  
    p1<-plot_ind(x,i_block=i_block,compx=compx,compy=compy,cex_sub=cex_sub,cex_main=cex_main,cex_lab=cex_lab,resp=resp,response_name=response_name,text=text_ind,title=title_ind,...)
    p2<-plot_var_2D(x,i_block=i_block,compx=compx,compy=compy,cex_sub=cex_sub,cex_main=cex_main,cex_lab=cex_lab,remove_var=remove_var,text=text_var,no_overlap=no_overlap,title=title_var,n_mark = n_mark,collapse=collapse)
    
   # p3<-plot_ave(x)
   # p4<-plot_network(x)
    #p5<-grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)
   
   titlePlot=toupper(names(x$call$blocks)[i_block])
  #  titlePlot=textGrob("Daily QC: Blue",gp=gpar(fontsize=20,font=3))
     #  p5<-grid.arrange(p1,p2,nrow=1,ncol=2,top = titlePlot)
   p5<-grid.arrange(p1,p2,nrow=1,ncol=2,top = titlePlot)
   
    invisible(p5)
}