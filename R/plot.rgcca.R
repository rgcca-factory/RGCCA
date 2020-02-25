#' plot 
#' Plots 
#' @title Regularized Generalized Canonical Correlation Analysis (RGCCA) 
#' @param x Result of rgcca function
#' @param type Type among c("ind","var","both","ave","top","network"). "ind" for individual graph, "var" for variable graph, "both" for both, "ave" for the variance average in each block, "net"for network
#' @param text_ind A bolean to represent the individuals with their row names (TRUE)
#' or with circles (FALSE)
#' @param text_var A bolean to represent the variables with their row names (TRUE)
#' or with circles (FALSE)
#' @param title_ind Character for the title of the individual space 
#' @param title_var Character for the title of the variable space 
#' @param colors representing a vector of the colors used in the graph. Either a vector of integers (each integer corresponding to a color) or of characters corresponding to names of colors (as "blue",see colors()) or RGB code ("#FFFFFF").
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
plot.rgcca=function(x,type="both",resp=rep(1, NROW(x$Y[[1]])),i_block=1,i_block_y=i_block,compx=1,compy=2,remove_var=FALSE,text_var=TRUE,text_ind=TRUE,response_name= "Response",no_overlap=FALSE,title=NULL,title_var="Variable correlations with",title_ind= "Sample space",n_mark=100,collapse=FALSE,cex=1,cex_sub=10,cex_main=14,cex_lab=12,colors=NULL,...)
{
    match.arg(type,c("ind","var","both","ave","top","network"))
    
    if(type=="both")
    {
        if(is.null(i_block)){i_block=length(x$call$blocks)}
        p1<-plot_ind(x,i_block=i_block,i_block_y=i_block_y,compx=compx,compy=compy,cex_sub=cex_sub,cex_main=cex_main,cex_lab=cex_lab,resp=resp,response_name=response_name,text=text_ind,title=title_ind,colors=colors,...)
        p2<-plot_var_2D(x,i_block=i_block,compx=compx,compy=compy,cex_sub=cex_sub,cex_main=cex_main,cex_lab=cex_lab,remove_var=remove_var,text=text_var,no_overlap=no_overlap,title=title_var,n_mark = n_mark,collapse=collapse,colors=colors,...)
        titlePlot=toupper(names(x$call$blocks)[i_block])
        p5<-grid.arrange(p1,p2,nrow=1,ncol=2,top = titlePlot)
        plot(p5)
    }
    if(type=="var")
    {
       p5 <- plot_var_2D(x,i_block=i_block,compx=compx,compy=compy,cex_sub=cex_sub,cex_main=cex_main,cex_lab=cex_lab,remove_var=remove_var,text=text_var,no_overlap=no_overlap,title=title_var,n_mark = n_mark,collapse=collapse,colors=colors)
        plot(p5)
     }
    if(type=="ind")
    {
        p5<-plot_ind(x,i_block=i_block,i_block_y=i_block_y,compx=compx,compy=compy,cex_sub=cex_sub,cex_main=cex_main,cex_lab=cex_lab,resp=resp,response_name=response_name,text=text_ind,title=title_ind,colors=colors,...)
        plot(p5)
     }
    if(type=="ave")
    {
        if(is.null(title)){title="Average Variance Explained"}
        p5 <- plot_ave (x,
            cex = cex,
            title = title,
            colors = colors,
            ...)
        plot(p5)
    }
    if(type=="network")
    {
        if(is.null(title)){title=paste0("Common rows between blocks : ",
                                        NROW(x$call$blocks[[1]]))}
        plot_network (
            x, 
            title = title)
        p5<-NULL
    }
      

   # p3<-plot_ave(x)
   # p4<-plot_network(x)
   #p5<-grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)
   
   
    invisible(p5)
}