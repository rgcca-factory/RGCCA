#'Plot cross-validation
#'
#'Plot a cross-validation object (tuning RGCA parameters in 'supervised' mode).
#'The parameters tuned for maximizing RMSE is displayed in the title. In
#'abscissa, the index of combination (number corresponding to the rownames of
#''penalties' in the tuning parameters object). In ordinate, the average of the
#'permuted RMSE criterion. The best parameters are in red by default.
#'@inheritParams plot2D
#'@inheritParams plot_permut_2D
#' @inheritParams plot.rgcca
#'@param x A rgcca_cv object (see \link{rgcca_cv})
#'@export
#'@examples
#'data("Russett")
#' blocks <- list(
#'     agriculture = Russett[, seq(3)],
#'     industry = Russett[, 4:5],
#'     politic = Russett[, 6:11])
#'     res=rgcca_cv(blocks, response=3,method="rgcca",par_type="tau",par_value=c(0,0.2,0.3),
#'     n_run=1,n_cores=1)
#'    plot(res)
#'@importFrom ggplot2 ggplot
plot.cval=function(x, bars="sd", cex = 1, cex_main = 14 * cex, cex_sub = 10 * cex, cex_lab = 10 * cex, colors = c("red", "black"), ...)
{

    stopifnot(is(x, "cval"))
    match.arg(bars,c("quantile","sd","stderr","points"))
    for (i in c("cex", "cex_main", "cex_sub", "cex_lab"))
        check_integer(i, get(i))
    check_colors(colors)
    if (length(colors) < 2)
        colors <- rep(colors, 2)

    configurations <- NULL -> y
    mat_cval=x$cv

    mean_b=apply(mat_cval,1,mean)

    if(x$call$type_cv=="regression")
    {
        main=paste0("RMSE according to the combinations \n (",
                    x$call$validation,ifelse(x$call$validation=="kfold",
                                             paste0(": with ",x$call$k," folds", ifelse(x$call$n_run>1,paste0(" and ",x$call$n_run," run",ifelse(x$call$n_run==1,"","s")),""),")"),
                                             ")\n "))
        main=paste0(main,"\nbest parameter: ",
                    paste(round(x$bestpenalties,digits=2), collapse = ", "))
        y_lab="Mean RMSE"
    }
    if(x$call$type_cv=="classification")
    {
        main=paste0("Cross-validated prediction error according to the combinations \n (",
                    x$call$validation,
                    ifelse(x$call$validation=="kfold",
                           paste0(": with ",x$call$k," folds", ifelse(x$call$n_run>1,
                                                                      paste0(" and ",x$call$n_run," run",ifelse(x$call$n_run==1,"","s")),""),")"),
                                        ")\n "))
        main=paste0(main,"\nbest parameter: ",
                    paste(round(x$bestpenalties,digits=2), collapse = ", "))
        y_lab="Mean error rate"
    }

    axis <- function(margin){
        element_text(
            size = cex_lab * 0.75,
            face = "italic",
            margin=margin
        )
    }
    if(bars!="none"&&dim(mat_cval)[2]<3){bars=="none"; warning("Standard deviations can not be calculated with less than 3 columns in mat_cval")}
    if(bars!="none")
    {
        if(bars=="quantile")
        {
            inf_b=apply(mat_cval,1,function(y){return(quantile(y,0.05))})
            sup_b=apply(mat_cval,1,function(y){return(quantile(y,0.95))})
        }
        if(bars=="sd")
        {
            inf_b=mean_b-apply(mat_cval,1,sd)
            sup_b=mean_b+apply(mat_cval,1,sd)
        }
        if(bars=="stderr")
        {
            inf_b=mean_b-apply(mat_cval,1,function(y){sd(y)/sqrt(length(y))})
            sup_b=mean_b+apply(mat_cval,1,function(y){sd(y)/sqrt(length(y))})
        }
        # if(bars=="cim")
        # {
        #     stat=qt(1-alpha/2,df=dim(mat_cval)[2]-1)
        #     inf_b=mean_b-apply(mat_cval,1,function(y){stat*sd(y)/sqrt(length(y))})
        #     sup_b=mean_b+apply(mat_cval,1,function(y){stat*sd(y)/sqrt(length(y))})
        # }
        # if(bars=="ci")
        # {
        #     stat=qt(1-alpha/2,df=dim(mat_cval)[2]-1)
        #     inf_b=mean_b-apply(mat_cval,1,function(y){stat*sd(y)})
        #     sup_b=mean_b+apply(mat_cval,1,function(y){stat*sd(y)})
        # }
        if(bars=="points")
        {
            inf_b=apply(mat_cval,1,function(y){min(y)})
            sup_b=apply(mat_cval,1,function(y){max(y)})
        }

    }

    df=data.frame(configurations=1:nrow(mat_cval),mean=mean_b,inf=inf_b,sup=sup_b)
    p<- ggplot(data=df,aes(x=configurations,y=mean))+geom_point(colour=colors[2])+theme_classic() + xlab("Combinations")+ylab(y_lab)
    if(bars!="none"&& bars!="points")
    {
        p<-p+geom_segment(data=df,aes(x=configurations,y=inf_b,xend=configurations,yend=sup_b),colour=colors[2])
    }
    if(bars=="points")
    {
        df2=data.frame(configurations=rep(1:nrow(mat_cval),ncol(mat_cval)),y=as.vector(mat_cval))
        p<-p+geom_point(data=df2,aes(x=configurations,y=y),colour=colors[2])
    }
    optimal_x=df[which.min(df[,"mean"]),"configurations"]
    optimal_y=df[which.min(df[,"mean"]),"mean"]
    decalage= (max(df[,"mean"])-min(df[,"mean"]))/10
    p<- p+geom_point(x=optimal_x,y=optimal_y,colour= colors[1]) +
        geom_segment(data=df,aes(x=optimal_x,y=inf_b[optimal_x],xend=optimal_x,yend=sup_b[optimal_x]),colour= colors[1], size = 0.4)
   # p<-p+geom_text(label=rownames(x)[which.min(mean_b)],x=optimal_x,y=optimal_y+decalage,colour="green")
    p<-p+ggtitle(main)+theme_perso(cex, cex_main, cex_sub) +
        theme(
            axis.text = element_text(size = 10, face = "bold"),
            axis.title.y = axis(margin(0, 20, 0, 0)),
            axis.title.x = axis(margin(20, 0, 0, 0)),
            axis.line = element_line(size = 0.5),
            axis.ticks  = element_line(size = 0.5),
            axis.ticks.length = unit(2, "mm"),
            legend.position = "none"
        )
#    p<- p+ scale_x_continuous(breaks=1:nrow(x),  labels=rownames(x))
#    p<-p + theme(axis.text.x = element_text(angle=45))
    return(p)
}
