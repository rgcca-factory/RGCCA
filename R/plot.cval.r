#' plot.cval
#'
#'@inheritParams plot2D
#'@param x A rgcca_cv object (see \link{rgcca_cv})
#'@param bars A character among "sd" for standard deviations, "stderr" for standard error (standard deviations divided by sqrt(n), "points" or "quantile" for the 0.05-0.95 quantiles
# #'@param alpha A numeric value giving the risk for the confidence interval bars (ci or cim).
#'@param ... Further plot options
#'@export
#'@examples
#'data("Russett")
#' blocks <- list(
#'     agriculture = Russett[, seq(3)],
#'     industry = Russett[, 4:5],
#'     politic = Russett[, 6:11])
#'     res=rgcca_cv(blocks, response=3,type="rgcca",par_type="tau",par_value=c(0,0.2,0.3),
#'     n_run=1,n_cores=1)
#'    plot(res)
#'@importFrom ggplot2 ggplot
plot.cval=function(x,bars="sd",cex = 1, cex_main = 14 * cex, cex_sub = 10 * cex,...)
{
    configurations <- NULL -> y
    mat_cval=x$cv
    match.arg(bars,c("quantile","sd","stderr","points"))
    mean_b=apply(mat_cval,1,mean)
    main=paste0("RMSE according to the combinations \n (",x$call$validation,ifelse(x$call$validation=="kfold", paste0(": with ",x$call$k," folds", ifelse(x$call$n_run>1,paste0(" and ",x$call$n_run," run",ifelse(x$call$n_run==1,"","s")),""),")"),";,\n "))
    main=paste0(main,"\nbest value, in green : ",
            paste(round(x$bestpenalties,digits=2), collapse = ", "))
            
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
    p<- ggplot(data=df,aes(x=configurations,y=mean))+geom_point()+theme_classic() +xlab("Combinations")+ylab("Mean RMSE")
    if(bars!="none"&& bars!="points")
    {
        p<-p+geom_segment(data=df,aes(x=configurations,y=inf_b,xend=configurations,yend=sup_b),colour="grey")
    }
    if(bars=="points")
    {
        df2=data.frame(configurations=rep(1:nrow(mat_cval),ncol(mat_cval)),y=as.vector(mat_cval))
        p<-p+geom_point(data=df2,aes(x=configurations,y=y),colour="grey")
    }
    optimal_x=df[which.min(df[,"mean"]),"configurations"]
    optimal_y=df[which.min(df[,"mean"]),"mean"]
    decalage= (max(df[,"mean"])-min(df[,"mean"]))/10
    p<- p+geom_point(x=optimal_x,y=optimal_y,colour="green")
   # p<-p+geom_text(label=rownames(x)[which.min(mean_b)],x=optimal_x,y=optimal_y+decalage,colour="green")
    p<-p+ggtitle(main)+theme_perso(cex, cex_main, cex_sub)
#    p<- p+ scale_x_continuous(breaks=1:nrow(x),  labels=rownames(x))
#    p<-p + theme(axis.text.x = element_text(angle=45))
    return(p)
}
