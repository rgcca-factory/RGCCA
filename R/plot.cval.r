#' plot.cval
#' 
#'@param x result of rgcca_cv function (see \link{rgcca_cv})
#'@param bars among "sd" for standard deviations, "stderr" for standard error (standard deviations divided by sqrt(n), ci for confidence interal, cim for confidence interval of the mean.
#'@param alpha used for confidence interval bars (ci or cim), risk. Default to 0.05
#'@param ... Further plot options
#'@export 
#'@examples
#'data("Russett")
#' blocks <- list(
#'     agriculture = Russett[, seq(3)],
#'     industry = Russett[, 4:5],
#'     politic = Russett[, 6:11])
#'     res=rgcca_cv(blocks, type="rgcca",par="tau",par_value=c(0,0.2,0.3),n_cv=1,n_cores=1)
#'    plot(res)
plot.cval=function(x,bars="sd",alpha=0.05,...)
{
    config <- NULL -> y
    match.arg(bars,c("sd","stderr","ci","cim","points"))
    mean_b=apply(x,1,mean)
    main=paste0("Mean CV criterion according to the configuration set\n (",ncol(x)," runs)")
    if(bars!="none"&&dim(x)[2]<3){bars=="none"; warning("Standard deviations can not be calculated with less than 3 columns in x")}
    if(bars!="none")
    {
        if(bars=="sd")
        {
            inf_b=mean_b-apply(x,1,sd)
            sup_b=mean_b+apply(x,1,sd)
        }
        if(bars=="stderr")
        {
            inf_b=mean_b-apply(x,1,function(y){sd(y)/sqrt(length(y))})
            sup_b=mean_b+apply(x,1,function(y){sd(y)/sqrt(length(y))})
        }
        if(bars=="cim")
        {
            stat=qt(1-alpha/2,df=dim(x)[2]-1)
            inf_b=mean_b-apply(x,1,function(y){stat*sd(y)/sqrt(length(y))})
            sup_b=mean_b+apply(x,1,function(y){stat*sd(y)/sqrt(length(y))})
        }
        if(bars=="ci")
        {
            stat=qt(1-alpha/2,df=dim(x)[2]-1)
            inf_b=mean_b-apply(x,1,function(y){stat*sd(y)})
            sup_b=mean_b+apply(x,1,function(y){stat*sd(y)})
        }
        if(bars=="points")
        {
            inf_b=apply(x,1,function(y){min(y)})
            sup_b=apply(x,1,function(y){max(y)})
        }
        
    }
    
    df=data.frame(config=1:nrow(x),mean=mean_b,inf=inf_b,sup=sup_b)
    p<- ggplot(data=df,aes(x=config,y=mean))+geom_point()
    if(bars!="none"&& bars!="points")
    {
        p<-p+geom_segment(data=df,aes(x=config,y=inf_b,xend=config,yend=sup_b),colour="grey")
    }
    if(bars=="points")
    {
        df2=data.frame(config=rep(1:nrow(x),ncol(x)),y=as.vector(x))
        p<-p+geom_point(data=df2,aes(x=config,y=y),colour="grey")
    }
    optimal_x=df[which.min(df[,"mean"]),"config"]
    optimal_y=df[which.min(df[,"mean"]),"mean"]
    decalage= (max(df[,"mean"])-min(df[,"mean"]))/10
    p<- p+geom_point(x=optimal_x,y=optimal_y,colour="green")
   # p<-p+geom_text(label=rownames(x)[which.min(mean_b)],x=optimal_x,y=optimal_y+decalage,colour="green")
    p<-p+ggtitle(main)
#    p<- p+ scale_x_continuous(breaks=1:nrow(x),  labels=rownames(x))
#    p<-p + theme(axis.text.x = element_text(angle=45))
    plot(p)
    return(df)
}