#' print.cval
#' 
#'@param x result of rgcca_cv function (see \link{rgcca_cv})
#'@param bars among "sd" for standard deviations, "stderr" for standard error (standard deviations divided by sqrt(n), ci for confidence interal, cim for confidence interval of the mean.
#'@param alpha used for confidence interval bars (ci or cim), risk. Default to 0.05
#'@param ... Further print options
#'@export 
#'@examples
#'data("Russett")
#' blocks <- list(
#'     agriculture = Russett[, seq(3)],
#'     industry = Russett[, 4:5],
#'     politic = Russett[, 6:11])
#'     res=rgcca_cv(blocks, type="rgcca",par="tau",par_value=c(0,0.2,0.3),n_cv=1,n_cores=1)
#'    print(res)
print.cval=function(x,bars="sd",alpha=0.05,...)
{
    mean_b=apply(x,1,mean)

    match.arg(bars,c("sd","stderr","ci","cim"))
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

    }

    df=data.frame(config=1:nrow(x),mean=mean_b,inf=inf_b,sup=sup_b)

    optimal_ind=which.min(df[,"mean"])
    optimal_x=df[optimal_ind,"config"]
    optimal_y=df[optimal_ind,"mean"]
    cat(paste0(nrow(x)," configurations were tested, with ", ncol(x)," runs each \n"))
    print(df)
    
    cat(paste("The best configuration was:", rownames(x)[optimal_ind],"for a value of ", round(optimal_y,digits=2)),"\n",...)


}