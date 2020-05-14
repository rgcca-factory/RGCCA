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
    
    cat("Call: ")
    dput(x$call,control=c())
    cat("\n")
    c1s <- round(x$penalties, 4)
    rownames(c1s) = 1:NROW(c1s)
    cat(fill = TRUE)
    cat("Tuning parameters used: ", fill = TRUE)
    print(c1s, quote = FALSE,...)
    cat("\n")
    
    mat_cval=x$cv
    mean_b=apply(mat_cval,1,mean)

    match.arg(bars,c("sd","stderr","ci","quantile"))
    if(bars!="none"&&dim(mat_cval)[2]<3){bars=="none"; warning("Standard deviations can not be calculated with less than 3 columns in x")}
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
        if(bars=="cim")
        {
            stat=qt(1-alpha/2,df=dim(mat_cval)[2]-1)
            inf_b=mean_b-apply(mat_cval,1,function(y){stat*sd(y)/sqrt(length(y))})
            sup_b=mean_b+apply(mat_cval,1,function(y){stat*sd(y)/sqrt(length(y))})
        }
        if(bars=="ci")
        {
            stat=qt(1-alpha/2,df=dim(mat_cval)[2]-1)
            inf_b=mean_b-apply(mat_cval,1,function(y){stat*sd(y)})
            sup_b=mean_b+apply(mat_cval,1,function(y){stat*sd(y)})
        }

    }

    df=data.frame(config=1:nrow(mat_cval),mean=mean_b,inf=inf_b,sup=sup_b)

    optimal_ind=which.min(df[,"mean"])
    optimal_x=df[optimal_ind,"config"]
    optimal_y=df[optimal_ind,"mean"]
    cat(paste0(nrow(mat_cval)," configurations were tested. \n"))
    
   cat(paste0("Validation: ",x$call$validation,ifelse(x$call$validation=="kfold", paste0(" with ",x$call$k," folds and ",x$call$n_cv," runs)"),")")),"\n")
    
    print(df)
    
    cat(paste("The best configuration was:", paste(round(x$bestpenalties,digits=3),collapse=" "),"for a value of ", round(optimal_y,digits=2)),"\n",...)


}