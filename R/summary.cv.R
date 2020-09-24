summary_cval <- function(x, bars="quantile") {
    mat_cval=x$cv
    mean_b=apply(mat_cval,1,mean)
    
    match.arg(bars,c("sd","stderr","quantile"))
    if(bars!="none"&&dim(mat_cval)[2]<3){bars=="none"; warning("Standard deviations can not be calculated with less than 3 columns in x")}
    if(bars!="none")
    {
        if(bars=="quantile")
        {
            inf_b=apply(mat_cval,1,function(y){return(quantile(y,0.05))})
            sup_b=apply(mat_cval,1,function(y){return(quantile(y,0.95))})
            lowlim="Quant. 0.05"
            uplim="Quant. 0.95"
        }
        if(bars=="sd")
        {
            inf_b=mean_b-apply(mat_cval,1,sd)
            sup_b=mean_b+apply(mat_cval,1,sd)
            lowlim="Mean - Sd"
             uplim="Mean + Sd"
        }
        if(bars=="stderr")
        {
            inf_b=mean_b-apply(mat_cval,1,function(y){sd(y)/sqrt(length(y))})
            sup_b=mean_b+apply(mat_cval,1,function(y){sd(y)/sqrt(length(y))})
            lowlim="Mean - Std Error"
            uplim="Mean + Std Error"
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
        
    }
    df <- round(data.frame(config=1:nrow(mat_cval),mean=mean_b,inf=inf_b,sup=sup_b), 3)
    if(x$call$type_cv=="regression")
    {
        colnames(df) <- c("Combination", "Mean RMSE", lowlim, uplim)
    }
    if(x$call$type_cv=="classification")
    {
        colnames(df) <- c("Combination", "Mean Error Prediction Rate", lowlim, uplim)
    }
     return(df)
}
