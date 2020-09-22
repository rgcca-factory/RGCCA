#' print.cval
#' 
#'@inheritParams plot.cval
#'@param ... Further print options
#'@export 
#'@examples
#'data("Russett")
#' blocks <- list(
#'     agriculture = Russett[, seq(3)],
#'     industry = Russett[, 4:5],
#'     politic = Russett[, 6:11])
#'     res=rgcca_cv(blocks, response=3,type="rgcca",par_type="tau",par_value=c(0,0.2,0.3),
#'     n_run=1,n_cores=1)
#'    print(res)
print.cval=function(x,bars="sd",alpha=0.05,...)
{
    
    cat("Call: ")
    names_call=c("type","nperm","method","tol","scale","scale_block")
    char_to_print=""
    for(name in names_call)
    {
        if(name=="ncomp"){if(length(x$call$ncomp)>1){value=(paste(x$call$ncomp,sep="",collapse=","));value=paste0("c(",value,")")}}
        if(name!="ncomp"){value=x$call[[name]]}
        quo=ifelse(is.character(value)&name!="ncomp","'","")
        vir=ifelse(name==names_call[length(names_call)],"",", ")
        char_to_print=paste(char_to_print,name,'=',quo,value,quo,vir, collapse="",sep="")
    }
    cat(char_to_print)
    
    
    c1s <- round(x$penalties, 4)
    rownames(c1s) = 1:NROW(c1s)
    cat(fill = TRUE)
    cat("Tuning parameters used: ", fill = TRUE)
    print(c1s, quote = FALSE,...)
    cat("\n")
    
    df <- summary.cv(x, bars, alpha)

    optimal_ind=which.min(df[,"Mean RMSE"])
    optimal_x=df[optimal_ind,"Combination"]
    optimal_y=df[optimal_ind,"Mean RMSE"]
    cat(paste0(nrow(x$cv)," configurations were tested. \n"))
    
   cat(paste0("Validation: ",x$call$validation,ifelse(x$call$validation=="kfold", paste0(" with ",x$call$k," folds and ",x$call$n_run," run(s))"),")")),"\n")
    
    print(df)
    
    cat(paste("The best combination was:", paste(round(x$bestpenalties,digits=3),collapse=" "),"for a mean CV criterion (RMSE) of ", round(optimal_y,digits=2)),"\n",...)


}
