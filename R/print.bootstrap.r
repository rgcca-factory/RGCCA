#'Print bootstrap
#'@param x results of a bootstrap function
#'@param ... Further arguments in print
#'@export
print.bootstrap=function(x,...)
{
    print(paste(length(x),"bootstrap(s) were run with the following RGCCA parameters"),...)
    b=get_bootstrap(b=x,bars="ci")
    ncompmax=min(x$rgcca$call$ncomp)
    res=lapply(1:ncompmax, function(comp)
    {
        
        return(Reduce(rbind,lapply(1:length(x$rgcca$call$blocks),function(block)
            { b=get_bootstrap(b=x,i_block=block,comp=comp,bars="ci");return(b)}
            )))
    })
    
}