#'Print bootstrap
#'@param x results of a bootstrap function
#'@param ... Further arguments in print
#'@export
print.bootstrap=function(x,...)
{
    print(paste(dim(x$bootstrap[[1]][[1]])[2],"bootstrap(s) were run"),...)
    ncompmax=min(x$rgcca$call$ncomp)
    for(comp in 1:ncompmax)
    {
        cat(paste("Dimension:",comp,"\n"))
        print(Reduce(rbind,lapply(1:length(x$rgcca$call$blocks),
                                  function(block)
                                  { b=get_bootstrap(b=x,
                                                      i_block=block,
                                                    comp=comp,
                                                    bars="ci",
                                                    display_order =TRUE)
                                  ;return(b)}
            )))
    }
    
}