#'Print bootstrap
#'@param x results of a bootstrap function
#'@param ... Further arguments in print
#'@export
print.bootstrap=function(x,...)
{
    print(paste0(length(x),"bootstrap(s) were run with the following RGCCA parameters"),...)
    print(x$rgcca$call)
    
}