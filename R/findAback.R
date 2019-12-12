#' findAback
#' finds the original values of A (scaled)
#'@param A results of scale
#'@example A=
findAback=function(A,sameBlockWeight=TRUE)
{
    Ares=lapply(A,function(x)
        {
            n=dim(x)[1]
            standardDev= matrix(rep(attributes(x)$'scaled:scale',n),dim(x)[1],dim(x)[2],byrow=TRUE)
            means=matrix(rep( attributes(x)$'scaled:center',dim(x)[1]),dim(x)[1],dim(x)[2],byrow=TRUE)
            if(!sameBlockWeight)
            {
                resx=x*(standardDev)+means
            }
            if(sameBlockWeight)
            {
                resx=x*(standardDev)*sqrt(NCOL(x))+means
            }
        })
   return(Ares)
}