# #'@export
scaling = function(
    blocks,
    scale = TRUE,
    bias = TRUE,
    sameBlockWeight=TRUE
    ) {

    # if (scale) {
    #     lapply(blocks, function(x)
    #     scale2(x, bias = bias) / sqrt(NCOL(x)))
    # }else{
    #     lapply(blocks, function(x)
    #     scale2(x, scale = FALSE))
    # }
    
    
    if (scale == TRUE) 
    {
        
        A = lapply(blocks, function(x) scale2(x,scale=TRUE, bias = bias)) # le biais indique si on recherche la variance biaisee ou non
        if(sameBlockWeight)
        {
            A = lapply(A, function(x) 
            {
                y=x/sqrt(NCOL(x));
                return(y)
            } 
            )
        }
        # on divise chaque bloc par la racine du nombre de variables pour avoir chaque poids pour le meme bloc
    }
    if (scale == FALSE)
    { 
        
        A = lapply(blocks, function(x) scale2(x, scale=FALSE, bias = bias)) 
        
        if(sameBlockWeight)
        {
            A = lapply(A, function(x) 
            {
                covarMat=cov2(x,bias=bias);
                varianceBloc=sum(diag(covarMat)); 
                res=x/sqrt(varianceBloc)
                attr(res,'scaled:scale')=rep(sqrt(sum(diag(covarMat))),dim(x)[2])
                return(res);
            })
        }
        
    } 
    return(A)
    
}
