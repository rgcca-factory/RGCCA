listOfMatrices=function(A)
{#
    #
    #
    B=A
    i=1;continue=TRUE
    while(continue & i<length(A)+1)
    {
        if(!is.null(rownames(A[[i]])))
        {
            rowNames=rownames(A[[i]])
            continue=FALSE
        }
        else{i<-i+1}
    }
    for(i in 1:length(A))
    {
        if(is.null(dim(A[[i]])))
        {
            print("one !")
            B[[i]]=matrix(A[[i]],ncol=1)
            rownames(B[[i]])=names(A[[i]])
            colnames(B[[i]])=names(A)[i]
        }
        else
        { 
            if(dim(A[[i]])[2]==1)
            {
                if(is.null(rownames(A[[i]])))
                {
                    warning(paste("No rownames on the ",i,"block. The same rownames as the first block with rownames were attributed"))
                }
                rownames(B[[i]])=rowNames
                colnames(B[[i]])=names(A)[i]  
            }
            else
            {
                rownames(B[[i]])=rownames(A[[i]])
                colnames(B[[i]])=colnames(A[[i] ])
            }
            
        }
       
    }
    names(B)=names(A)
    return(B)
}
