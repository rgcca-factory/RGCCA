# Keeps only subject without missing values
# @inheritParams rgccad
# @return The intersected list from matrices
# @title intersection_list
# @examples
# set.seed(42);X1=matrix(rnorm(35),7,5);
# set.seed(22);X2=matrix(rnorm(28),7,4);
# set.seed(2);X3=matrix(rnorm(49),7,7);
## usual test 
#X1[1,]=NA
#X2[7,1]=NA
#X2[5,1]=NA
#A=list(X1,X2)
#intersection_list(A=A)
## too many subjects with missing values
#X3[3,1:2]=NA
#intersection_list(A=list(X1,X2,X3))
# @export intersection_list
intersection_list=function(A)
{
	A=lapply(A,as.matrix)
    centering=lapply(A,function(x){attributes(x)$'scaled:center'})
    scaling=lapply(A,function(x){attributes(x)$'scaled:scale'})
    
#	lapply(A,function(x){if(is.null(rownames(x))){print("one matrix has no colnames")}})
	newList=lapply(A,
		function(x) 
		{
		  if(is.null(rownames(x)))
		  {
		    warnings("No rownames in matrix, they were initialized as S1,...Sn \n");
		    rownames(x)=paste0("S",1:dim(x)[1])
		  }
			vecNA=apply(x,1,sum);
			if(dim(x)[1]>1)
			{return(x[which(!is.na(vecNA)),,drop=FALSE])}
			else
			{
				y=x[which(!is.na(x)),,drop=FALSE];
				rownames(y)=rownames(x)[which(!is.na(x),)];
		
				return(y)
			}
		}
		)
	final = Reduce(intersect, lapply(newList, function(x) {
	    if (!is.null(dim(x))) {
	        return(rownames(x))
	    } else{
	        return(names(x))
	    }
	}))
	if(length(final)>3)
	{
  	interlist=lapply(newList,function(x){if(!is.null(dim(x))){x[final,,drop=FALSE]}else{x[final,,drop=FALSE]}})
	}
	else
	{print("less than 3 subjects in the intersection_list");interlist=NA}
    for (i in 1:length(interlist))
    {
        if(!is.null(centering[[i]]))
        {
            attr(interlist[[i]],"scaled:center")=centering[[i]]
        }
        if(!is.null(scaling[[i]]))
        {
            attr(interlist[[i]],"scaled:scale")=scaling[[i]]
        }
    }
	return(interlist)
}
