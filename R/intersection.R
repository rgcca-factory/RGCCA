#' Keeps only subject without missing values
#' @param A A list of matrices
#' @return The intersected list from matrices
#' @title Intersection
#' @examples
#' set.seed(42);X1=matrix(rnorm(35),7,5);
#'set.seed(22);X2=matrix(rnorm(28),7,4);
#'set.seed(2);X3=matrix(rnorm(49),7,7);
#'# usual test 
#'X1[1,]=NA
#'X2[7,1]=NA
#'X2[5,1]=NA
#'A=list(X1,X2)
#'intersection(A=A)
#'# too many subjects with missing values
#'X3[3,1:2]=NA
#'intersection(A=list(X1,X2,X3))
#' @export intersection
intersection=function(A)
{
	A=lapply(A,as.matrix)
#	lapply(A,function(x){if(is.null(rownames(x))){print("one matrix has no colnames")}})
	newList=lapply(A,
		function(x) 
		{
		  if(is.null(rownames(x)))
		  {
		    warnings("No rownames in matrix, they were initialized as S1,...Sn");
		    rownames(x)=paste0("S",1:dim(x)[1])
		  }
			vecNA=apply(x,1,sum);
			if(dim(x)[1]>1)
			{return(x[which(!is.na(vecNA)),])}
			else
			{
				y=x[which(!is.na(x)),];
				rownames(y)=rownames(x)[which(!is.na(x),)];return(y)
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
  	interlist=lapply(newList,function(x){if(!is.null(dim(x))){x[final,]}else{x[final]}})
	}
	else
	{print("less than 3 subjects in the intersection");interlist=NA}
	return(interlist)
}