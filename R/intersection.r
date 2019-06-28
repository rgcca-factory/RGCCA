#' Keeps only subject without missing values
#' @param A A list of matrices
#' @return The intersected list from matrices
#' @title Intersection
#' @examples
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
	final=Reduce(intersect, lapply(newList,function(x){if(!is.null(dim(x))){return(rownames(x))}else{print("one column");return(names(x))}}))
	if(length(final)>3)
	{
  	interlist=lapply(newList,function(x){if(!is.null(dim(x))){x[final,]}else{x[final]}})
	}
	else
	{print("less than 3 subjects in the intersection");interlist=NA}
	return(interlist)
}