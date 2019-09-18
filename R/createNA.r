#' Create a dataset with missing values
#' @param A A list of complete matrices with the same number of lines
#' @param option "bloc" if the structure of missing data is by block or "ponc" if it is random across variables
#' @param pNA A logical value. If scale = TRUE, each column is transformed to have unit variance (default = TRUE).
#' @param nAllRespondant If the option "ponc" is chosen,nAllRespondant corresponds to the number of complete individuals required
#' @return \item{A}{a list of NA}
#' @title Create a list with missing data (for simulation)
#' @export scale2

createNA=function(A,option="block",pNA=0.1,nAllRespondants=4,output="list")
{
  if(length(pNA)==1){pNA=rep(pNA,length(A))}
  if(length(pNA)!=1 & length(pNA)!=length(A)){stop("pNA should be a number between 0 and 1 or a vector of the same size as A ")}
  if(!option%in% c("block","rand","ponc")){stop("option should be chosen as 'block' or 'rand' or 'ponc'")}
  if(is.null(rownames(A[[1]]))){rownames(A[[1]])=paste("S",1:dim(A[[1]])[1],sep="")}
  if(option=="block")
	{
	  # For all the one but the last one, 
	  nbloc=length(A)
	  indToRemove=list()
		for(i in 1:(length(A)-1))
		{
			n=nrow(A[[i]])
			p=ncol(A[[i]])
			nbNAparBloc=round(pNA[i]*n)	
			if(nbNAparBloc!=0)
			{
		  	indToRemove[[i]]=sample(1:n, nbNAparBloc, replace = FALSE)	
		  	 A[[i]][indToRemove[[i]],]=NA
		  }
			else
			{
				indToRemove=NA
			}
		}
	   n=nrow(A[[nbloc]])
	  p=ncol(A[[nbloc]])
	  nbNAparBloc=round(pNA[length(A)]*n)	
	  # The last block can not be a missing line if all other blocks are missing
	  #---------------------------------------------------------------------------
	  W=do.call(cbind,A[1:(nbloc-1)])
    indToTake=which(apply(W,1,function(x){sum(!is.na(x))})==0) #index with missing data for all first blocks
    if(length(indToTake)==0){indicesToUse=1:n}else{indicesToUse=(1:n)[!(1:n)%in%indToTake]}
	  if(nbNAparBloc!=0)
	  {
	    indToRemove[[nbloc]]=sample(indicesToUse, min(length(indicesToUse),nbNAparBloc), replace = FALSE)	
	    A[[nbloc]][indToRemove[[nbloc]],]=NA
	   }
	  else
	  {
	    indToRemove=NA
	  }
    W2=do.call(cbind,A)
    subjectKept=rownames(A[[1]])[which(apply(W2,1,function(x){sum(is.na(x))})==0)]
    
  }
   
	if(option=="rand")
	{
		
		for(i in 1:length(A))
		{
			n=nrow(A[[i]])
			p=ncol(A[[i]])
			nbNA=round(pNA[i]*n*p)
			if(nbNA!=0)
			{
				listeIndicePossible=merge(1:n,1:p)
				indToRemove=sample(1:(n*p), nbNA, replace = FALSE)
				listeIndicesToRemove=listeIndicePossible[indToRemove,]		
				for(j in 1:length(indToRemove)){A[[i]][listeIndicesToRemove[j,"x"],listeIndicesToRemove[j,"y"]]=NA}
			}
			else{indToRemove=NA}
		}
	}
	
	if(option=="ponc") # valeurs ponctuelles en conservant un nombre d'individu complets
	{
		n=nrow(A[[i]])
		allResp=sample(1:n,nAllRespondants)
		restData=(1:n)[!(1:n)%in%allResp]
	
		for(i in 1:length(A))
		{
			p=ncol(A[[i]])
			nbNA=round(pNA[i]*n*p)
			if(nbNA!=0)
			{
				listeIndicePossible=merge(1:n,1:p)
				indToRemove=sample(1:(n*p), nbNA, replace = FALSE)
				listeIndicesToRemove=listeIndicePossible[indToRemove,]		
				for(j in 1:length(indToRemove)){A[[i]][listeIndicesToRemove[j,"x"],listeIndicesToRemove[j,"y"]]=NA}
			}
			else{indToRemove=NA}
		}
		W2=do.call(cbind,A)
		subjectKept=rownames(A[[1]])[which(apply(W2,1,function(x){sum(is.na(x))})==0)]
		
	}
		
  if(output=="list"){res=list(dat=A,naPos=indToRemove,subjectKept=subjectKept)}else{res=do.call(cbind,A)}
	return(res)
}
