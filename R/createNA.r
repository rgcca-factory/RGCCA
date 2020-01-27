#' Create a dataset with missing values
#' @param A A list of complete matrices with the same number of lines
#' @param output list by default, else, can return the concatened list.
#' @param seed if NULL, the creation of missing data is totally random (and consequently unreproducible), if numeric value, select randomness according to a seed (and reproducible).
#' @param option "bloc" if the structure of missing data is by block or "ponc" if it is random across variables
#' @param pNA A logical value. If scale = TRUE, each column is transformed to have unit variance (default = TRUE).
#' @param nAllRespondants If the option "ponc" is chosen,nAllRespondant corresponds to the number of complete individuals required
#' @return \item{A}{a list of NA}
#' @title Create a list with missing data (for simulation)
#' @export
#' @examples 
#' data(Russett)
#' X1=Russett[,1:3]
#' X2=Russett[,4:5]
#' X3=Russett[,8:11]
#' A=list(agri=X1,ind=X2,polit=X3)
#' createNA(A,pNA=0.2)
createNA=function(A,option="block",pNA=0.1,nAllRespondants=4,output="list",seed=NULL)
{
   
  A= listOfMatrices(A)
  if(length(pNA)==1){pNA=rep(pNA,length(A))}
  if(length(pNA)!=1 & length(pNA)!=length(A)){stop("pNA should be a number between 0 and 1 or a vector of the same size as A ")}
  if(!option%in% c("block","rand","ponc","byvar")){stop("option should be chosen as 'block' or 'rand' or 'ponc' or 'byvar'")}
  
    if(is.null(rownames(A[[1]]))){rownames(A[[1]])=paste("S",1:dim(A[[1]])[1],sep="")}
    if(is.list(pNA)){warnings("The percentage of missing data is chosen by variable. 'byvar' option is chosen")}
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
			 if(!is.null(seed)){set.seed(seed+i)}
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
	      if(!is.null(seed)){set.seed(seed+length(A))}
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
				if(!is.null(seed)){set.seed(seed+i)}
				indToRemove=sample(1:(n*p), nbNA, replace = FALSE)
				listeIndicesToRemove=listeIndicePossible[indToRemove,]		
				for(j in 1:length(indToRemove)){A[[i]][listeIndicesToRemove[j,"x"],listeIndicesToRemove[j,"y"]]=NA}
			}
			else{indToRemove=NA}
		}
	}
	
	if(option=="ponc") # valeurs ponctuelles en conservant un nombre d'individu complets
	{
	  n=nrow(A[[1]])
		allResp=sample(1:n,nAllRespondants)
		restData=(1:n)[!(1:n)%in%allResp]
	
		for(i in 1:length(A))
		{
			p=ncol(A[[i]])
			nbNA=round(pNA[i]*n*p)
			if(nbNA!=0)
			{
				listeIndicePossible=merge(restData,1:p)
				if(!is.null(seed)){set.seed(seed+i)}
				indToRemove=sample(1:dim(listeIndicePossible)[1], nbNA, replace = FALSE)
				listeIndicesToRemove=listeIndicePossible[indToRemove,]		
				for(j in 1:length(indToRemove)){A[[i]][listeIndicesToRemove[j,"x"],listeIndicesToRemove[j,"y"]]=NA}
			}
			else{indToRemove=NA}
		}
		W2=do.call(cbind,A)
		subjectKept=rownames(A[[1]])[which(apply(W2,1,function(x){sum(is.na(x))})==0)]
	}
    if(option=="byvar")
    {
        
        indToRemove=list()
        if(length(pNA)!=length(A)){stop("If pNA is a list, it should have the same length as A")}
        for(i in 1:length(A))
        {
            n=nrow(A[[i]])
            
            indToRemove[[i]]=list()
            if(length(pNA[[i]])!=dim(A[[i]])[2]){stop(paste("pNA[",i,"] should have the same length that the number of variable in A[",i,"]"))}
            for(j in 1:dim(A[[i]])[2])
            {
               
                nbNA=round(pNA[[i]][j]*n)
                indToRemove[[i]][[j]]=sample(1:n,nbNA)
             
                A[[i]][indToRemove[[i]][[j]],j]=NA
            }
        }
        W2=do.call(cbind,A)
        subjectKept=rownames(A[[1]])[which(apply(W2,1,function(x){sum(is.na(x))})==0)]
    }
		
  if(output=="list"){res=list(dat=A,naPos=indToRemove,subjectKept=subjectKept)}else{res=do.call(cbind,A)}
	return(res)
}
