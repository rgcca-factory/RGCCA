#' Create a dataset with missing values
#' @param blocks A list of complete matrices with the same number of lines
#' @param output list by default, else, can return the concatened list.
#' @param seed if NULL, the creation of missing data is totally random (and consequently unreproducible), if numeric value, select randomness according to a seed (and reproducible).
#' @param typeNA "bloc" if the structure of missing data is by block or "ponc" if it is random across variables
#' @param pNA A logical value. If scale = TRUE, each column is transformed to have unit variance (default = TRUE).
#' @param nAllRespondants If the typeNA "ponc" is chosen,nAllRespondant corresponds to the number of complete individuals required
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
createNA=function(blocks,typeNA="block",pNA=0.1,nAllRespondants=4,output="list",seed=NULL)
{
  blocks= listOfMatrices(blocks)
  if(length(pNA)==1){pNA=rep(pNA,length(blocks))}
  if(length(pNA)!=1 & length(pNA)!=length(blocks)){stop_rgcca("pNA should be a number between 0 and 1 or a vector of the same size as blocks ")}
   match.arg(typeNA,c("block","rand","ponc","byvar"))
   match.arg(output,c("list","superblock"))
   if(!is.null(seed))
   {
       check_integer("seed",seed)
   }
   check_integer(nAllRespondants)
  if(is.null(rownames(blocks[[1]]))){rownames(blocks[[1]])=paste("S",1:dim(blocks[[1]])[1],sep="")}
    if(is.list(pNA)){message("The percentage of missing data is chosen by variable. 'byvar' typeNA is chosen")}
  if(typeNA=="block")
	{
	  # For all the one but the last one, 
	  nbloc=length(blocks)
	  indToRemove=list()
		for(i in 1:(length(blocks)-1))
		{
			n=nrow(blocks[[i]])
			p=ncol(blocks[[i]])
			nbNAparBloc=round(pNA[i]*n)	
			if(nbNAparBloc!=0)
			{
			 if(!is.null(seed)){set.seed(seed+i)}
		  	indToRemove[[i]]=sample(1:n, nbNAparBloc, replace = FALSE)	
		  	 blocks[[i]][indToRemove[[i]],]=NA
		    }
			else
			{
				indToRemove[[i]]=NA
			}
		}
	   n=nrow(blocks[[nbloc]])
	  p=ncol(blocks[[nbloc]])
	  nbNAparBloc=round(pNA[length(blocks)]*n)	
	  # The last block can not be a missing line if all other blocks are missing
	  #---------------------------------------------------------------------------
	  W=do.call(cbind,blocks[1:(nbloc-1)])
    indToTake=which(apply(W,1,function(x){sum(!is.na(x))})==0) #index with missing data for all first blocks
    if(length(indToTake)==0){indicesToUse=1:n}else{indicesToUse=(1:n)[!(1:n)%in%indToTake]}
	  if(nbNAparBloc!=0)
	  {
	      if(!is.null(seed)){set.seed(seed+length(blocks))}
	    indToRemove[[nbloc]]=sample(indicesToUse, min(length(indicesToUse),nbNAparBloc), replace = FALSE)	
	    blocks[[nbloc]][indToRemove[[nbloc]],]=NA
	   }
	  else
	  {
	    indToRemove[[nbloc]]=NA
	  }
    W2=do.call(cbind,blocks)
    subjectKept=rownames(blocks[[1]])[which(apply(W2,1,function(x){sum(is.na(x))})==0)]
    
  }
   
	if(typeNA=="rand")
	{
		
		for(i in 1:length(blocks))
		{
			n=nrow(blocks[[i]])
			p=ncol(blocks[[i]])
			nbNA=round(pNA[i]*n*p)
			if(nbNA!=0)
			{
				listeIndicePossible=merge(1:n,1:p)
				if(!is.null(seed)){set.seed(seed+i)}
				indToRemove=sample(1:(n*p), nbNA, replace = FALSE)
				listeIndicesToRemove=listeIndicePossible[indToRemove,]		
				for(j in 1:length(indToRemove)){blocks[[i]][listeIndicesToRemove[j,"x"],listeIndicesToRemove[j,"y"]]=NA}
			}
			else{indToRemove=NA}
		}
	}
	
	if(typeNA=="ponc") # valeurs ponctuelles en conservant un nombre d'individu complets
	{
	  n=nrow(blocks[[1]])
		allResp=sample(1:n,nAllRespondants)
		restData=(1:n)[!(1:n)%in%allResp]
	
		for(i in 1:length(blocks))
		{
			p=ncol(blocks[[i]])
			nbNA=round(pNA[i]*n*p)
			if(nbNA!=0)
			{
				listeIndicePossible=merge(restData,1:p)
				if(!is.null(seed)){set.seed(seed+i)}
				indToRemove=sample(1:dim(listeIndicePossible)[1], nbNA, replace = FALSE)
				listeIndicesToRemove=listeIndicePossible[indToRemove,]		
				for(j in 1:length(indToRemove)){blocks[[i]][listeIndicesToRemove[j,"x"],listeIndicesToRemove[j,"y"]]=NA}
			}
			else{indToRemove=NA}
		}
		W2=do.call(cbind,blocks)
		subjectKept=rownames(blocks[[1]])[which(apply(W2,1,function(x){sum(is.na(x))})==0)]
	}
    if(typeNA=="byvar")
    {
        
        indToRemove=list()
        if(length(pNA)!=length(blocks)){stop_rgcca("If pNA is a list, it should have the same length as blocks")}
        for(i in 1:length(blocks))
        {
            n=nrow(blocks[[i]])
            
            indToRemove[[i]]=list()
            if(length(pNA[[i]])!=dim(blocks[[i]])[2]){stop_rgcca(paste("pNA[",i,"] should have the same length that the number of variable in blocks[",i,"]"))}
            for(j in 1:dim(blocks[[i]])[2])
            {
               
                nbNA=round(pNA[[i]][j]*n)
                indToRemove[[i]][[j]]=sample(1:n,nbNA)
             
                blocks[[i]][indToRemove[[i]][[j]],j]=NA
            }
        }
        W2=do.call(cbind,blocks)
        subjectKept=rownames(blocks[[1]])[which(apply(W2,1,function(x){sum(is.na(x))})==0)]
    }
		
  if(output=="list"){res=list(dat=blocks,naPos=indToRemove,subjectKept=subjectKept)}else{res=do.call(cbind,blocks)}
	return(res)
}
