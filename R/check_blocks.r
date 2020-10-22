# rand_mat <- function(x) matrix(runif(9), 3, 3)
# A = lapply(1:3, rand_mat)
# check_blocks(A)
# names(A) <- LETTERS[1:3]
# check_blocks(A[1])
# check_blocks(A)
# row.names(A[[1]]) <- letters[1:3]
# check_blocks(A)
# for(i in 1:3)
#   row.names(A[[i]]) <- letters[(0+i*3):(2+i*3)]
# check_blocks(A)
# for(i in 1:3)
#   row.names(A[[i]]) <- letters[1:3]
# A[[1]][2, 3] <- NA
# for(i in 1:3)
#   colnames(A[[i]]) <- letters[(0+i*3):(2+i*3)]
# check_blocks(A,add_NAlines=TRUE)
# A[[1]][2, 3] <- "character"
# check_blocks(A)
# A[[1]][2, 3] <- runif(1)
# init : boolean (FALSE by default) for the first block checking

check_blocks <- function(blocks, init = FALSE, n = 2, add_NAlines=FALSE, allow_unnames =  TRUE,quiet=FALSE, no_character = FALSE) {
    
  
    msg <- ""
    if(is.matrix(blocks))
    {
        blocks=list(blocks)
    }
    if (!is.list(blocks))
        stop_rgcca(paste(msg, "is not a list."))
    
    if (!init && length(blocks) < n)
        stop_rgcca(paste(msg, "should at least have two elements."))
    
    # Completing block names
    if (is.null(names(blocks))){
        names(blocks)=paste0("block",1:length(blocks))
        message("Blocks of the list had no names. The blocks were named block1,... blockJ in the order of the entered list")
    }
    # Gestion of the case of one variable only
    blocks=lapply(blocks,as.matrix)
    nameBlocks=names(blocks)

    # Dealing with rownames (if they are all missing)
    if (all(sapply(blocks, function(x) is.null(row.names(x))))){
        if(sd(sapply(blocks,function(x)dim(x)[1]))==0 && allow_unnames){
            blocks=lapply(blocks,function(x){rownames(x)=paste0("S",1:(dim(x)[1]));return(x)})
            message("Elements of the list have no rownames. They were named as S1,...Sn \n ") #TODO : verify
        }
        else
            stop_rgcca(paste(msg, "elements of the list should have rownames.\n "))
    }
  
    
    # Dealing with colnames (in case of one variable onaly)
     blocks1=lapply(1:length(blocks),function(i)
    {
             if(dim(blocks[[i]])[2]==1&is.null(colnames(blocks[[i]])))
                 {
                     colnames(blocks[[i]])=nameBlocks[i]
                     return(blocks[[i]])
             }
         else{ return(blocks[[i]])}
         })
    blocks=blocks1
    names(blocks)=nameBlocks
    if (any(sapply(blocks, function(x) is.null(colnames(x)))))
    {
        stop_rgcca(paste(msg, "elements of the list should have colnames.\n "))
    }
    # if one of the colnames is identical in one block and another one
    if(sum(duplicated(unlist(sapply(blocks,colnames))))!=0)
    {
        if(!quiet)
        {
            message("At least one variable name is duplicated: the block names are added for avoiding confusion \n")
        }
       blocks_i=lapply(1:length(blocks),function(i){x=blocks[[i]];colnames(x)=paste(names(blocks)[i],colnames(blocks[[i]]),sep="_");return(x)})
        names(blocks_i)=names(blocks)
        blocks=blocks_i
    }
 
    # Dealing with rownames (if one of them is missing but the block sizes are the sames)
    if(any(sapply(blocks, function(x) is.null(row.names(x)))))
    {
        matrixOfRownames=Reduce(cbind,lapply(blocks,row.names))
        if(sum(!apply(matrixOfRownames,2,function(x) x==matrixOfRownames[,1]))==0)
        {
            blocks=lapply(blocks,function(x){row.names(x)=matrixOfRownames[,1];return(x)})
        }
    }
    
    lapply(blocks,function(x){
        resdup=duplicated(rownames(x));
        if(sum(resdup)!=0)
            {
                if(!quiet)
                {
                    warning(paste0("Rownames are duplicated and were removed : ", rownames(x)[resdup],"\n"))
                    
                }
         }
        }
        )
    inters_rows <- Reduce(intersect, lapply(blocks, row.names))
 #   if(length(inters_rows)<min(sapply(blocks,function(x){nrow(x)})))
    
    if (length(inters_rows) == 0)
        stop_rgcca(paste(msg, "elements of the list should have at least a common rowname.\n "))

    equal_rows <- Reduce(identical, lapply(blocks, row.names))

    # If add_NAlines=FALSE, taking the intersection_list
    if(!add_NAlines)
    {
       
        if (length(blocks) > 1 && !equal_rows)
            blocks <- common_rows(blocks)
    }

    if (init) {
   
        blocks <- remove_null_sd(blocks)
        for (i in seq(length(blocks)))
            attributes(blocks[[i]])$nrow <- nrow(blocks[[i]])
    }
   

    if (no_character) {
        if (any(sapply(blocks, is.character2)))
            stop(paste(msg, "an element contains non-numeric data."))
        
        for (i in seq(length(blocks)))
           if (is.character(blocks[[i]]))
               blocks[[i]] <- to_numeric(blocks[[i]])
    }
   # Add lines if subjects are missing

    if(add_NAlines)
    {
      
        union_rows <- Reduce(union, lapply(blocks,row.names))
        blocks2=lapply(nameBlocks,function(name)
        {
            if(sum(!union_rows%in%rownames(blocks[[name]]))!=0) # if some subjects are missing (in the rownames)
            {
                message("Some subjects are not present in some blocks. NA lines were added to have blocks with same dimensions")
                y=matrix(NA,length(union_rows),ncol=ifelse(is.null(dim(blocks[[name]])),1,dim(blocks[[name]])[2]));
                if(is.null(dim(blocks[[name]]))){colnames(y)=name}else{colnames(y)=colnames(blocks[[name]])};rownames(y)=union_rows
                y[rownames(blocks[[name]]),]=blocks[[name]]
                return(y)
            }
            else
            {
                if(dim( blocks[[name]])[2]==1) # if the block has 1 variable
                {
                   y=matrix( blocks[[name]][union_rows,],ncol=1);
                    rownames( y)=union_rows;
                    colnames(y)=name
                    
                 }
                else
                {
                   y=blocks[[name]][union_rows,]
                }
                return(y)
            }
        })

        names(blocks2)=nameBlocks
       blocks=blocks2 
       
    }
   
    # if any(sapply(blocks, is.character)) # optimization ?
    # if (any(is.na(unlist(blocks)))) {
    #     warning(paste(msg, "an element contains NA that will be imputed by mean."))
    #     for (i in seq(length(blocks)))
    #         blocks[[i]] <- impute_mean(blocks[[i]])
    # }
    
    # if (type != "pca")
    # if (length(blocks) > 1 && length(Reduce(intersect, lapply(blocks, colnames))))
    #     stop_rgcca(paste(msg, "elements of the list should have different colnames."))
    # TODO: automatic conversation and warning

    invisible(blocks)
}