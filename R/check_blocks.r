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
# check_blocks(A)
# A[[1]][2, 3] <- "character"
# check_blocks(A)
# A[[1]][2, 3] <- runif(1)
# init : boolean (FALSE by default) for the first block checking

check_blocks <- function(blocks, init = FALSE, add_NAlines=FALSE) {
    
    msg <- "In blocks arg:"
    
    if (!is.list(blocks))
        stop(paste(msg, "is not a list."))
    
    if (!init && length(blocks) < 2)
        stop(paste(msg, "should at least have two elements."))
    
    if (is.null(names(blocks)))
    {
        #        stop(paste(msg, "elements of the list should have names."))
        names(blocks)=paste0("block",1:length(blocks))
        warning("Warnings in check_blocks(A):\n blocks of the list had no names. The blocks were named block1,... blockJ in the order of the entered list")
    }
    # Gestion of the case of one variable only
    blocks=lapply(blocks,as.matrix)
    
    
    if (all(sapply(blocks, function(x) is.null(row.names(x)))))
    {
        if(sd(sapply(blocks,function(x)dim(x)[1]))==0)
        {
            blocks=lapply(blocks,function(x){rownames(x)=paste0("S",1:(dim(x)[1]));return(x)})
            print("Warnings in check_blocks(A):\n Elements of the list have no rownames. They were named as S1,...Sn") #TODO : verify
        }
        else
        {
            stop(paste(msg, "elements of the list should have rownames."))  
        }
       
    }
    if (any(sapply(blocks, function(x) is.null(colnames(x)))))
        stop(paste(msg, "elements of the list should have colnames."))
    
    inters_rows <- Reduce(intersect, lapply(blocks, row.names))
 
    if (length(inters_rows) == 0)
        warnings(paste(msg, "elements of the list should have at least a common rowname."))
    
    equal_rows <- Reduce(identical, lapply(blocks, row.names))
    
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
    
    if (any(sapply(blocks, is.character2)))
        stop(paste(msg, "an element contains non-numeric data."))
    
    for (i in seq(length(blocks)))
        if (is.character(blocks[[i]]))
            blocks[[i]] <- to_numeric(blocks[[i]])
   # Add lines if subjects are missing
    if(add_NAlines)
    {
        union_rows <- Reduce(union, lapply(blocks,row.names))
        blocks=lapply(blocks,function(x)
        {
            if(sum(!union_rows%in%rownames(x))!=0)
            {
                warnings("Some subjects are not present in some blocks. NA lines were added to have blocks with same dimensions") 
                y=matrix(NA,length(union_rows),dim(x)[2]);colnames(y)=colnames(x);rownames(y)=union_rows
                y[rownames(x),]=x
                return(y)
            }
            else
            {
                x=x[union_rows,]
            }
        })
    }
   
    # if any(sapply(blocks, is.character)) # optimization ?
    # if (any(is.na(unlist(blocks)))) {
    #     warning(paste(msg, "an element contains NA that will be imputed by mean."))
    #     for (i in seq(length(blocks)))
    #         blocks[[i]] <- impute_mean(blocks[[i]])
    # }
    
    # if (type != "pca")
    # if (length(blocks) > 1 && length(Reduce(intersect, lapply(blocks, colnames))))
    #     stop(paste(msg, "elements of the list should have different colnames."))
    # TODO: automatic conversation and warning
    
    invisible(blocks)
}