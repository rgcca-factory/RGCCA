#'Impute with k-Nearest Neighbors
#'
#'This method is used for the implementation of EM algorithm for missing data
#'
#' @param A A list of J blocks
#' @param k number of nearest neighbors. Can also be "all" or "auto"
#' @param output "mean","random" or "weightedMean". Corresponds to the kind of output required by the user. If "random" is chosen, the imputation will be done by selecting one neighbor among the k nearests. If "mean" is chosen, the imputation will be done by averaging all k-neighbors scores. If "weightedMean" is chosen, this average is weighted by the inverse of the distance.
#' @param klim Vector of two integers with klim(1)<klim(2). if k=auto, it is optimised between klim(1) and klim(2)
#' @param  scale  If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
#' @param sameBlockWeight A logical value indicating if the different blocks should have the same weight in the analysis (default, sameBlockWeight=TRUE)
#' @param superblock if TRUE the distance between two subjects is calculated on the superblock. If FALSE the distance is calculated by blocks. 
#' @return \item{A}{A list of the imputed blocks}
#' @title imputeNN: Impute with k-Nearest Neighbors

imputeNN <- function(
  A,
  k = 1,
  output = "mean",
  klim = NULL,
  scale = TRUE,
  sameBlockWeight = TRUE,
  superblock=TRUE
  ) {

    if(!superblock)
    {
        B=lapply(A,function(x){imputeNN(
            x,
            k = k,
            output = output,
            klim = klim,
            scale = scale,
            sameBlockWeight = FALSE,
            superblock=TRUE
        )})
        return(B)
    }
  # TODO: no case for non euclidian
  distance <- function(vec, mat, method = "euclidean") {
    if (method == "euclidean") {
      res <- apply(mat, 1, function(v) {
        return(sqrt(sum((v - vec)^2)))
      })
    }
    return(res)
  }
  
  # Each variable is centered in anycase, divided by its standard deviations, if scale =TRUE
 
     if(!is.matrix(A))
    { superblockNA <- do.call(cbind, A)}
     else{ superblockNA=A}

   superblockNAs <- scale3(superblockNA, scale = scale)
  
  # Each block is divided by its standard deviations, if sameBlockWeight=TRUE

  if (sameBlockWeight&!is.matrix(A)) {
 
    group <- unlist(lapply(A, "NCOL"))
    D <- matrix(0, dim(superblockNA)[1], dim(superblockNA)[2])
    debutBlock <- c(1, 1 + cumsum(group)[1:(length(group) - 1)])
    finBlock <- cumsum(group)
    
    for (u in 1:length(finBlock)) {
        
      if(finBlock[u]-debutBlock[u]!=0)
      {
          variances <- apply(superblockNAs[, debutBlock[u]:finBlock[u]], 2, var, na.rm = T)
      }
      else
      {
          variances=var(superblockNAs[, debutBlock[u]:finBlock[u]], na.rm = T)
      }
           D[, debutBlock[u]:finBlock[u]] <- 1/sqrt(sum(variances))
    }
    
  } else
    D <- matrix(1, dim(superblockNA)[1], dim(superblockNA)[2])
  
  superblockNAs2 <- superblockNAs * D
 
   # dectection of subjects with missing data
  if(is.matrix(A)){  namesInd <- rownames(A)}else{  namesInd <- rownames(A[[1]])}

  posNA <- which(is.na(superblockNAs2), arr.ind = T)
  nbNA <- dim(posNA)[1]
  nbLineNA <- namesInd[unique(as.vector(posNA[, "row"]))]
  if(is.matrix(A)){J=1}else{  J <- length(A)}


  # for each subject with missing values
  for (i in nbLineNA) {
    
    distances <- NULL
    # select the not-NA variables for this subject
    colForComparison <- colnames(superblockNAs2)[!is.na(superblockNAs2[as.character(i), ])]
    # select all the complete subjects
    linesWithoutNaInCFC <- (apply(superblockNAs2, 1, sum))
    rowForComparison <- names(linesWithoutNaInCFC[!is.na(linesWithoutNaInCFC)])
    rowForComparison=rowForComparison[rowForComparison!=i]
    
    if (length(rowForComparison) <= 4)
      stop("Not enough subjects with complete data (<=4)")
    if (k == "all") 
      k=length(rowForComparison)
    if(length(colForComparison)==0){stop("One subject is missing for all blocks. Please change your imputation function, or choose superblock=TRUE")}
    if (length(rowForComparison) > 4)
    {
       
          # calculation of distances between the subject and all the complete subjects
    # if(superblock)
   #  {
         distances <- distance(
             superblockNAs2[as.character(i), 
                            colForComparison], 
             superblockNAs2[rowForComparison, colForComparison]
         )

   #  }
    # else
    # {
    #     distances=NULL
    #     for(j in 1:J)
    #     {
    #         if(scale){blockNAs=scale(blockNAs,scale=scale)}
    #         colForComparisonj <- colnames(blockNAs)[!is.na(blockNAs[as.character(i), ])]
    #         # select all the complete subjects
    #         linesWithoutNaInCFCj <- (apply(blockNAs, 1, sum))
    #         rowForComparisonj <- names(linesWithoutNaInCFCj[!is.na(linesWithoutNaInCFCj)])
    #         rowForComparisonj=rowForComparisonj[rowForComparisonj!=i]
    #         distances[[j]] <- distance(
    #         blockNAs[as.character(i), 
    #                        colForComparisonj], 
    #        blockNAs[rowForComparisonj, colForComparisonj])
    #     }
    # }
    # 
      
      # Imputation if k=1
      if (k == 1) {
        nameMin <- names(distances)[which.min(distances)]
        for (j in 1:J) {
            if(is.matrix(A))
            {
                A[as.character(i), is.na(A[as.character(i), ])] <- A[nameMin, 
                                                                                    is.na(A[as.character(i), ])] 
            }
            else
            {
                A[[j]][as.character(i), is.na(A[[j]][as.character(i), ])] <- A[[j]][nameMin, 
                                                                                    is.na(A[[j]][as.character(i), ])]
            }
          
        }
      }

      if (k == "all")
        k <- length(rowForComparison)
      
      if (k == "auto" || (k > 1 & k <= length(rowForComparison))) 
      {
        # sorting the distances
        orderedDistance <- sort(distances)
       # /!\ To be uncommented 
        if (k == "auto")
          knb <- kChoice(orderedDistance, klim = klim)
        else
          knb <- k
          # /!\end  To be uncommented 
        contributors <- names(orderedDistance)[1:knb]
        
        # calculating a vector of weights for each subject: each weight component
        # correspond to the inverse of proportion of distances between the subject and
        # one complete subject
        
        if (output == "weightedMean")
          w <- (1/(orderedDistance/sum(orderedDistance)))[1:knb]
        if (output == "random")
          randomIndex <- contributors[sample(knb, 1)]
        
        # imputation for each block
       
        for (j in 1:J) {
          if(J>1){naCol <- which(is.na(A[[j]][as.character(i), ]))}
          if(J==1){naCol <- which(is.na(A[as.character(i), ]))}
          nameContributors <- as.character(contributors)
          
          if (length(naCol) > 0) {
            if(J>1){ mat <- A[[j]][nameContributors, naCol]}
            if(J==1){ mat <- A[nameContributors, naCol]}
           
            if (length(naCol) == 1) {
              if (output == "weightedMean") {
                names(mat) <- nameContributors
                lineToImpute <- weighted.mean(mat, w)
              }
              if (output == "mean") {
                names(mat) <- nameContributors
                lineToImpute <- mean(mat)
              }
              if (output == "random") {
                names(mat) <- nameContributors
                lineToImpute <- mat[randomIndex]
              }
            }
            if (length(naCol) > 1) {
              if (output == "weightedMean") {
                rownames(mat) <- nameContributors
                lineToImpute <- apply(mat, 2, "weighted.mean", w)
              }
              if (output == "mean") {
                rownames(mat) <- nameContributors
                lineToImpute <- apply(mat, 2, "mean")
              }
              if (output == "random") {
                rownames(mat) <- nameContributors
                lineToImpute <- mat[randomIndex, ]
              }
            }
            if(J==1){ A[as.character(i), naCol] <- lineToImpute}
            if(J>1){ A[[j]][as.character(i), naCol] <- lineToImpute}
              
            
          }
          
        }
        
      }
    }
  }
  return(A)
}

