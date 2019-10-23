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
#' @return \item{A}{A list of the imputed blocks}
#' @title imputeNN: Impute with k-Nearest Neighbors
#' @examples 
#'  set.seed(42);X1=matrix(rnorm(500),100,5);
#'  set.seed(22);X2=matrix(rnorm(400),100,4);
#'  set.seed(2);X3=matrix(rnorm(700),100,7);
#'  X1[1,]=NA;X2[7,1]=NA;X2[5,1]=NA;X3[3,]=NA;X3[4,]=NA
#'  rownames(X1)=rownames(X2)=rownames(X3)=paste("S",1:100,sep="")
#'  colnames(X1)=paste("A",1:5,sep="");colnames(X2)=paste("A",6:9,sep="");
#'  colnames(X3)=paste("A",10:16,sep="")
#' A=list(X1,X2,X3)  
#' Ares=imputeNN(A,k=1,output="mean", klim=NULL,scale=TRUE,sameBlockWeight=TRUE)


imputeNN <- function(
  A,
  k = 1,
  output = "mean",
  klim = NULL,
  scale = TRUE,
  sameBlockWeight = TRUE) {
  
  # TODO: no case for non euclidian
  distance <- function(vec, mat, method = "euclidean") {
    if (method == "euclidean") {
      res <- apply(mat, 1, function(v) {
        return(sqrt(sum((v - vec)^2)))
      })
    }
    return(res)
  }
  
  # Each variable is divided by its standard deviations, if scale =TRUE
  superblockNA <- do.call(cbind, A)
  superblockNAs <- scale3(superblockNA, scale = scale)
  
  # Each block is divided by its standard deviations, if sameBlockWeight=TRUE
  if (sameBlockWeight) {
    
    group <- unlist(lapply(A, "NCOL"))
    D <- matrix(0, dim(superblockNA)[1], dim(superblockNA)[2])
    debutBlock <- c(1, 1 + cumsum(group)[1:(length(group) - 1)])
    finBlock <- cumsum(group)
    
    for (u in 1:length(finBlock)) {
      variances <- apply(superblockNAs[, debutBlock[u]:finBlock[u]], 2, var, na.rm = T)
      D[, debutBlock[u]:finBlock[u]] <- 1/sqrt(sum(variances))
    }
    
  } else
    D <- matrix(1, dim(superblockNA)[1], dim(superblockNA)[2])
  
  superblockNAs2 <- superblockNAs * D
  # dectection of subjects with missing data
  namesInd <- rownames(A[[1]])
  posNA <- which(is.na(superblockNAs2), arr.ind = T)
  nbNA <- dim(posNA)[1]
  nbLineNA <- namesInd[unique(as.vector(posNA[, "row"]))]
  J <- length(A)
  
  # for each subject with missing values
  for (i in nbLineNA) {
    
    distances <- NULL
    # select the not-NA variables for this subject
    colForComparison <- colnames(superblockNAs2)[!is.na(superblockNAs2[as.character(i), ])]
    # select all the complete subjects
    linesWithoutNaInCFC <- (apply(superblockNAs2, 1, sum))
    rowForComparison <- names(linesWithoutNaInCFC[!is.na(linesWithoutNaInCFC)])
    
    if (length(rowForComparison) <= 5)
      stop("Not enough subjects with complete data (<=5)")
    if (k == "all") 
      k <- length(rowForComparison)
    
    if (length(rowForComparison) > 5) {
      
      # calculation of distances between the subject and all the complete subjects
      distances <- distance(
        superblockNAs2[as.character(i), 
                       colForComparison], 
        superblockNAs2[rowForComparison, colForComparison]
      )
      
      # group=cutree(hclust(dist(distances)),2)
      if (k == 1) {
        nameMin <- names(distances)[which.min(distances)]
        for (j in 1:J) {
          A[[j]][as.character(i), is.na(A[[j]][as.character(i), ])] <- A[[j]][nameMin, 
                                                                              is.na(A[[j]][as.character(i), ])]
        }
      }
      
      if (k == "all")
        k <- length(rowForComparison)
      if (k == "auto" || (k > 1 & k <= length(rowForComparison))) {
        # sorting the distances
        orderedDistance <- sort(distances)
        
        if (k == "auto")
          knb <- kChoice(orderedDistance, klim = klim)
        else
          knb <- k
        
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
          naCol <- which(is.na(A[[j]][as.character(i), ]))
          nameContributors <- as.character(contributors)
          
          if (length(naCol) > 0) {
            mat <- A[[j]][nameContributors, naCol]
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
            
            A[[j]][as.character(i), naCol] <- lineToImpute
          }
          
        }
        
      }
    }
  }
  return(A)
}

# bloc1=scale(matrix(rbind(c(1,1,2),c(1,1.1,2),c(3,5,1),c(2,2,5)),4,3),scale=FALSE);rownames(bloc1)=c('21','1','4','6');colnames(bloc1)=c('A','B','C')
# bloc2=scale(matrix(rbind(c(0,3,0),c(1,12,0),c(NA,NA,NA),c(1,1,2)),4,3),scale=FALSE);rownames(bloc2)=c('21','1','4','6');colnames(bloc2)=c('D','E','F')
# bloc3=scale(matrix(rbind(c(0,3,1),c(1,NA,0),c(1,1,1),c(1,2,2)),4,3),scale=FALSE);rownames(bloc3)=c('21','1','4','6');colnames(bloc3)=c('G','I','H')
# listRes=list(bloc1,bloc2,bloc3) knn2(listRes,k=2) listRes
# superblockNA=do.call(cbind,listRes) group=unlist(lapply(listRes,'NCOL'))
# D=matrix(0,dim(superblockNA)[1],dim(superblockNA)[2])
# debutBlock=c(1,1+cumsum(group)[1:(length(group)-1)]) finBlock=cumsum(group)
# for(k in 1:length(finBlock)) { D[,debutBlock[k]:finBlock[k]]=1/sqrt(group[k]) }
# # on divise chaque colonne par son ecart type superblockNAs=scale(superblockNA)
# # on divise par racine de p superblockNAs2=superblockNAs*D
# res=kNN(superblockNAs2)
