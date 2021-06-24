# Internal function which does not have to be used by the users
# @author Arnaud Gloaguen
# @keywords internal
BinarySearch <-
function(argu,sumabs){
  if(norm2(argu)==0 || sum(abs(argu/norm2(argu)))<=sumabs) return(0)
  lam_max = max(abs(argu))
  lam1      <- 0
  lam2      <- lam_max
  iter <- 1
  while(iter < 500){
    su <- soft(argu,(lam1+lam2)/2)
    if(sum(abs(su/norm2(su)))<sumabs){
      lam2 <- (lam1+lam2)/2
    } else {
      lam1 <- (lam1+lam2)/2
    }
    if((lam2-lam1)/lam1 < 1e-10){
      if (lam2 != lam_max){
        return(lam2)
      }else{
        return(lam1)
      }
    }
    iter <- iter+1
  }
  warning("Didn't quite converge")
  return((lam1+lam2)/2)
}

BinarySearch2 <- function(argu, sumabs,coefs) {
  normalize = function(X) {
    apply(X, 2, function(x) {
      mu = norm(x, type = "2")
      if (mu == 0) return(x)
      x / mu
    })
  }
  if (any(apply(argu, 2, norm, type = "2") == 0) || sum(abs(normalize(argu) %*% coefs)) <= sumabs) return(0)
  lam_max = max(abs(argu))
  lam1      <- 0
  lam2      <- lam_max
  iter <- 1
  while (iter < 500) {
    su <- soft(argu, (lam1 + lam2) / 2)
    if (sum(abs(normalize(su) %*% coefs)) < sumabs) {
      lam2 <- (lam1 + lam2)/2
    } else {
      lam1 <- (lam1 + lam2)/2
    }
    if ((lam2 - lam1)/lam1 < 1e-15) {
      if (lam2 != lam_max) {
        return(lam2)
      }else{
        return(lam1)
      }
    }
    iter <- iter + 1
  }
  stop("Didn't quite converge")
  return((lam1 + lam2)/2)
}

BinarySearch3 <- function(argu, sumabs, coefs) {
  if (norm2(argu) == 0 || sum(abs(coefs * argu/(norm2(argu) * sqrt(length(argu))))) <= sumabs) return(0)
  lam_max = max(abs(argu))
  lam1      <- 0
  lam2      <- lam_max
  iter <- 1
  while (iter < 500) {
    su <- soft(argu, (lam1 + lam2) / 2)
    if (sum(abs(coefs * su / (norm2(su) * sqrt(length(argu))))) < sumabs) {
      lam2 <- (lam1 + lam2)/2
    } else {
      lam1 <- (lam1 + lam2)/2
    }
    if ((lam2 - lam1)/lam1 < 1e-10) {
      if (lam2 != lam_max) {
        return(lam2)
      }else{
        return(lam1)
      }
    }
    iter <- iter + 1
  }
  warning("Didn't quite converge")
  return((lam1 + lam2)/2)
}
