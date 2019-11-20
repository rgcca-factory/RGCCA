#' The function defl.select() computes residual matrices \eqn{X_{1,h+1}, \ldots, X_{J,h+1}}. These 
#' residual matrices are determined according to the following formula: \eqn{X_{j,h+1} = X_{jh} - y_{jh} p_{jh}^t}.
#' @param yy  A matrix that contains the SGCCA block components of each block: \eqn{y_{1h}, \ldots, y_{Jh}}
#' @param rr  A list that contains the residual matrices \eqn{X_{1h}, \ldots, X_{Jh}} 
#' @param nncomp A \eqn{1 \times J} vector that contains the number of components to compute for each block.
#' @param nn  A \eqn{1 \times J} vector that contains the numbers of deflations for each block
#' @param nbloc Number of blocks.
#' @return \item{resdefl}{A list of \eqn{J} elements that contains \eqn{X_{1,h+1}, \ldots, X_{J,h+1}}.}
#' @return \item{pdefl}{A list of \eqn{J} elements that contains \eqn{p_{1h}, \ldots, p_{Jh}}.}
#' @title deflation function
#' @examples
#' 
#' set.seed(42);X1=matrix(rnorm(15),3,5);
#'set.seed(22);X2=matrix(rnorm(12),3,4);
#'set.seed(2);X3=matrix(rnorm(12),3,7);
#'A=list(X1,X2,X3)
#'yy=cbind(c(1,0,0),c(0,0,1),c(1/sqrt(2),1/sqrt(2),0)) # projection vectors are identity
#'res=defl.select (yy=yy, rr=A, nncomp=c(1,1,1), nn=1, nbloc=3)  
#' 
#' @export defl.select
defl.select=function (yy, rr, nncomp, nn, nbloc) 
{
  resdefl <- NULL
  pdefl <- list()
  for (q in 1:nbloc) {
    if (nn <= nncomp[q]) {
      defltmp <- deflation(as.matrix(rr[[q]]), as.matrix(yy[, q]))
      resdefl[[q]] <- defltmp$R
      pdefl[[q]]=as.matrix(defltmp$p)
    }
    else {
      resdefl[[q]] <- rr[[q]]
      pdefl[[q]] <- rep(0, NCOL(rr[[q]]))
    }
  }
  return(list(resdefl = resdefl, pdefl = pdefl))
}