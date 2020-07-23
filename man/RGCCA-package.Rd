\name{RGCCA-package}
\alias{RGCCA-package}
\alias{RGCCA}
\docType{package}
\title{
\packageTitle{RGCCA}
}
\description{
\packageDescription{RGCCA}
}

\author{
\packageAuthor{RGCCA}

Maintainer: \packageMaintainer{RGCCA}
}

\keyword{ package }
\seealso{\link{rgcca},\link{bootstrap},\link{rgcca_predict},\link{rgcca_crossvalidation},\link{rgcca_permutation},\link{get_patternNA}, \link{whichNAmethod},\link{naEvolution}
}
\examples{
data(Russett)
blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
    politic = Russett[, 6:11] )
res_rgcca=rgcca(blocks, response=3)

t0=Sys.time()
res_rgcca=rgcca(blocks, tau=0)
print(Sys.time()-t0)

cv=rgcca_cv(blocks,response=3,n_cores=1)
plot(cv)
b=bootstrap(res_rgcca, n_cores=1,nboot=10)
plot(b, n_cores=1)

res_rgcca=rgcca(blocks, response=3)

perm=rgcca_permutation(blocks, n_cores=1,nperm=10)
b=bootstrap(res_rgcca, n_cores=1,n_boot=10)

blockNA=blocks
blocks[[1]][1,2]=NA
get_patternNA(blocks)
#whichNAmethod(blocks)
naEvolution(blocks)


}

