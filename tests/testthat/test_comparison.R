# '# ImputeSB test
# 
# '''
#setwd("/home/caroline.peltier/Bureau/RGCCA")
# library(RGCCA)
# library(MASS)
# data(Russett)
# X_agric =as.matrix(Russett[,c("gini","farm","rent")])
# X_ind = as.matrix(Russett[,c("gnpr","labo")])
# X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
# A = list(agri=X_agric, ind=X_ind, polit=X_polit)
# 
# refRgcca=rgccad(A, connection = 1 - diag(length(A)), tau = rep(1, length(A)), ncomp = rep(2, length(A)), scheme = "centroid", scale = TRUE,   init = "svd", bias = TRUE, tol = 1e-08, verbose = TRUE,sameBlockWeight=TRUE,na.rm=FALSE)
# listNAdataset=  createNA(blocks=A,option="block",pNA=rep(0.2,3),nAllRespondants=10,output="list")
# selectCompletePatient=listNAdataset$subjectKept
# methodRgcca=rgccaNa(blocks=listNAdataset$dat,method="em",verbose=FALSE,ncomp=rep(2,length(A)))
# resComp=NULL
# resComp=comparison(rgcca1=refRgcca,rgcca2=methodRgcca$rgcca,selectPatient=selectCompletePatient,indNA=methodRgcca$indNA)
# test_that("comparison_1",{expect_true(!is.null(resComp))})