# # setwd("/home/caroline.peltier/Bureau/RGCCA")
# set.seed(42);X1=matrix(rnorm(350),70,5);X2=matrix(rnorm(280),70,4)
# A=list(X1,X2)
# listResults=naEvolution(A=A,prctNA=c(0.1,0.2),listMethods=c("mean","complete","nipals","em","sem1","knn4"))
# res=1
# #graphics.off()
# res=plotEvol(listFinale=listResults,ylim=c(0,1),output="a")
# test_that("naEvolution_1",{expect_true(is.null(res))})

