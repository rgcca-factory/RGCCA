# # setwd("/home/caroline.peltier/Bureau/RGCCA")
set.seed(42);X1=matrix(rnorm(350),70,5);X2=matrix(rnorm(280),70,4);
colnames(X1)=c("A","B","C","D","E")
colnames(X2)=c("a","b","c","d")
rownames(X1)=rownames(X2)=paste("S",1:70)
 A=list(X1,X2)
 listResults=naEvolution(blocks=A,prctNA=c(0.1,0.2),listMethods=c("mean","complete","nipals","em","sem1","knn4"),nDatasets = 1)
 
 
 listResults=naEvolution(blocks=A,prctNA=c(0.1,0.2),listMethods=c("knn4","knn5"),nDatasets = 1)
 
 res=1
 #graphics.off()
 res=plot.naEvolution(x=listResults,ylim=c(0,1),output="a")
 test_that("naEvolution_1",{expect_true(is.null(res))})

