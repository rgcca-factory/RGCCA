set.seed(42);X1=matrix(rnorm(350),70,5);X2=matrix(rnorm(280),70,4)
A=list(X1,X2)
listResults=naEvolution(refData=A,prctNA=c(0.1,0.2,0.3,0.4),listMethods=c("mean","complete","nipals","em","sem1","knn4","complete"))
plotEvol(listFinale=listResults,ylim=c(0,1),output="a")

#--------------------
# test on russets
#---------------------
X_agric =as.matrix(Russett[,c("gini","farm","rent")])
X_ind = as.matrix(Russett[,c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("demostab", "dictator")])
X_agric[c(2,4),]=NA
X_ind[1,]=NA
X_polit[5,1]=NA
A = list(agri=X_agric, ind=X_ind, polit=X_polit)
listResults=naEvolution(A=A,listMethods=c("complete","mean","em","superblockEM","knn4","nipals"),prctNA=c(0.1,0.2,0.3))
# bug pour plus grand que 04
plotEvol(listResults,ylim=c(0.9,1),output="a")
plotEvol(listResults,ylim=NULL,output="rv")
plotEvol(listResults,ylim=NULL,output="rvComplete")
plotEvol(listResults,ylim=NULL,output="bm")
#------------------
# test on biosca
#------------------
setwd("/home/caroline.peltier/Bureau/EtudeNA/Datasets/Biosca/Reference")
refData=readDataset(c("CLI","MRS","VOL"))
listResults=naEvolution(A=refData,listMethods=c("complete","mean","em","sem1","knn4","nipals"),prctNA=c(0.1,0.2,0.3),typeNA="ponc")
A=refData
plotEvol(listResults,ylim=NULL,output="a")
plotEvol(listResults,ylim=NULL,output="rv")
plotEvol(listResults,ylim=NULL,output="rvComplete")
plotEvol(listResults,ylim=NULL,output="bm")
plotEvol(listResults,ylim=c(0.4,1.5),output="rmse")
