library(parallel)
library(FactoMineR)

# test on random data
#---------------------
set.seed(42);X1=matrix(rnorm(350),70,5);X2=matrix(rnorm(280),70,4)
A=list(X1,X2)
resultComparison=whichNAmethod(A,listMethods=c("nipals","mean"),patternNA=rep(0.1,2))

plotAnalysis(resultComparison,output="a",fileName=NULL,ylim=c(0,1),block="all",barType="sd",namePlot=NULL,width=480,height=480)

resultComparison=whichNAmethod(A,listMethods=c("complete","mean","nipals","knn","em"),patternNA=rep(0.3,2))
plotAnalysis(resultComparison,output="rv",fileName=NULL,ylim=c(0,1),block="all",barType="sd",namePlot=NULL,width=480,height=480)


plotAnalysis(resultComparison,output="a",fileName=NULL,ylim=c(0,1),block="all",barType="sd",namePlot=NULL,width=480,height=480)
plotAnalysis(resultComparison,output="bm",fileName=NULL,ylim=c(0,1),block="all",barType="sd",namePlot=NULL,width=480,height=480)


# Attention ! ! Ne fonctionne pas si trop de valeurs manquantes
resultComparison=whichNAmethod(A,listMethods=c("mean"),patternNA=rep(0.5,2))

# test on biosca
#---------------
setwd("/home/caroline.peltier/Bureau/EtudeNA/Datasets/Biosca/Reference")
refData=readDataset(c("CLI","MRS","VOL"))
resultComparison=whichNAmethod(refData,listMethods=c("complete","mean","nipals","knn","em"),patternNA=c(0.3,0.3,0.3))
plotAnalysis(resultComparison,output="rv",fileName=NULL,ylim=c(0,1),block="all",barType="sd",namePlot=NULL,width=480,height=480)
plotAnalysis(resultComparison,output="a",fileName=NULL,ylim=c(0.6,1),block="all",barType="sd",namePlot=NULL,width=480,height=480)




