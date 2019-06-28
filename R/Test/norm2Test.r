'# Test svdinit

'''
library(devtools)
set.seed(42)
X=matrix(rnorm(15),3,5)
# saving  the reference output
library(RGCCA)
T1=Sys.time()
resInitsvdRgcca=RGCCA:::initsvd(X)
T2=Sys.time()
T2-T1
save(resRgccaInitsvd,file="resRgccaInitsvd")

# loading the reference output
load(file="resRgccaInitsvd")

# loading the new function output
source("initsvd2.r")
resInitsvdTest=initsvd2(X)

# testing the same results than previously
all.equal(resInitsvdRgcca,resInitsvdTest)
        
# testing the additional features
X2=matrix(rnorm(15),3,5)
X2[2,]=NA
initsvd2(X2)
