#' Test initsvd
#' 
#' 

wd="/home/caroline.peltier/Bureau/RGCCAtoPush"
setwd(wd)


# saving the reference output
setwd("./Test/Results")
library(RGCCA)
set.seed(42)
X=matrix(rnorm(15),3,5)
T1=Sys.time()
resInitsvdRgcca=RGCCA:::initsvd(X)
T2=Sys.time()
tInitsvdRgcca=T2-T1 #Time difference of 0.001314878 secs
save(resInitsvdRgcca,file="resInitsvdRgcca")
setwd("./../..")


# loading the reference output
setwd("./Test/Results")
load(file="resInitsvdRgcca")

# loading the new function output
setwd("./../..")
source("initsvd2.r")
T1=Sys.time()
resInitsvdTest=initsvd2(X)
T2=Sys.time()
tInitsvdTest=T2-T1
# testing the same results than previously
all.equal(resInitsvdRgcca,resInitsvdTest)
tInitsvdRgcca
tInitsvdTest

# testing the additional features
X2=matrix(rnorm(15),3,5)
X2[2,]=NA
initsvd2(X2)
