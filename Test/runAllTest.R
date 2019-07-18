# setwd("/home/caroline.peltier/Bureau/RGCCA")
rm(list=ls())
library(RGCCA)
library(MASS)
library(nipals)
namesFiles=dir("./R")
# loading functions in R directory
sapply(namesFiles,function(x){source(paste0("./R/",x))})

# loading test functions in Test directory
# nameTest=dir("./Test")
# nameTest=nameTest[-which(nameTest=="runAllTest.R")]
nameTest=c("cov3Test.R","createNATest.R","defl.selectTest.r","imputeColmeansTest.R","imputeNNTest.R","imputeSBTest.R","initsvdTest.r","intersectionTest.r","MIRGCCATest.R","pmTest.r","rgccakTest.r","rgccaTest.r","scale3Test.R")
for(i in 1:length(nameTest))
{
  print(nameTest[i])
  source(paste0("./Test/",nameTest[[i]]),echo=TRUE)
}
