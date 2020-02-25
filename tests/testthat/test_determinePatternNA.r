#'# get_patternNA

#'''
X1=matrix(rnorm(150),30,5)
X2=matrix(rnorm(150),30,5)
X3=matrix(rnorm(150),30,5)
X1[1:10,1]=NA
X1[2,2]=NA
X1[11,]=NA
X2[1:10,1]=NA
colnames(X1)=paste0("A",1:5)
colnames(X2)=paste0("B",1:5)
colnames(X3)=paste0("C",1:5)
A=list(bloc1=X1,bloc2=X2,bloc3=X3)
res=get_patternNA(A)
# Check old functionalities
test_that("test_get_patternNA_1",{expect_true(sum(res$pctNA$bloc3)==0)})
test_that("test_get_patternNA_2",{expect_true(length(res$pctNAbyBlock)==3)})
test_that("test_get_patternNA_3",{expect_true(dim(res$completeSubjectByBlock)[2]==3)})

# setwd("/home/caroline.peltier/Bureau/RGCCA/R")
# lapply(list.files(),source)
# setwd("/home/caroline.peltier/Bureau/DATA/blocks_full")
# refData=readDataset(c("clinic_full","lipidomic","metabolomic","transcriptomic"),type="tsv",sep="\t",header=TRUE)
# get_patternNA(refData)
